from __future__ import annotations

import json
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from cmdstanpy import CmdStanModel

from Stan_CodeGenerator_extended import (
    DATASET_ORDER,
    DiagnosticSpec,
    diagnostic_specs,
    generate_from_config,
)


PROJECT_DIR = Path(__file__).resolve().parent

# Physics-model controls.
ANALYSIS_NAME = "CasePDG-noUB_SM_BGLnf1_extended"
FF_MODEL = "BGL"
FF_NF = 1
NP_MODEL = "SM"
HQET_MODEL = "3/2/1"
UB_MODEL = "no-UB"  # "no-UB", "soft-UB", or "hard-UB"
CUT_CONTROL = 1.0e-4
BELLE_FIT_OPTION = "sc"  # "sc", "re_small", "re_large", "de", "pdg", or "rms"

# Run controls.  The defaults perform only the all-data fit.  LODO, node
# splitting, and covariance scans can be enabled independently.
RUN_ALL_DATA = True
RUN_LODO = False
LODO_DATASETS = list(DATASET_ORDER)

RUN_NODE_SPLIT = False
NODE_SPLIT_DATASETS: list[str] = []

RUN_COVARIANCE_SCAN = False
COVARIANCE_SCAN_DATASETS: list[str] = []
COVARIANCE_ALPHAS = [0.0, 0.25, 0.5, 0.75, 1.0]

# Scalar or vector outputs used for LODO influence summaries.
PHYSICS_OUTPUTS = ["Vcb_mean"]

# Physical prediction vectors used as the default node-split separator.  The
# list can be shortened to a pre-declared common FF grid before a production
# analysis.
NODE_SPLIT_SEPARATOR_PREFIXES = [
    "MILC15func",
    "MILC21func",
    "JLQCD23func",
    "HPQCD23func",
    "LCSR18func",
    "LCSR23func",
]
NODE_SPLIT_RANK_TOL = 1.0e-10

INPUT_ROOT = PROJECT_DIR / "Stan_inputs"
OUTPUT_ROOT = PROJECT_DIR / "Stan_outputs"
FIT_FUNCTION_CACHE_ROOT = INPUT_ROOT / "StanFunctionCache"

UB_MAX = np.ones(4)
UB_SIGMA = np.full(4, 0.05)

PRE_FIT_SETTINGS = {
    "chains": 10,
    "iter_warmup": 500,
    "iter_sampling": 200,
    "adapt_delta": 0.95,
    "max_treedepth": 12,
}
PRE_FIT_CHAIN_TIMEOUT = 600.0
PRE_FIT_PARALLEL_CHAINS = None
PRE_FIT_MIN_SUCCESSFUL_CHAINS = 1
PRE_FIT_LP_THRESHOLD_DROP = 200.0
STAN_SIG_FIGS = 18

FIT_SETTINGS = {
    "chains": 3,
    "iter_warmup": 3000,
    "iter_sampling": 6000,
    "adapt_delta": 0.99,
    "max_treedepth": 15,
    "seed": 42,
}


@dataclass(frozen=True)
class FitCase:
    name: str
    kind: str
    selected_data: tuple[str, ...]
    target_dataset: str | None = None
    alpha: float | None = None


@dataclass
class FitResult:
    case: FitCase
    draws: pd.DataFrame
    stan_data: dict[str, Any]
    output_file: Path


def effective_ub_model() -> str:
    return UB_MODEL if FF_MODEL == "HQET" else "no-UB"


def canonical_dataset_list(values: list[str], label: str) -> list[str]:
    unknown = sorted(set(values) - set(DATASET_ORDER))
    if unknown:
        raise ValueError(f"Unknown {label} dataset(s): {', '.join(unknown)}")
    return [key for key in DATASET_ORDER if key in set(values)]


def build_generator_config() -> dict[str, Any]:
    return {
        "WorkingDirectory": str(PROJECT_DIR),
        "AnalysisName": ANALYSIS_NAME,
        "FFModel": FF_MODEL,
        "FFnf": FF_NF,
        "NPModel": NP_MODEL,
        "HQETModel": HQET_MODEL,
        "UBModel": effective_ub_model(),
        "CutControl": CUT_CONTROL,
        # The master source always contains all datasets.  This list only
        # supplies the default run-time mask in the generated JSON.
        "SelectedData": DATASET_ORDER,
        "GeneratedOutputs": ["VcbFromBr"],
        "OutputRoot": str(OUTPUT_ROOT),
        "InputDirectory": str(INPUT_ROOT),
        "DataPointFilePrefix": "standata",
        "BelleFitOption": BELLE_FIT_OPTION,
        "FitFunctionCacheDirectory": str(FIT_FUNCTION_CACHE_ROOT),
        "HPQCD23FitFile": str(
            PROJECT_DIR.parent
            / "mycode_FitAnalysisCode_v1.2"
            / "Stan_output"
            / "StanFit_HPQCD23_full_1001_1917.csv"
        ),
    }


def run_python_generator() -> dict[str, Any]:
    config = build_generator_config()
    config_path = INPUT_ROOT / f"config_{ANALYSIS_NAME}.json"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
    return generate_from_config(config)


def load_stan_data(data_file: Path) -> dict[str, Any]:
    with data_file.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)

    stan_data: dict[str, Any] = {
        key: np.asarray(value) if isinstance(value, list) else value
        for key, value in raw.items()
    }
    ub_model = effective_ub_model()
    if ub_model == "soft-UB":
        stan_data.update({"UBmax": UB_MAX, "UB_sigma": UB_SIGMA})
    elif ub_model == "hard-UB":
        stan_data.update({"UBmax": UB_MAX})
    return stan_data


def copy_stan_data(stan_data: dict[str, Any]) -> dict[str, Any]:
    copied: dict[str, Any] = {}
    for key, value in stan_data.items():
        copied[key] = value.copy() if isinstance(value, np.ndarray) else value
    return copied


def dataset_mask(selected_data: tuple[str, ...]) -> np.ndarray:
    selected = set(selected_data)
    return np.asarray([int(key in selected) for key in DATASET_ORDER], dtype=int)


def shrink_covariance(covariance: np.ndarray, alpha: float) -> np.ndarray:
    if not 0.0 <= alpha <= 1.0:
        raise ValueError(f"Covariance shrinkage alpha must be in [0, 1], got {alpha}.")
    sigma = np.sqrt(np.diag(covariance))
    if np.any(sigma <= 0):
        raise ValueError("Covariance shrinkage requires positive marginal variances.")
    correlation = covariance / np.outer(sigma, sigma)
    shrunk_correlation = (1.0 - alpha) * np.eye(len(sigma)) + alpha * correlation
    shrunk = np.outer(sigma, sigma) * shrunk_correlation
    return 0.5 * (shrunk + shrunk.T)


def prepare_case_data(base_data: dict[str, Any], case: FitCase) -> dict[str, Any]:
    stan_data = copy_stan_data(base_data)
    stan_data["use_dataset"] = dataset_mask(case.selected_data)
    if case.kind == "covariance_scan":
        if case.target_dataset is None or case.alpha is None:
            raise ValueError(f"Incomplete covariance-scan case: {case}")
        spec = diagnostic_specs(NP_MODEL)[case.target_dataset]
        nominal = np.asarray(base_data[spec.covariance_name], dtype=float)
        stan_data[spec.covariance_name] = shrink_covariance(nominal, case.alpha)
    return stan_data


def alpha_tag(alpha: float) -> str:
    return f"{alpha:.3f}".rstrip("0").rstrip(".").replace(".", "p")


def build_fit_cases() -> list[FitCase]:
    lodo_targets = canonical_dataset_list(LODO_DATASETS if RUN_LODO else [], "LODO")
    node_targets = canonical_dataset_list(
        NODE_SPLIT_DATASETS if RUN_NODE_SPLIT else [], "node-split"
    )
    scan_targets = canonical_dataset_list(
        COVARIANCE_SCAN_DATASETS if RUN_COVARIANCE_SCAN else [],
        "covariance-scan",
    )

    cases: list[FitCase] = []
    need_baseline = RUN_ALL_DATA or bool(lodo_targets or node_targets or scan_targets)
    if need_baseline:
        cases.append(FitCase("AllData", "all_data", tuple(DATASET_ORDER)))

    for target in canonical_dataset_list(
        list(set(lodo_targets) | set(node_targets) | set(scan_targets)),
        "LODO/node-split/covariance-scan",
    ):
        selected = tuple(key for key in DATASET_ORDER if key != target)
        cases.append(FitCase(f"LODO_{target}", "lodo", selected, target))

    for target in node_targets:
        cases.append(FitCase(f"Only_{target}", "dataset_only", (target,), target))

    for target in scan_targets:
        for alpha in COVARIANCE_ALPHAS:
            if np.isclose(alpha, 1.0):
                continue
            cases.append(
                FitCase(
                    f"CovShrink_{target}_a{alpha_tag(alpha)}",
                    "covariance_scan",
                    tuple(DATASET_ORDER),
                    target,
                    float(alpha),
                )
            )

    names = [case.name for case in cases]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate fit-case names were generated.")
    return cases


def good_chain_ids(draws: pd.DataFrame, threshold_drop: float = 200.0) -> list[int]:
    threshold = draws["lp__"].max() - threshold_drop
    chains: list[int] = []
    for chain in sorted(draws["chain__"].unique()):
        sub = draws[draws["chain__"] == chain]["lp__"]
        print(f"Chain {int(chain)}: sample_length = {len(sub)}, lp_max = {sub.max():.2f}")
        if sub.mean() > threshold:
            chains.append(int(chain))
    print(f"Good chains: {chains}")
    return chains


def run_pre_fit(
    model: CmdStanModel,
    stan_data: dict[str, Any],
    output_dir: Path,
) -> pd.DataFrame:
    settings = dict(PRE_FIT_SETTINGS)
    chains = int(settings.pop("chains", 1))
    timeout = settings.pop("timeout", PRE_FIT_CHAIN_TIMEOUT)
    parallel_chains = settings.pop("parallel_chains", PRE_FIT_PARALLEL_CHAINS)
    if parallel_chains is None:
        parallel_chains = min(chains, os.cpu_count() or chains)
    parallel_chains = max(1, min(chains, int(parallel_chains)))

    output_dir.mkdir(parents=True, exist_ok=True)
    print(
        "Pre-fit settings: "
        f"chains = {chains}, parallel_chains = {parallel_chains}, "
        f"timeout_per_chain = {timeout}"
    )

    def run_chain(chain_id: int) -> pd.DataFrame | None:
        chain_dir = output_dir / f"chain_{chain_id}"
        chain_dir.mkdir(parents=True, exist_ok=True)
        try:
            fit = model.sample(
                data=stan_data,
                chains=1,
                parallel_chains=1,
                chain_ids=[chain_id],
                output_dir=chain_dir,
                timeout=timeout,
                sig_figs=STAN_SIG_FIGS,
                **settings,
            )
        except Exception as exc:
            print(f"Pre-fit chain {chain_id} discarded: {type(exc).__name__}: {exc}")
            return None

        draws = fit.draws_pd()
        lp_max = draws["lp__"].max() if "lp__" in draws else float("nan")
        print(
            f"Pre-fit chain {chain_id} accepted: "
            f"sample_length = {len(draws)}, lp_max = {lp_max:.2f}"
        )
        return draws

    accepted: list[pd.DataFrame] = []
    with ThreadPoolExecutor(max_workers=parallel_chains) as executor:
        futures = {
            executor.submit(run_chain, chain_id): chain_id
            for chain_id in range(1, chains + 1)
        }
        for future in as_completed(futures):
            draws = future.result()
            if draws is not None and not draws.empty:
                accepted.append(draws)

    if len(accepted) < PRE_FIT_MIN_SUCCESSFUL_CHAINS:
        raise RuntimeError(
            "Pre-fit did not produce enough successful chains: "
            f"{len(accepted)} < {PRE_FIT_MIN_SUCCESSFUL_CHAINS}"
        )
    return pd.concat(accepted, ignore_index=True)


def select_initial_values_from_prefit(
    pre_df: pd.DataFrame,
    fit_params: list[str],
) -> dict[str, float]:
    good_df = pre_df[
        pre_df["lp__"] > pre_df["lp__"].max() - PRE_FIT_LP_THRESHOLD_DROP
    ]
    if good_df.empty:
        raise RuntimeError("Pre-fit produced no draw within the configured lp__ threshold.")

    selected = good_df.sort_values("lp__", ascending=False).iloc[0]
    print(
        "Selected pre-fit initial point: "
        f"lp__ = {selected['lp__']:.2f}, "
        f"candidate_draws = {len(good_df)} / {len(pre_df)}"
    )
    return {param: float(selected[param]) for param in fit_params}


def prefixed_columns(
    df: pd.DataFrame,
    prefix: str,
    dimension: int | None = None,
) -> list[str]:
    if dimension is None and prefix in df.columns:
        return [prefix]
    if dimension is not None:
        columns = [f"{prefix}[{index}]" for index in range(1, dimension + 1)]
        missing = [column for column in columns if column not in df.columns]
        if missing:
            raise KeyError(f"Missing output columns for {prefix}: {', '.join(missing)}")
        return columns
    columns = [
        column
        for column in df.columns
        if column.startswith(prefix + "[")
    ]
    return sorted(columns, key=lambda name: int(name.rsplit("[", 1)[1][:-1]))


def selected_output_columns(df: pd.DataFrame, fit_params: list[str]) -> list[str]:
    columns = list(fit_params)
    for spec in diagnostic_specs(NP_MODEL).values():
        columns.extend(prefixed_columns(df, spec.prediction_name, spec.dimension))
        columns.extend(
            [
                f"T_obs_{spec.key}",
                f"T_rep_{spec.key}",
                f"log_lik_{spec.key}",
            ]
        )
    for prefix in ["rawBrRatio_func", "BrRatio_func", "BrRatio_data", "Vcb_Br"]:
        columns.extend(prefixed_columns(df, prefix))
    for scalar in ["Vcb_mean", "lp__", "chain__", "iter__", "draw__"]:
        if scalar in df.columns:
            columns.append(scalar)
    return list(dict.fromkeys(column for column in columns if column in df.columns))


def run_fit_case(
    model: CmdStanModel,
    case: FitCase,
    stan_data: dict[str, Any],
    fit_params: list[str],
    run_root: Path,
) -> FitResult:
    case_dir = run_root / case.name
    case_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== Fit case: {case.name} ===")
    print(f"Selected datasets: {list(case.selected_data)}")
    if case.alpha is not None:
        print(f"Covariance target: {case.target_dataset}, alpha = {case.alpha}")

    prefit_dir = case_dir / "cmdstan_prefit"
    pre_df = run_pre_fit(model, stan_data, output_dir=prefit_dir)
    good_init = select_initial_values_from_prefit(pre_df, fit_params)
    inits_set = [good_init] * int(FIT_SETTINGS["chains"])

    fit_dir = case_dir / "cmdstan_fit"
    fit_dir.mkdir(parents=True, exist_ok=True)
    fit = model.sample(
        data=stan_data,
        inits=inits_set,
        output_dir=fit_dir,
        sig_figs=STAN_SIG_FIGS,
        **FIT_SETTINGS,
    )
    draws = fit.draws_pd()

    print("Bad chain check:")
    good_chains = good_chain_ids(draws)
    if not good_chains:
        raise RuntimeError(f"No good chains survived for case {case.name}.")
    if len(good_chains) < 2:
        print("Warning: fewer than 2 good chains survived the lp__ check.")

    draws_good = draws[draws["chain__"].isin(good_chains)].copy()
    saved = draws_good[selected_output_columns(draws_good, fit_params)]
    output_file = case_dir / f"StanFit_{ANALYSIS_NAME}_{case.name}.csv"
    saved.to_csv(output_file, index=False)
    print(f"Fit case done: {output_file}")
    return FitResult(case, saved, stan_data, output_file)


def posterior_summary(values: np.ndarray, prefix: str = "") -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        f"{prefix}mean": float(np.mean(array)),
        f"{prefix}sd": float(np.std(array, ddof=1)),
        f"{prefix}median": float(np.median(array)),
        f"{prefix}q16": float(np.quantile(array, 0.16)),
        f"{prefix}q84": float(np.quantile(array, 0.84)),
        f"{prefix}q025": float(np.quantile(array, 0.025)),
        f"{prefix}q975": float(np.quantile(array, 0.975)),
    }


def log_mean_exp(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    maximum = float(np.max(values))
    return maximum + float(np.log(np.mean(np.exp(values - maximum))))


def prediction_matrix(draws: pd.DataFrame, spec: DiagnosticSpec) -> np.ndarray:
    columns = prefixed_columns(draws, spec.prediction_name, spec.dimension)
    return draws[columns].to_numpy(dtype=float)


def case_diagnostics(
    result: FitResult,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    group_rows: list[dict[str, Any]] = []
    bin_rows: list[dict[str, Any]] = []
    mode_rows: list[dict[str, Any]] = []
    specs = diagnostic_specs(NP_MODEL)
    active = set(result.case.selected_data)
    standardized_obs: list[np.ndarray] = []
    standardized_rep: list[np.ndarray] = []

    for key in DATASET_ORDER:
        spec = specs[key]
        observed = np.asarray(result.stan_data[spec.observed_name], dtype=float)
        covariance = np.asarray(result.stan_data[spec.covariance_name], dtype=float)
        prediction = prediction_matrix(result.draws, spec)
        residual = observed[None, :] - prediction

        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        if np.any(eigenvalues <= 0):
            raise ValueError(
                f"Non-positive covariance eigenvalue in {result.case.name}/{key}: "
                f"{eigenvalues.min()}"
            )
        ordinary_pull = residual / np.sqrt(np.diag(covariance))[None, :]
        eigen_pull = residual @ eigenvectors / np.sqrt(eigenvalues)[None, :]
        contribution = np.square(eigen_pull)
        t_from_modes = np.sum(contribution, axis=1)

        random_seed = int(FIT_SETTINGS.get("seed", 42)) + DATASET_ORDER.index(key)
        random_seed += sum(ord(char) for char in result.case.name)
        rng = np.random.default_rng(random_seed)
        standard_normal = rng.standard_normal((len(result.draws), spec.dimension))
        cholesky = np.linalg.cholesky(covariance)
        replicated_residual = standard_normal @ cholesky.T
        replicated_ordinary_pull = (
            replicated_residual / np.sqrt(np.diag(covariance))[None, :]
        )
        replicated_eigen_pull = (
            replicated_residual @ eigenvectors / np.sqrt(eigenvalues)[None, :]
        )

        t_obs = result.draws[f"T_obs_{key}"].to_numpy(dtype=float)
        t_rep = result.draws[f"T_rep_{key}"].to_numpy(dtype=float)
        max_t_error = float(np.max(np.abs(t_obs - t_from_modes)))
        scale = max(1.0, float(np.max(np.abs(t_obs))))
        if max_t_error > 1.0e-7 * scale:
            raise RuntimeError(
                f"Mahalanobis/eigenmode mismatch in {result.case.name}/{key}: "
                f"{max_t_error}"
            )

        log_lik = result.draws[f"log_lik_{key}"].to_numpy(dtype=float)
        standardized_obs.append(
            (t_obs - spec.dimension) / np.sqrt(2.0 * spec.dimension)
        )
        standardized_rep.append(
            (t_rep - spec.dimension) / np.sqrt(2.0 * spec.dimension)
        )
        group_row: dict[str, Any] = {
            "case": result.case.name,
            "kind": result.case.kind,
            "target_dataset": result.case.target_dataset,
            "alpha": result.case.alpha,
            "dataset": key,
            "in_likelihood": int(key in active),
            "dimension": spec.dimension,
            "pp_tail": float(np.mean(t_rep >= t_obs)),
            "pp_tail_maxbin": float(
                np.mean(
                    np.max(np.abs(replicated_ordinary_pull), axis=1)
                    >= np.max(np.abs(ordinary_pull), axis=1)
                )
            ),
            "pp_tail_maxeig": float(
                np.mean(
                    np.max(np.abs(replicated_eigen_pull), axis=1)
                    >= np.max(np.abs(eigen_pull), axis=1)
                )
            ),
            "lpd": log_mean_exp(log_lik),
            "condition_number": float(eigenvalues[-1] / eigenvalues[0]),
            "eigenvalue_min": float(eigenvalues[0]),
            "eigenvalue_max": float(eigenvalues[-1]),
            "T_identity_max_abs_error": max_t_error,
        }
        group_row.update(posterior_summary(t_obs, "T_obs_"))
        group_row.update(posterior_summary(t_rep, "T_rep_"))
        group_rows.append(group_row)

        for index in range(spec.dimension):
            row = {
                "case": result.case.name,
                "kind": result.case.kind,
                "target_dataset": result.case.target_dataset,
                "alpha": result.case.alpha,
                "dataset": key,
                "bin_index": index + 1,
            }
            row.update(posterior_summary(ordinary_pull[:, index], "pull_"))
            bin_rows.append(row)

        mean_total = float(np.mean(t_from_modes))
        for index in range(spec.dimension):
            loadings = eigenvectors[:, index]
            row = {
                "case": result.case.name,
                "kind": result.case.kind,
                "target_dataset": result.case.target_dataset,
                "alpha": result.case.alpha,
                "dataset": key,
                "mode_index": index + 1,
                "eigenvalue": float(eigenvalues[index]),
                "relative_eigenvalue": float(eigenvalues[index] / eigenvalues[-1]),
                "dominant_bin_index": int(np.argmax(np.abs(loadings))) + 1,
                "loadings": json.dumps(loadings.tolist()),
                "mean_chi_sq_fraction": (
                    float(np.mean(contribution[:, index]) / mean_total)
                    if mean_total > 0
                    else np.nan
                ),
            }
            row.update(posterior_summary(eigen_pull[:, index], "eigen_pull_"))
            row.update(posterior_summary(contribution[:, index], "contribution_"))
            mode_rows.append(row)

    max_obs = np.max(np.column_stack(standardized_obs), axis=1)
    max_rep = np.max(np.column_stack(standardized_rep), axis=1)
    global_row: dict[str, Any] = {
        "case": result.case.name,
        "kind": result.case.kind,
        "target_dataset": result.case.target_dataset,
        "alpha": result.case.alpha,
        "dataset": "GLOBAL_MAX",
        "in_likelihood": np.nan,
        "dimension": len(DATASET_ORDER),
        "pp_tail": float(np.mean(max_rep >= max_obs)),
        "pp_tail_maxbin": np.nan,
        "pp_tail_maxeig": np.nan,
        "lpd": np.nan,
        "condition_number": np.nan,
        "eigenvalue_min": np.nan,
        "eigenvalue_max": np.nan,
        "T_identity_max_abs_error": np.nan,
    }
    global_row.update(posterior_summary(max_obs, "T_obs_"))
    global_row.update(posterior_summary(max_rep, "T_rep_"))
    group_rows.append(global_row)

    return group_rows, bin_rows, mode_rows


def output_components(
    draws: pd.DataFrame,
    names: list[str],
) -> tuple[np.ndarray, list[str]]:
    arrays: list[np.ndarray] = []
    labels: list[str] = []
    for name in names:
        if name in draws.columns:
            arrays.append(draws[[name]].to_numpy(dtype=float))
            labels.append(name)
            continue
        columns = prefixed_columns(draws, name)
        if not columns:
            raise KeyError(f"No output columns found for {name}.")
        arrays.append(draws[columns].to_numpy(dtype=float))
        labels.extend(columns)
    if not arrays:
        return np.empty((len(draws), 0)), []
    return np.concatenate(arrays, axis=1), labels


def interval_width(values: np.ndarray) -> float:
    return float(np.quantile(values, 0.84) - np.quantile(values, 0.16))


def safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if abs(denominator) > 0 else np.nan


def lodo_influence_rows(
    all_result: FitResult,
    lodo_results: list[FitResult],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    all_values, labels = output_components(all_result.draws, PHYSICS_OUTPUTS)
    for result in lodo_results:
        lodo_values, lodo_labels = output_components(result.draws, PHYSICS_OUTPUTS)
        if labels != lodo_labels:
            raise RuntimeError("All-data and LODO physics-output columns differ.")
        for index, label in enumerate(labels):
            all_component = all_values[:, index]
            lodo_component = lodo_values[:, index]
            rows.append(
                {
                    "dataset": result.case.target_dataset,
                    "output": label,
                    "all_median": float(np.median(all_component)),
                    "lodo_median": float(np.median(lodo_component)),
                    "delta_median": float(
                        np.median(lodo_component) - np.median(all_component)
                    ),
                    "delta_in_all_sd": safe_ratio(
                        float(np.median(lodo_component) - np.median(all_component)),
                        float(np.std(all_component, ddof=1)),
                    ),
                    "all_width68": interval_width(all_component),
                    "lodo_width68": interval_width(lodo_component),
                    "width_ratio": safe_ratio(
                        interval_width(lodo_component),
                        interval_width(all_component),
                    ),
                }
            )
    return rows


def node_split_rows(
    target: str,
    only_result: FitResult,
    lodo_result: FitResult,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    only_values, labels = output_components(
        only_result.draws, NODE_SPLIT_SEPARATOR_PREFIXES
    )
    lodo_values, lodo_labels = output_components(
        lodo_result.draws, NODE_SPLIT_SEPARATOR_PREFIXES
    )
    if labels != lodo_labels:
        raise RuntimeError(f"Node-split separator columns differ for {target}.")

    sample_size = min(len(only_values), len(lodo_values))
    rng = np.random.default_rng(FIT_SETTINGS.get("seed", 42))
    only_indices = rng.choice(len(only_values), size=sample_size, replace=False)
    lodo_indices = rng.choice(len(lodo_values), size=sample_size, replace=False)
    delta_draws = only_values[only_indices] - lodo_values[lodo_indices]

    component_rows: list[dict[str, Any]] = []
    for index, label in enumerate(labels):
        delta = delta_draws[:, index]
        probability_positive = float(np.mean(delta >= 0))
        probability_negative = float(np.mean(delta <= 0))
        sd = float(np.std(delta, ddof=1))
        row = {
            "dataset": target,
            "separator": label,
            "only_median": float(np.median(only_values[:, index])),
            "lodo_median": float(np.median(lodo_values[:, index])),
            "Z_split": float(np.mean(delta) / sd) if sd > 0 else np.nan,
            "p_split_two": min(
                1.0, 2.0 * min(probability_positive, probability_negative)
            ),
        }
        row.update(posterior_summary(delta, "delta_"))
        component_rows.append(row)

    delta_mean = np.mean(only_values, axis=0) - np.mean(lodo_values, axis=0)
    only_cov = np.atleast_2d(np.cov(only_values, rowvar=False, ddof=1))
    lodo_cov = np.atleast_2d(np.cov(lodo_values, rowvar=False, ddof=1))
    delta_cov = only_cov + lodo_cov
    eigenvalues = np.linalg.eigvalsh(delta_cov)
    threshold = NODE_SPLIT_RANK_TOL * max(float(eigenvalues[-1]), 1.0)
    rank = int(np.count_nonzero(eigenvalues > threshold))
    inverse = np.linalg.pinv(delta_cov, rcond=NODE_SPLIT_RANK_TOL)
    q_split = float(delta_mean @ inverse @ delta_mean)
    try:
        from scipy.stats import chi2

        p_multi = float(chi2.sf(q_split, rank)) if rank > 0 else np.nan
    except ImportError:
        p_multi = np.nan

    summary = {
        "dataset": target,
        "separator_dimension": len(labels),
        "retained_rank": rank,
        "Q_split": q_split,
        "p_split_multi_gaussian": p_multi,
    }
    return component_rows, summary


def covariance_scan_rows(
    all_result: FitResult,
    scan_results: list[FitResult],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    all_values, labels = output_components(all_result.draws, PHYSICS_OUTPUTS)
    for result in scan_results:
        scan_values, scan_labels = output_components(result.draws, PHYSICS_OUTPUTS)
        if labels != scan_labels:
            raise RuntimeError("All-data and covariance-scan output columns differ.")
        for index, label in enumerate(labels):
            nominal = all_values[:, index]
            scan = scan_values[:, index]
            rows.append(
                {
                    "dataset": result.case.target_dataset,
                    "alpha": result.case.alpha,
                    "output": label,
                    "nominal_median": float(np.median(nominal)),
                    "scan_median": float(np.median(scan)),
                    "delta_median": float(np.median(scan) - np.median(nominal)),
                    "delta_in_nominal_sd": safe_ratio(
                        float(np.median(scan) - np.median(nominal)),
                        float(np.std(nominal, ddof=1)),
                    ),
                    "nominal_width68": interval_width(nominal),
                    "scan_width68": interval_width(scan),
                    "width_ratio": safe_ratio(
                        interval_width(scan),
                        interval_width(nominal),
                    ),
                }
            )
    return rows


def covariance_scan_lodo_rows(
    target: str,
    lodo_result: FitResult,
    base_data: dict[str, Any],
) -> list[dict[str, Any]]:
    spec = diagnostic_specs(NP_MODEL)[target]
    observed = np.asarray(base_data[spec.observed_name], dtype=float)
    nominal = np.asarray(base_data[spec.covariance_name], dtype=float)
    prediction = prediction_matrix(lodo_result.draws, spec)
    residual = observed[None, :] - prediction
    rows: list[dict[str, Any]] = []

    for alpha in COVARIANCE_ALPHAS:
        covariance = shrink_covariance(nominal, float(alpha))
        cholesky = np.linalg.cholesky(covariance)
        whitened = np.linalg.solve(cholesky, residual.T).T
        t_obs = np.sum(np.square(whitened), axis=1)

        random_seed = int(FIT_SETTINGS.get("seed", 42))
        random_seed += DATASET_ORDER.index(target) + int(round(1000 * alpha))
        rng = np.random.default_rng(random_seed)
        t_rep = rng.chisquare(spec.dimension, size=len(t_obs))

        sign, log_determinant = np.linalg.slogdet(covariance)
        if sign <= 0:
            raise ValueError(f"Non-positive covariance determinant for {target}.")
        log_lik = -0.5 * (
            spec.dimension * np.log(2.0 * np.pi) + log_determinant + t_obs
        )
        eigenvalues = np.linalg.eigvalsh(covariance)
        row: dict[str, Any] = {
            "dataset": target,
            "alpha": float(alpha),
            "p_LODO": float(np.mean(t_rep >= t_obs)),
            "lpd": log_mean_exp(log_lik),
            "condition_number": float(eigenvalues[-1] / eigenvalues[0]),
            "eigenvalue_min": float(eigenvalues[0]),
            "eigenvalue_max": float(eigenvalues[-1]),
        }
        row.update(posterior_summary(t_obs, "T_obs_"))
        row.update(posterior_summary(t_rep, "T_rep_"))
        rows.append(row)
    return rows


def write_metadata(
    path: Path,
    generated: dict[str, Any],
    cases: list[FitCase],
) -> None:
    metadata = {
        "analysis_name": ANALYSIS_NAME,
        "ff_model": FF_MODEL,
        "ff_nf": FF_NF,
        "np_model": NP_MODEL,
        "hqet_model": HQET_MODEL,
        "ub_model": effective_ub_model(),
        "dataset_order": DATASET_ORDER,
        "cases": [
            {
                "name": case.name,
                "kind": case.kind,
                "selected_data": list(case.selected_data),
                "target_dataset": case.target_dataset,
                "alpha": case.alpha,
            }
            for case in cases
        ],
        "physics_outputs": PHYSICS_OUTPUTS,
        "node_split_separator_prefixes": NODE_SPLIT_SEPARATOR_PREFIXES,
        "pre_fit_settings": PRE_FIT_SETTINGS,
        "fit_settings": FIT_SETTINGS,
        "generated": generated,
    }
    path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def run_extended_analysis() -> dict[str, Path]:
    generated = run_python_generator()
    fit_dir = Path(generated["FitDirectory"])
    stan_file = Path(generated["StanFile"])
    param_file = Path(generated["ParamFile"])
    data_file = Path(generated["DataFile"])

    base_data = load_stan_data(data_file)
    fit_params = json.loads(param_file.read_text(encoding="utf-8"))
    cases = build_fit_cases()
    if not cases:
        raise RuntimeError("No fit cases are enabled.")

    # This is the only CmdStanModel construction.  Every case below reuses
    # the same compiled executable and changes only standata.
    print("Loading/compiling the extended master Stan model...")
    model = CmdStanModel(stan_file=stan_file)

    timestamp = datetime.today().strftime("%m%d_%H%M%S")
    run_root = fit_dir / f"extended_runs_{timestamp}"
    run_root.mkdir(parents=True, exist_ok=True)
    write_metadata(run_root / "run_metadata.json", generated, cases)

    results: list[FitResult] = []
    group_rows: list[dict[str, Any]] = []
    bin_rows: list[dict[str, Any]] = []
    mode_rows: list[dict[str, Any]] = []

    for case in cases:
        case_data = prepare_case_data(base_data, case)
        result = run_fit_case(model, case, case_data, fit_params, run_root)
        results.append(result)
        groups, bins, modes = case_diagnostics(result)
        group_rows.extend(groups)
        bin_rows.extend(bins)
        mode_rows.extend(modes)

    result_by_name = {result.case.name: result for result in results}
    output_files: dict[str, Path] = {}

    dataset_file = run_root / "dataset_diagnostics.csv"
    pd.DataFrame(group_rows).to_csv(dataset_file, index=False)
    output_files["dataset_diagnostics"] = dataset_file

    bin_file = run_root / "bin_diagnostics.csv"
    pd.DataFrame(bin_rows).to_csv(bin_file, index=False)
    output_files["bin_diagnostics"] = bin_file

    mode_file = run_root / "eigenmode_diagnostics.csv"
    pd.DataFrame(mode_rows).to_csv(mode_file, index=False)
    output_files["eigenmode_diagnostics"] = mode_file

    all_result = result_by_name.get("AllData")
    lodo_results = [result for result in results if result.case.kind == "lodo"]
    if all_result is not None and lodo_results:
        influence_file = run_root / "lodo_influence.csv"
        pd.DataFrame(lodo_influence_rows(all_result, lodo_results)).to_csv(
            influence_file, index=False
        )
        output_files["lodo_influence"] = influence_file

    node_component_rows: list[dict[str, Any]] = []
    node_summary_rows: list[dict[str, Any]] = []
    for target in canonical_dataset_list(
        NODE_SPLIT_DATASETS if RUN_NODE_SPLIT else [], "node-split"
    ):
        components, summary = node_split_rows(
            target,
            result_by_name[f"Only_{target}"],
            result_by_name[f"LODO_{target}"],
        )
        node_component_rows.extend(components)
        node_summary_rows.append(summary)
    if node_component_rows:
        node_component_file = run_root / "node_split_components.csv"
        pd.DataFrame(node_component_rows).to_csv(node_component_file, index=False)
        output_files["node_split_components"] = node_component_file
        node_summary_file = run_root / "node_split_summary.csv"
        pd.DataFrame(node_summary_rows).to_csv(node_summary_file, index=False)
        output_files["node_split_summary"] = node_summary_file

    scan_results = [
        result for result in results if result.case.kind == "covariance_scan"
    ]
    if all_result is not None and scan_results:
        scan_file = run_root / "covariance_scan_outputs.csv"
        pd.DataFrame(covariance_scan_rows(all_result, scan_results)).to_csv(
            scan_file, index=False
        )
        output_files["covariance_scan_outputs"] = scan_file

        scan_lodo_rows: list[dict[str, Any]] = []
        for target in canonical_dataset_list(
            COVARIANCE_SCAN_DATASETS, "covariance-scan"
        ):
            scan_lodo_rows.extend(
                covariance_scan_lodo_rows(
                    target,
                    result_by_name[f"LODO_{target}"],
                    base_data,
                )
            )
        scan_lodo_file = run_root / "covariance_scan_lodo.csv"
        pd.DataFrame(scan_lodo_rows).to_csv(scan_lodo_file, index=False)
        output_files["covariance_scan_lodo"] = scan_lodo_file

    output_files["metadata"] = run_root / "run_metadata.json"
    print(f"\nExtended analysis done: {run_root}")
    return output_files


def run_fit() -> Path:
    """Compatibility entry point returning the all-data draw CSV."""
    global RUN_ALL_DATA, RUN_LODO, RUN_NODE_SPLIT, RUN_COVARIANCE_SCAN
    RUN_ALL_DATA = True
    RUN_LODO = False
    RUN_NODE_SPLIT = False
    RUN_COVARIANCE_SCAN = False
    outputs = run_extended_analysis()
    metadata_path = outputs["metadata"]
    run_root = metadata_path.parent
    candidates = sorted((run_root / "AllData").glob("StanFit_*.csv"))
    if not candidates:
        raise RuntimeError("All-data fit output was not created.")
    return candidates[-1]


if __name__ == "__main__":
    run_extended_analysis()
