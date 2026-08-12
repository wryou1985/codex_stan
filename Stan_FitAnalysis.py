from __future__ import annotations

import json
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from cmdstanpy import CmdStanModel

from Stan_CodeGenerator import generate_from_config


PROJECT_DIR = Path(__file__).resolve().parent

# Main analysis controls.  CutControl is intentionally configured here.
ANALYSIS_NAME = "CasePDG-noUB_SM_BGLnf1"
FF_MODEL = "BGL"
FF_NF = 1
NP_MODEL = "SM"
HQET_MODEL = "3/2/1"
UB_MODEL = "no-UB"  # "no-UB", "soft-UB", or "hard-UB"
CUT_CONTROL = 1.0e-4
BELLE_FIT_OPTION = "sc"  # "sc", "re_small", "re_large", "de", "pdg", or "rms"

SELECTED_DATA = [
    "MILC15",
    "MILC21",
    "JLQCD23",
    "HPQCD23",
    "LCSR18",
    "LCSR23",
    "BelleFITD",
    "BelleFITDstW",
    "BelleFITDstCosL",
    "BelleFITDstCosV",
    "BelleFITDstChi",
]

INPUT_ROOT = PROJECT_DIR / "Stan_inputs"
OUTPUT_ROOT = PROJECT_DIR / "Stan_outputs"
FIT_FUNCTION_CACHE_ROOT = INPUT_ROOT / "StanFunctionCache"
OUTPUT_MODE = "full"

UB_MAX = np.ones(4)
UB_SIGMA = np.full(4, 0.05)

PRE_FIT_SETTINGS = {
    "chains": 10,
    "iter_warmup": 500,
    "iter_sampling": 200,
    "adapt_delta": 0.95,
    "max_treedepth": 12,
}
PRE_FIT_CHAIN_TIMEOUT = 600.0  # seconds per chain; set to None to disable.
PRE_FIT_PARALLEL_CHAINS = None  # None uses min(chains, os.cpu_count()).
PRE_FIT_MIN_SUCCESSFUL_CHAINS = 1
PRE_FIT_LP_THRESHOLD_DROP = 200.0

FIT_SETTINGS = {
    "chains": 3,
    "iter_warmup": 3000,
    "iter_sampling": 6000,
    "adapt_delta": 0.99,
    "max_treedepth": 15,
    "seed": 42,
}


def effective_ub_model() -> str:
    return UB_MODEL if FF_MODEL == "HQET" else "no-UB"


def build_generator_config() -> dict:
    return {
        "WorkingDirectory": str(PROJECT_DIR),
        "AnalysisName": ANALYSIS_NAME,
        "FFModel": FF_MODEL,
        "FFnf": FF_NF,
        "NPModel": NP_MODEL,
        "HQETModel": HQET_MODEL,
        "UBModel": effective_ub_model(),
        "CutControl": CUT_CONTROL,
        "SelectedData": SELECTED_DATA,
        "GeneratedOutputs": "Automatic",
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


def run_python_generator() -> dict:
    config = build_generator_config()
    config_path = INPUT_ROOT / f"config_{ANALYSIS_NAME}.json"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
    return generate_from_config(config)


def load_stan_data(data_file: Path) -> dict:
    with data_file.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)

    stan_data = {key: np.asarray(value) for key, value in raw.items()}
    ub_model = effective_ub_model()
    if ub_model == "soft-UB":
        stan_data.update({"UBmax": UB_MAX, "UB_sigma": UB_SIGMA})
    elif ub_model == "hard-UB":
        stan_data.update({"UBmax": UB_MAX})
    return stan_data


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


def selected_observable_columns(df: pd.DataFrame) -> list[str]:
    obs_option = {
        "full": [
            "MILC15func",
            "MILC21func",
            "HPQCD23func",
            "JLQCD23func",
            "LCSR18func",
            "LCSR23func",
            "rawBrRatio_func",
            "rawBelleD_func",
            "rawBelleDst_w_func",
            "rawBelleDst_cosL_func",
            "rawBelleDst_cosV_func",
            "rawBelleDst_chi_func",
            "BrRatio_func",
            "BelleD_func",
            "BelleDst_w_func",
            "BelleDst_cosL_func",
            "BelleDst_cosV_func",
            "BelleDst_chi_func",
            "RD_func",
            "chi_sq",
            "Vcb_Br",
        ],
        "light": ["MILC15func", "MILC21func", "HPQCD23func", "JLQCD23func", "LCSR18func", "LCSR23func"],
    }
    wanted = obs_option.get(OUTPUT_MODE, [])
    return [col for col in df.columns if any(col == obs or col.startswith(obs + "[") for obs in wanted)]


def run_pre_fit(model: CmdStanModel, stan_data: dict, output_dir: Path) -> pd.DataFrame:
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
                **settings,
            )
        except Exception as exc:
            print(f"Pre-fit chain {chain_id} discarded: {type(exc).__name__}: {exc}")
            return None

        draws = fit.draws_pd()
        lp_max = draws["lp__"].max() if "lp__" in draws else float("nan")
        print(f"Pre-fit chain {chain_id} accepted: sample_length = {len(draws)}, lp_max = {lp_max:.2f}")
        return draws

    accepted: list[pd.DataFrame] = []
    with ThreadPoolExecutor(max_workers=parallel_chains) as executor:
        futures = {executor.submit(run_chain, chain_id): chain_id for chain_id in range(1, chains + 1)}
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


def select_initial_values_from_prefit(pre_df: pd.DataFrame, fit_params: list[str]) -> dict:
    good_df = pre_df[pre_df["lp__"] > pre_df["lp__"].max() - PRE_FIT_LP_THRESHOLD_DROP]
    if good_df.empty:
        raise RuntimeError("Pre-fit produced no draw within the configured lp__ threshold.")

    selected = good_df.sort_values("lp__", ascending=False).iloc[0]
    print(
        "Selected pre-fit initial point: "
        f"lp__ = {selected['lp__']:.2f}, "
        f"candidate_draws = {len(good_df)} / {len(pre_df)}"
    )
    return selected[fit_params].to_dict()


def run_fit() -> Path:
    import arviz as az

    generated = run_python_generator()
    fit_dir = Path(generated["FitDirectory"])
    stan_file = Path(generated["StanFile"])
    param_file = Path(generated["ParamFile"])
    data_file = Path(generated["DataFile"])

    stan_data = load_stan_data(data_file)
    fit_params = json.loads(param_file.read_text(encoding="utf-8"))

    print("Compiling Stan model...")
    model = CmdStanModel(stan_file=stan_file)

    print("Pre-fit run...")
    pre_fit_output_dir = fit_dir / f"cmdstan_prefit_{datetime.today().strftime('%m%d_%H%M')}"
    pre_df = run_pre_fit(model, stan_data, output_dir=pre_fit_output_dir)
    good_init = select_initial_values_from_prefit(pre_df, fit_params)
    inits_set = [good_init] * FIT_SETTINGS["chains"]

    print("Fit analysis run...")
    fit = model.sample(data=stan_data, inits=inits_set, **FIT_SETTINGS)
    df = fit.draws_pd()

    print("Bad chain check:")
    good_chains = good_chain_ids(df)
    if len(good_chains) < 2:
        print("Warning: fewer than 2 good chains survived the lp__ check.")

    df_good = df[df["chain__"].isin(good_chains)]
    fit_all = fit_params + selected_observable_columns(df) + ["lp__", "chi_sq_total", "Vcb_mean"]
    fit_all = [col for col in fit_all if col in df_good.columns]
    params_df = df_good[fit_all]

    output_dir = OUTPUT_ROOT
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.today().strftime("%m%d_%H%M")
    output_file = output_dir / f"StanFit_{ANALYSIS_NAME}_{timestamp}.csv"
    params_df.to_csv(output_file, index=False)

    print("Display Fit Result:")
    print(params_df.head())

    idata = az.from_cmdstanpy(fit)
    good_chains_idata = [chain - 1 for chain in good_chains]
    good_idata = idata.sel(chain=good_chains_idata)
    print("R-hat values (FF parameters):")
    print(az.summary(good_idata, var_names=fit_params, round_to=4))
    if OUTPUT_MODE == "full":
        print("R-hat values (RD, Br, Vcb, chi_sq):")
        print(az.summary(good_idata, var_names=["RD_func", "BrRatio_func", "Vcb_Br", "Vcb_mean", "chi_sq_total"], round_to=4))

    print(f"Fit analysis done: {output_file}")
    print(f"Generated Stan files: {fit_dir}")
    return output_file


if __name__ == "__main__":
    run_fit()
