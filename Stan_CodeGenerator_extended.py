from __future__ import annotations

import json
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Stan_DataPoint import prepare_data_point_inputs
from Stan_FitFunction import FitFunctionSet, load_fit_functions
from Stan_IntFunctions import canonical_np_model


DEFAULT_SELECTED_DATA = [
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

DATASET_ORDER = list(DEFAULT_SELECTED_DATA)
DATASET_INDEX = {key: index for index, key in enumerate(DATASET_ORDER, start=1)}

FF_PARAMETER_DECLS = {
    ("BSZ", 2): ["afp0", "afp1", "afp2", "af01", "af02", "aV0", "aV1", "aV2", "aA10", "aA11", "aA12", "aA121", "aA122", "aA01", "aA02"],
    ("BSZ", 1): ["afp0", "afp1", "af01", "aV0", "aV1", "aA10", "aA11", "aA121", "aA01"],
    ("BGL", 2): ["bfp0", "bfp1", "bfp2", "bf01", "bf02", "bg0", "bg1", "bg2", "bf0", "bf1", "bf2", "bF11", "bF12", "bF21", "bF22"],
    ("BGL", 1): ["bfp0", "bfp1", "bf01", "bg0", "bg1", "bf0", "bf1", "bF11", "bF21"],
    ("HQET", 2): ["cxi1", "cxi2", "cxi3", "cet0", "cet1", "cet2", "cci20", "cci21", "cci22", "cci31", "cci32", "cel10", "cel11", "cel20", "cel21", "cel30", "cel31", "cel40", "cel41", "cel50", "cel51", "cel60", "cel61"],
    ("HQET", 1): ["cxi1", "cxi2", "cet0", "cet1", "cci20", "cci21", "cci31", "cel10", "cel20", "cel30", "cel40", "cel50", "cel60"],
}

FF_TENSOR_PARAMETER_DECLS = {
    ("BSZ", 2): ["afT0", "afT1", "afT2", "aT10", "aT11", "aT12", "aT21", "aT22", "aT231", "aT232"],
    ("BSZ", 1): ["afT0", "afT1", "aT10", "aT11", "aT21", "aT231"],
    ("BGL", 2): ["bfT0", "bfT1", "bfT2", "bT10", "bT11", "bT12", "bT21", "bT22", "bT231", "bT232"],
    ("BGL", 1): ["bfT0", "bfT1", "bT10", "bT11", "bT21", "bT231"],
}

HQET_QCDSR_PRIORS = {
    1: [
        "cxi1 ~ normal(0, 1);",
        "cxi2 ~ normal(0, 1);",
        "cet0 ~ normal(0.62, 0.12);",
        "cet1 ~ normal(0.04, 0.03);",
        "cci20 ~ normal(-0.06, 0.02);",
        "cci21 ~ normal(0.00, 0.02);",
        "cci31 ~ normal(0.04, 0.03);",
        "cel10 ~ normal(0, 1);",
        "cel20 ~ normal(0, 1);",
        "cel30 ~ normal(0, 1);",
        "cel40 ~ normal(0, 1);",
        "cel50 ~ normal(0, 1);",
        "cel60 ~ normal(0, 1);",
    ],
    2: [
        "cxi1 ~ normal(0, 1);",
        "cxi2 ~ normal(0, 1);",
        "cxi3 ~ normal(0, 1);",
        "cet0 ~ normal(0.62, 0.12);",
        "cet1 ~ normal(0.04, 0.03);",
        "cet2 ~ normal(0.05, 0.07);",
        "cci20 ~ normal(-0.06, 0.02);",
        "cci21 ~ normal(0.00, 0.02);",
        "cci22 ~ normal(-0.01, 0.02);",
        "cci31 ~ normal(0.04, 0.03);",
        "cci32 ~ normal(-0.03, 0.05);",
        "cel10 ~ normal(0, 1);",
        "cel11 ~ normal(0, 1);",
        "cel20 ~ normal(0, 1);",
        "cel21 ~ normal(0, 1);",
        "cel30 ~ normal(0, 1);",
        "cel31 ~ normal(0, 1);",
        "cel40 ~ normal(0, 1);",
        "cel41 ~ normal(0, 1);",
        "cel50 ~ normal(0, 1);",
        "cel51 ~ normal(0, 1);",
        "cel60 ~ normal(0, 1);",
        "cel61 ~ normal(0, 1);",
    ],
}


@dataclass(frozen=True)
class DataModule:
    key: str
    data_decl: str
    function_decl: str
    function_body_key: str | None
    model_term: str
    chi_sq_expr: str | None
    raw_output_key: str | None = None


@dataclass(frozen=True)
class DiagnosticSpec:
    key: str
    dimension: int
    observed_name: str
    covariance_name: str
    prediction_name: str


def unique(items: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for item in items:
        if item not in seen:
            seen.add(item)
            result.append(item)
    return result


def join_lines(parts: Iterable[str | list[str]]) -> str:
    lines: list[str] = []
    for part in parts:
        if isinstance(part, list):
            lines.extend(line for line in part if line)
        elif part:
            lines.append(part)
    return "\n".join(lines)


def effective_ub_model(ff_model: str, ub_model: str) -> str:
    return ub_model if ff_model == "HQET" else "no-UB"


def required_data_keys(selected_data: list[str]) -> list[str]:
    # The extended model always declares every dataset.  SelectedData changes
    # only the run-time use_dataset mask, so one compiled executable can be
    # reused for all all-data, LODO, and dataset-only fits.
    del selected_data
    return [*DATASET_ORDER, "BrRatio2025"]


def data_modules(np_model: str) -> dict[str, DataModule]:
    tensor = np_model == "T"
    hpqcd_decl = (
        "vector[35] bglHPQCD23cent;\ncov_matrix[35] bglHPQCD23cov;"
        if tensor
        else "vector[20] bglHPQCD23centNoTensor;\ncov_matrix[20] bglHPQCD23covNoTensor;"
    )
    hpqcd_model = (
        "bglHPQCD23cent ~ multi_normal(HPQCD23func, bglHPQCD23cov);"
        if tensor
        else "bglHPQCD23centNoTensor ~ multi_normal(HPQCD23func, bglHPQCD23covNoTensor);"
    )
    hpqcd_chi = (
        "dot_product(bglHPQCD23cent - HPQCD23func, inverse(bglHPQCD23cov) * (bglHPQCD23cent - HPQCD23func))"
        if tensor
        else "dot_product(bglHPQCD23centNoTensor - HPQCD23func, inverse(bglHPQCD23covNoTensor) * (bglHPQCD23centNoTensor - HPQCD23func))"
    )
    lcsr18_decl = (
        "vector[37] bszLCSR18cent;\ncov_matrix[37] bszLCSR18cov;"
        if tensor
        else "vector[22] bszLCSR18centNoTensor;\ncov_matrix[22] bszLCSR18covNoTensor;"
    )
    lcsr18_model = (
        "bszLCSR18cent ~ multi_normal(LCSR18func, bszLCSR18cov);"
        if tensor
        else "bszLCSR18centNoTensor ~ multi_normal(LCSR18func, bszLCSR18covNoTensor);"
    )
    lcsr18_chi = (
        "dot_product(bszLCSR18cent - LCSR18func, inverse(bszLCSR18cov) * (bszLCSR18cent - LCSR18func))"
        if tensor
        else "dot_product(bszLCSR18centNoTensor - LCSR18func, inverse(bszLCSR18covNoTensor) * (bszLCSR18centNoTensor - LCSR18func))"
    )

    return {
        "BrRatio2025": DataModule("BrRatio2025", "vector[2] BrRatio2025cent;\ncov_matrix[2] BrRatio2025cov;", "", None, "", None),
        "MILC15": DataModule("MILC15", "vector[6] MILC15cent;\ncov_matrix[6] MILC15cov;", "vector[6] MILC15func;", "MILC15", "MILC15cent ~ multi_normal(MILC15func, MILC15cov);", "dot_product(MILC15cent - MILC15func, inverse(MILC15cov) * (MILC15cent - MILC15func))"),
        "MILC21": DataModule("MILC21", "vector[12] bglMILC21cent;\ncov_matrix[12] bglMILC21cov;", "vector[12] MILC21func;", "MILC21", "bglMILC21cent ~ multi_normal(MILC21func, bglMILC21cov);", "dot_product(bglMILC21cent - MILC21func, inverse(bglMILC21cov) * (bglMILC21cent - MILC21func))"),
        "JLQCD23": DataModule("JLQCD23", "vector[12] bglJLQCD23cent;\ncov_matrix[12] bglJLQCD23cov;", "vector[12] JLQCD23func;", "JLQCD23", "bglJLQCD23cent ~ multi_normal(JLQCD23func, bglJLQCD23cov);", "dot_product(bglJLQCD23cent - JLQCD23func, inverse(bglJLQCD23cov) * (bglJLQCD23cent - JLQCD23func))"),
        "HPQCD23": DataModule("HPQCD23", hpqcd_decl, f"vector[{35 if tensor else 20}] HPQCD23func;", "HPQCD23", hpqcd_model, hpqcd_chi),
        "LCSR18": DataModule("LCSR18", lcsr18_decl, f"vector[{37 if tensor else 22}] LCSR18func;", "LCSR18", lcsr18_model, lcsr18_chi),
        "LCSR23": DataModule("LCSR23", "vector[34] bszLCSR23cent;\ncov_matrix[34] bszLCSR23cov;", "vector[34] LCSR23func;", "LCSR23", "bszLCSR23cent ~ multi_normal(LCSR23func, bszLCSR23cov);", "dot_product(bszLCSR23cent - LCSR23func, inverse(bszLCSR23cov) * (bszLCSR23cent - LCSR23func))"),
        "BelleFITD": DataModule("BelleFITD", "vector[10] BelleFITcent_D;\ncov_matrix[10] BelleFITcov_D;", "vector[10] BelleD_func;", "BelleFITD", "BelleFITcent_D ~ multi_normal(BelleD_func, BelleFITcov_D);", "dot_product(BelleFITcent_D - BelleD_func, inverse(BelleFITcov_D) * (BelleFITcent_D - BelleD_func))", "RawBelleD"),
        "BelleFITDstW": DataModule("BelleFITDstW", "vector[10] BelleFITcent_w;\ncov_matrix[10] BelleFITcov_w;", "vector[10] BelleDst_w_func;", "BelleFITDstW", "BelleFITcent_w ~ multi_normal(BelleDst_w_func, BelleFITcov_w);", "dot_product(BelleFITcent_w - BelleDst_w_func, inverse(BelleFITcov_w) * (BelleFITcent_w - BelleDst_w_func))", "RawBelleDstW"),
        "BelleFITDstCosL": DataModule("BelleFITDstCosL", "vector[10] BelleFITcent_cosL;\ncov_matrix[10] BelleFITcov_cosL;", "vector[10] BelleDst_cosL_func;", "BelleFITDstCosL", "BelleFITcent_cosL ~ multi_normal(BelleDst_cosL_func, BelleFITcov_cosL);", "dot_product(BelleFITcent_cosL - BelleDst_cosL_func, inverse(BelleFITcov_cosL) * (BelleFITcent_cosL - BelleDst_cosL_func))", "RawBelleDstCosL"),
        "BelleFITDstCosV": DataModule("BelleFITDstCosV", "vector[10] BelleFITcent_cosV;\ncov_matrix[10] BelleFITcov_cosV;", "vector[10] BelleDst_cosV_func;", "BelleFITDstCosV", "BelleFITcent_cosV ~ multi_normal(BelleDst_cosV_func, BelleFITcov_cosV);", "dot_product(BelleFITcent_cosV - BelleDst_cosV_func, inverse(BelleFITcov_cosV) * (BelleFITcent_cosV - BelleDst_cosV_func))", "RawBelleDstCosV"),
        "BelleFITDstChi": DataModule("BelleFITDstChi", "vector[10] BelleFITcent_chi;\ncov_matrix[10] BelleFITcov_chi;", "vector[10] BelleDst_chi_func;", "BelleFITDstChi", "BelleFITcent_chi ~ multi_normal(BelleDst_chi_func, BelleFITcov_chi);", "dot_product(BelleFITcent_chi - BelleDst_chi_func, inverse(BelleFITcov_chi) * (BelleFITcent_chi - BelleDst_chi_func))", "RawBelleDstChi"),
    }


def diagnostic_specs(np_model: str) -> dict[str, DiagnosticSpec]:
    tensor = np_model == "T"
    return {
        "MILC15": DiagnosticSpec("MILC15", 6, "MILC15cent", "MILC15cov", "MILC15func"),
        "MILC21": DiagnosticSpec("MILC21", 12, "bglMILC21cent", "bglMILC21cov", "MILC21func"),
        "JLQCD23": DiagnosticSpec("JLQCD23", 12, "bglJLQCD23cent", "bglJLQCD23cov", "JLQCD23func"),
        "HPQCD23": DiagnosticSpec(
            "HPQCD23",
            35 if tensor else 20,
            "bglHPQCD23cent" if tensor else "bglHPQCD23centNoTensor",
            "bglHPQCD23cov" if tensor else "bglHPQCD23covNoTensor",
            "HPQCD23func",
        ),
        "LCSR18": DiagnosticSpec(
            "LCSR18",
            37 if tensor else 22,
            "bszLCSR18cent" if tensor else "bszLCSR18centNoTensor",
            "bszLCSR18cov" if tensor else "bszLCSR18covNoTensor",
            "LCSR18func",
        ),
        "LCSR23": DiagnosticSpec("LCSR23", 34, "bszLCSR23cent", "bszLCSR23cov", "LCSR23func"),
        "BelleFITD": DiagnosticSpec("BelleFITD", 10, "BelleFITcent_D", "BelleFITcov_D", "BelleD_func"),
        "BelleFITDstW": DiagnosticSpec("BelleFITDstW", 10, "BelleFITcent_w", "BelleFITcov_w", "BelleDst_w_func"),
        "BelleFITDstCosL": DiagnosticSpec(
            "BelleFITDstCosL", 10, "BelleFITcent_cosL", "BelleFITcov_cosL", "BelleDst_cosL_func"
        ),
        "BelleFITDstCosV": DiagnosticSpec(
            "BelleFITDstCosV", 10, "BelleFITcent_cosV", "BelleFITcov_cosV", "BelleDst_cosV_func"
        ),
        "BelleFITDstChi": DiagnosticSpec(
            "BelleFITDstChi", 10, "BelleFITcent_chi", "BelleFITcov_chi", "BelleDst_chi_func"
        ),
    }


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)
        handle.write("\n")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.read_text(encoding="utf-8") == text:
        return
    path.write_text(text, encoding="utf-8")


def slice_matrix(matrix: list[list[float]], start: int, stop: int) -> list[list[float]]:
    return [row[start:stop] for row in matrix[start:stop]]


def build_standata(selected_data: list[str], input_directory: Path, np_model: str) -> dict[str, Any]:
    unknown = sorted(set(selected_data) - set(DATASET_ORDER))
    if unknown:
        raise ValueError(f"Unknown SelectedData key(s): {', '.join(unknown)}")
    all_data = required_data_keys(selected_data)
    result: dict[str, Any] = {}

    source = {
        "MILC15": load_json(input_directory / "standata_MILC15data.json"),
        "MILC21": load_json(input_directory / "standata_bglMILC21data.json"),
        "JLQCD23": load_json(input_directory / "standata_bglJLQCD23data.json"),
        "HPQCD23": load_json(input_directory / "standata_bglHPQCD23data.json"),
        "LCSR18": load_json(input_directory / "standata_bszLCSR18data.json"),
        "LCSR23": load_json(input_directory / "standata_bszLCSR23data.json"),
        "BrRatio2025": load_json(input_directory / "standata_BrRatio2025data.json"),
    }
    belle = load_json(input_directory / "standata_BellFITdata.json")

    for key in all_data:
        if key in {"MILC15", "MILC21", "JLQCD23", "HPQCD23", "LCSR18", "LCSR23", "BrRatio2025"}:
            if key == "HPQCD23" and np_model != "T":
                result["bglHPQCD23centNoTensor"] = source[key]["bglHPQCD23centNoTensor"]
                result["bglHPQCD23covNoTensor"] = source[key]["bglHPQCD23covNoTensor"]
            elif key == "HPQCD23":
                result["bglHPQCD23cent"] = source[key]["bglHPQCD23cent"]
                result["bglHPQCD23cov"] = source[key]["bglHPQCD23cov"]
            elif key == "LCSR18" and np_model != "T":
                result["bszLCSR18centNoTensor"] = source[key]["bszLCSR18centNoTensor"]
                result["bszLCSR18covNoTensor"] = source[key]["bszLCSR18covNoTensor"]
            elif key == "LCSR18":
                result["bszLCSR18cent"] = source[key]["bszLCSR18cent"]
                result["bszLCSR18cov"] = source[key]["bszLCSR18cov"]
            else:
                result.update(source[key])
        elif key == "BelleFITD":
            result["BelleFITcent_D"] = belle["BelleFITcent_D"]
            result["BelleFITcov_D"] = belle["BelleFITcov_D"]
        elif key == "BelleFITDstW":
            result["BelleFITcent_w"] = belle["BelleFITcent_Dst"][0:10]
            result["BelleFITcov_w"] = slice_matrix(belle["BelleFITcov_Dst"], 0, 10)
        elif key == "BelleFITDstCosL":
            result["BelleFITcent_cosL"] = belle["BelleFITcent_Dst"][10:20]
            result["BelleFITcov_cosL"] = slice_matrix(belle["BelleFITcov_Dst"], 10, 20)
        elif key == "BelleFITDstCosV":
            result["BelleFITcent_cosV"] = belle["BelleFITcent_Dst"][20:30]
            result["BelleFITcov_cosV"] = slice_matrix(belle["BelleFITcov_Dst"], 20, 30)
        elif key == "BelleFITDstChi":
            result["BelleFITcent_chi"] = belle["BelleFITcent_Dst"][30:40]
            result["BelleFITcov_chi"] = slice_matrix(belle["BelleFITcov_Dst"], 30, 40)
        else:
            raise KeyError(f"Unknown data key: {key}")

    selected = set(selected_data)
    result["use_dataset"] = [int(key in selected) for key in DATASET_ORDER]
    return result


def ub_data_decl(ub_model: str) -> str:
    if ub_model == "no-UB":
        return ""
    if ub_model == "soft-UB":
        return "vector[4] UBmax;\nvector<lower=0>[4] UB_sigma;"
    if ub_model == "hard-UB":
        return "vector[4] UBmax;"
    raise ValueError(f"Unknown UBModel: {ub_model}")


def ub_model_term(ub_model: str, fit_functions: FitFunctionSet) -> str:
    if not fit_functions.ub_body:
        return ""
    if ub_model == "no-UB":
        return ""
    if ub_model == "soft-UB":
        return join_lines(["for (k in 1:4) {", "  target += normal_lccdf(UBfunc[k] | UBmax[k], UB_sigma[k]);", "}"])
    if ub_model == "hard-UB":
        return join_lines(["for (k in 1:4) {", "  if (UBfunc[k] > UBmax[k]) {", "    target += negative_infinity();", "  }", "}"])
    raise ValueError(f"Unknown UBModel: {ub_model}")


def build_data_block(selected_data: list[str], np_model: str, ub_model: str) -> str:
    modules = data_modules(np_model)
    decls = [modules[key].data_decl for key in required_data_keys(selected_data)]
    decls.append(f"array[{len(DATASET_ORDER)}] int<lower=0, upper=1> use_dataset;")
    ub_decl = ub_data_decl(ub_model)
    return "data {\n" + join_lines([*decls, ub_decl]) + "\n}"


def cholesky_name(key: str) -> str:
    return f"L_{key}"


def build_transformed_data_block(np_model: str) -> str:
    lines = []
    for key in DATASET_ORDER:
        spec = diagnostic_specs(np_model)[key]
        lines.append(
            f"matrix[{spec.dimension}, {spec.dimension}] {cholesky_name(key)} "
            f"= cholesky_decompose({spec.covariance_name});"
        )
    return "transformed data {\n" + join_lines(lines) + "\n}"


def np_parameter_lines(np_model: str, ff_model: str, ff_nf: int) -> list[str]:
    if np_model == "SM":
        return []
    if np_model == "T":
        return [f"real<lower=-0.5, upper=0.5> C{np_model};", *[f"real {name};" for name in FF_TENSOR_PARAMETER_DECLS.get((ff_model, ff_nf), [])]]
    return [f"real<lower=-0.4, upper=0.3> C{np_model};"]


def build_parameter_block(ff_model: str, ff_nf: int, np_model: str, fit_functions: FitFunctionSet) -> str:
    params = FF_PARAMETER_DECLS.get((ff_model, ff_nf), fit_functions.param_list)
    lines = [*np_parameter_lines(np_model, ff_model, ff_nf), *[f"real {param};" for param in params]]
    return "parameters {\n" + join_lines(lines) + "\n\n}"


def ff_prior_lines(ff_model: str, ff_nf: int) -> list[str]:
    if ff_model == "HQET":
        return HQET_QCDSR_PRIORS[ff_nf]
    sigma = 100 if ff_model == "BSZ" else 1
    return [f"{param} ~ normal(0, {sigma});" for param in FF_PARAMETER_DECLS[(ff_model, ff_nf)]]


def np_prior_lines(np_model: str, ff_model: str, ff_nf: int) -> list[str]:
    if np_model == "SM":
        return []
    lines = [f"C{np_model} ~ normal(0, 0.5);"]
    if np_model == "T":
        sigma = 100 if ff_model == "BSZ" else 1
        lines.extend(f"{param} ~ normal(0, {sigma});" for param in FF_TENSOR_PARAMETER_DECLS.get((ff_model, ff_nf), []))
    return lines


def build_transformed_parameter_block(selected_data: list[str], np_model: str, ub_model: str, fit_functions: FitFunctionSet) -> str:
    modules = data_modules(np_model)
    function_data = [key for key in required_data_keys(selected_data) if key != "BrRatio2025"]
    decls = [modules[key].function_decl for key in function_data]
    if fit_functions.ub_body and ub_model != "no-UB":
        decls.append("vector[4] UBfunc;")

    bodies: list[str] = []
    for key in function_data:
        body_key = modules[key].function_body_key
        if body_key is not None:
            bodies.extend(fit_functions.function_bodies[body_key])
    if fit_functions.ub_body and ub_model != "no-UB":
        bodies.extend(fit_functions.ub_body)

    return "transformed parameters {\n" + join_lines(decls) + "\n" + join_lines(bodies) + "\n}"


def build_model_block(selected_data: list[str], ff_model: str, ff_nf: int, np_model: str, ub_model: str, fit_functions: FitFunctionSet) -> str:
    del selected_data
    terms: list[str] = []
    specs = diagnostic_specs(np_model)
    for key in DATASET_ORDER:
        spec = specs[key]
        index = DATASET_INDEX[key]
        terms.extend(
            [
                f"if (use_dataset[{index}] == 1) {{",
                f"  target += multi_normal_cholesky_lpdf({spec.observed_name} | "
                f"{spec.prediction_name}, {cholesky_name(key)});",
                "}",
            ]
        )
    ub_term = ub_model_term(ub_model, fit_functions)
    return "model {\n" + join_lines([*np_prior_lines(np_model, ff_model, ff_nf), *ff_prior_lines(ff_model, ff_nf), "", *terms, ub_term]) + "\n}"


def resolve_output_keys(selected_data: list[str], generated_outputs: Any, np_model: str) -> list[str]:
    del selected_data, generated_outputs, np_model
    return ["VcbFromBr"]


def output_decl(key: str) -> str:
    if key != "VcbFromBr":
        raise KeyError(key)
    return (
        "vector[2] rawBrRatio_func;\n"
        "vector[2] BrRatio_func;\n"
        "vector[2] BrRatio_data;\n"
        "vector[2] Vcb_Br;\n"
        "real Vcb_mean;"
    )


def output_body(key: str, fit_functions: FitFunctionSet) -> list[str]:
    if key != "VcbFromBr":
        raise KeyError(key)
    return [
        *fit_functions.output_bodies[key],
        "BrRatio_data = multi_normal_rng(BrRatio2025cent, BrRatio2025cov);",
        "Vcb_Br = sqrt(BrRatio_data ./ rawBrRatio_func);",
        "Vcb_mean = mean(sqrt(BrRatio2025cent ./ rawBrRatio_func));",
        "BrRatio_func = rawBrRatio_func * square(Vcb_mean);",
    ]


def diagnostic_output_decls(np_model: str) -> list[str]:
    del np_model
    decls: list[str] = []
    for key in DATASET_ORDER:
        decls.extend(
            [
                f"real T_obs_{key};",
                f"real T_rep_{key};",
                f"real log_lik_{key};",
            ]
        )
    return decls


def diagnostic_output_bodies(np_model: str) -> list[str]:
    bodies: list[str] = []
    specs = diagnostic_specs(np_model)
    for key in DATASET_ORDER:
        spec = specs[key]
        bodies.extend(
            [
                "{",
                f"  vector[{spec.dimension}] y_rep = "
                f"multi_normal_cholesky_rng({spec.prediction_name}, {cholesky_name(key)});",
                f"  vector[{spec.dimension}] z_obs = mdivide_left_tri_low("
                f"{cholesky_name(key)}, {spec.observed_name} - {spec.prediction_name});",
                f"  vector[{spec.dimension}] z_rep = mdivide_left_tri_low("
                f"{cholesky_name(key)}, y_rep - {spec.prediction_name});",
                f"  T_obs_{key} = dot_self(z_obs);",
                f"  T_rep_{key} = dot_self(z_rep);",
                f"  log_lik_{key} = multi_normal_cholesky_lpdf("
                f"{spec.observed_name} | {spec.prediction_name}, {cholesky_name(key)});",
                "}",
            ]
        )
    return bodies


def build_generated_quantities_block(selected_data: list[str], generated_outputs: Any, np_model: str, fit_functions: FitFunctionSet) -> str:
    output_keys = resolve_output_keys(selected_data, generated_outputs, np_model)
    decls = [output_decl(key) for key in output_keys]
    decls.extend(diagnostic_output_decls(np_model))
    bodies: list[str] = []
    for key in output_keys:
        bodies.extend(output_body(key, fit_functions))
    bodies.extend(diagnostic_output_bodies(np_model))
    return "generated quantities {\n" + join_lines(decls) + "\n" + join_lines(bodies) + "\n}"


def build_stan_code(
    *,
    selected_data: list[str],
    generated_outputs: Any,
    ff_model: str,
    ff_nf: int,
    np_model: str,
    ub_model: str,
    fit_functions: FitFunctionSet,
) -> str:
    known = set(DATASET_ORDER)
    unknown = sorted(set(selected_data) - known)
    if unknown:
        raise ValueError(f"Unknown SelectedData key(s): {', '.join(unknown)}")

    blocks = [
        build_data_block(selected_data, np_model, ub_model),
        build_transformed_data_block(np_model),
        build_parameter_block(ff_model, ff_nf, np_model, fit_functions),
        build_transformed_parameter_block(selected_data, np_model, ub_model, fit_functions),
        build_model_block(selected_data, ff_model, ff_nf, np_model, ub_model, fit_functions),
        build_generated_quantities_block(selected_data, generated_outputs, np_model, fit_functions),
    ]
    return "\n\n".join(blocks)


def generate_from_config(config: dict[str, Any]) -> dict[str, Any]:
    working_directory = Path(config.get("WorkingDirectory", Path.cwd()))
    input_directory = Path(config.get("InputDirectory", working_directory / "Stan_inputs"))
    output_root = Path(config.get("OutputRoot", working_directory / "Stan_outputs"))
    analysis_name = config["AnalysisName"]
    ff_model = config.get("FFModel", "BGL")
    ff_nf = int(config.get("FFnf", 2))
    np_model = canonical_np_model(config.get("NPModel", "SM"))
    hqet_model = config.get("HQETModel", "3/2/1")
    cut_control = float(config.get("CutControl", 1.0e-4))
    ub_model = effective_ub_model(ff_model, config.get("UBModel", "no-UB"))
    selected_data = config.get("SelectedData", DEFAULT_SELECTED_DATA)
    generated_outputs = config.get("GeneratedOutputs", "Automatic")
    file_prefix = config.get("DataPointFilePrefix", "standata")
    belle_fit_option = config.get("BelleFitOption", "sc")
    hpqcd23_fit_file = Path(config["HPQCD23FitFile"]) if config.get("HPQCD23FitFile") else None
    cache_directory = Path(config.get("FitFunctionCacheDirectory", input_directory / "StanFunctionCache"))

    prepare_data_point_inputs(
        input_directory=input_directory,
        output_directory=input_directory,
        hpqcd23_fit_file=hpqcd23_fit_file,
        belle_fit_option=belle_fit_option,
        file_prefix=file_prefix,
    )

    fit_functions = load_fit_functions(
        cache_directory=cache_directory,
        ff_model=ff_model,
        ff_nf=ff_nf,
        np_model=np_model,
        hqet_model=hqet_model,
        cut_control=cut_control,
        input_directory=input_directory,
    )

    stan_code = build_stan_code(
        selected_data=selected_data,
        generated_outputs=generated_outputs,
        ff_model=ff_model,
        ff_nf=ff_nf,
        np_model=np_model,
        ub_model=ub_model,
        fit_functions=fit_functions,
    )
    standata = build_standata(selected_data, input_directory, np_model)

    fit_dir = output_root / analysis_name
    data_file = fit_dir / f"standata_{analysis_name}.json"
    stan_file = fit_dir / f"stancode_{analysis_name}.stan"
    param_file = fit_dir / f"stancode_{analysis_name}_param.json"

    write_json(data_file, standata)
    write_text(stan_file, stan_code)
    write_json(param_file, fit_functions.param_list)

    return {
        "AnalysisName": analysis_name,
        "FitDirectory": str(fit_dir),
        "DataFile": str(data_file),
        "StanFile": str(stan_file),
        "ParamFile": str(param_file),
        "SelectedData": selected_data,
        "DatasetOrder": DATASET_ORDER,
        "UseDataset": standata["use_dataset"],
        "GeneratedOutputs": resolve_output_keys(selected_data, generated_outputs, np_model),
        "UBModel": ub_model,
        "Generator": "python-extended",
        "FitFunctionCache": str(cache_directory),
        "IntegratedFunctionSource": str(input_directory / f"IntFunctions_{ff_model}.ma"),
    }


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate the extended master Stan model without Wolfram."
    )
    parser.add_argument("config", type=Path)
    args = parser.parse_args()
    result = generate_from_config(load_json(args.config))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
