from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


DATA_POINT_SYMBOLS = [
    "MILC15data",
    "HPQCD15data",
    "bglMILC21data",
    "bglJLQCD23data",
    "bglHPQCD23data",
    "hqetHPQCD23data",
    "bszLCSR18data",
    "bszLCSR23data",
    "BrRatio2025data",
    "BellFITdata",
]

BELLE_FIT_OPTIONS = {"sc", "re_small", "re_large", "de", "pdg", "rms"}


def _json_load(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _json_dump(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(_jsonable(data), handle, indent=2)
        handle.write("\n")


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _jsonable(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def block_diagonal(mats: list[list[list[float]]]) -> list[list[float]]:
    arrays = [np.asarray(mat, dtype=float) for mat in mats]
    size = sum(mat.shape[0] for mat in arrays)
    result = np.zeros((size, size), dtype=float)
    offset = 0
    for mat in arrays:
        n = mat.shape[0]
        result[offset : offset + n, offset : offset + n] = mat
        offset += n
    return result.tolist()


def block_mask(block_sizes: list[int]) -> np.ndarray:
    size = sum(block_sizes)
    mask = np.zeros((size, size), dtype=float)
    offset = 0
    for block_size in block_sizes:
        mask[offset : offset + block_size, offset : offset + block_size] = 1.0
        offset += block_size
    return mask


def _load_existing_data(input_directory: Path, prefix: str, symbol: str) -> dict[str, Any] | None:
    path = input_directory / f"{prefix}_{symbol}.json"
    if not path.exists():
        return None
    return _json_load(path)


def _load_belle_fit_input(input_directory: Path, dist: str, option: str) -> dict[str, Any]:
    path = input_directory / f"BelleFit_{dist}_inputs.json"
    if not path.exists():
        raise FileNotFoundError(f"BelleFit input JSON was not found: {path}")
    data = _json_load(path)
    key = f"BelleFit_{dist}_{option}"
    cov_key = f"{key}_cov"
    if key not in data or cov_key not in data:
        raise KeyError(f"Required BelleFit keys are missing in {path}: {key}, {cov_key}")
    return {"cent": data[key], "cov": data[cov_key]}


def load_belle_fit_data(input_directory: Path, belle_fit_option: str = "sc") -> dict[str, Any]:
    if belle_fit_option not in BELLE_FIT_OPTIONS:
        available = ", ".join(sorted(BELLE_FIT_OPTIONS))
        raise ValueError(f"Unknown BelleFit option {belle_fit_option!r}. Available options: {available}")

    dst_dists = ["w", "cosL", "cosV", "chi"]
    dst_data = {dist: _load_belle_fit_input(input_directory, dist, belle_fit_option) for dist in dst_dists}
    d_data = _load_belle_fit_input(input_directory, "Dmode", belle_fit_option)

    dst_cent: list[float] = []
    dst_covs: list[list[list[float]]] = []
    for dist in dst_dists:
        dst_cent.extend(dst_data[dist]["cent"])
        dst_covs.append(dst_data[dist]["cov"])

    dst_cov = block_diagonal(dst_covs)
    d_cent = d_data["cent"]
    d_cov = d_data["cov"]

    return {
        "BelleFITcent_Dst": dst_cent,
        "BelleFITcov_Dst": dst_cov,
        "BelleFITcent_D": d_cent,
        "BelleFITcov_D": d_cov,
        "BelleRawFITcent_Dst": dst_cent,
        "BelleRawFITcov_Dst": dst_cov,
        "BelleRawFITcent_D": d_cent,
        "BelleRawFITcov_D": d_cov,
    }


def load_hpqcd23_bgl_fit(fit_file: Path) -> dict[str, Any]:
    if not fit_file.exists():
        raise FileNotFoundError(f"HPQCD23 BGL fit CSV was not found: {fit_file}")

    fit_result = pd.read_csv(fit_file)
    bgl_columns = list(fit_result.columns)[35:70]
    no_tensor_columns = bgl_columns[:20]
    bgl_values = fit_result[bgl_columns]
    no_tensor_values = fit_result[no_tensor_columns]

    full_cov = bgl_values.cov().to_numpy() * block_mask([5, 5, 5, 5, 5, 5, 5])
    no_tensor_cov = no_tensor_values.cov().to_numpy() * block_mask([5, 5, 5, 5])

    return {
        "bglHPQCD23cent": bgl_values.mean().to_list(),
        "bglHPQCD23cov": full_cov.tolist(),
        "bglHPQCD23centNoTensor": no_tensor_values.mean().to_list(),
        "bglHPQCD23covNoTensor": no_tensor_cov.tolist(),
    }


def prepare_data_point_inputs(
    *,
    input_directory: Path,
    output_directory: Path | None = None,
    hpqcd23_fit_file: Path | None = None,
    belle_fit_option: str = "sc",
    file_prefix: str = "standata",
) -> dict[str, Path]:
    output_directory = output_directory or input_directory
    data: dict[str, dict[str, Any]] = {}

    for symbol in DATA_POINT_SYMBOLS:
        existing = _load_existing_data(input_directory, file_prefix, symbol)
        if existing is not None:
            data[symbol] = existing

    if hpqcd23_fit_file is not None and hpqcd23_fit_file.exists():
        data["bglHPQCD23data"] = load_hpqcd23_bgl_fit(hpqcd23_fit_file)
    elif "bglHPQCD23data" not in data:
        raise FileNotFoundError("bglHPQCD23data is missing and HPQCD23FitFile was not available.")

    data["BellFITdata"] = load_belle_fit_data(input_directory, belle_fit_option)

    missing = [symbol for symbol in DATA_POINT_SYMBOLS if symbol not in data]
    if missing:
        raise FileNotFoundError(f"Missing data-point inputs: {', '.join(missing)}")

    exported: dict[str, Path] = {}
    for symbol in DATA_POINT_SYMBOLS:
        path = output_directory / f"{file_prefix}_{symbol}.json"
        _json_dump(path, data[symbol])
        exported[symbol] = path

    _json_dump(output_directory / f"{file_prefix}_available_data.json", list(exported))
    return exported
