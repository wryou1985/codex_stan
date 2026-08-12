from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from Stan_IntFunctions import (
    Polynomial,
    add_poly,
    canonical_np_model,
    load_intfunctions,
    load_simplified_intfunction_polys,
    polynomial_to_stan,
    scale_poly,
    simplify_int_poly,
)


FUNCTION_PREFIXES = {
    "MILC15func": "MILC15",
    "MILC21func": "MILC21",
    "JLQCD23func": "JLQCD23",
    "HPQCD23func": "HPQCD23",
    "LCSR18func": "LCSR18",
    "LCSR23func": "LCSR23",
    "BelleD_func": "BelleFITD",
    "BelleDst_w_func": "BelleFITDstW",
    "BelleDst_cosL_func": "BelleFITDstCosL",
    "BelleDst_cosV_func": "BelleFITDstCosV",
    "BelleDst_chi_func": "BelleFITDstChi",
}

OUTPUT_PREFIXES = {
    "rawBrRatio_func": "VcbFromBr",
    "RD_func": "RD",
    "rawBelleD_func": "RawBelleD",
    "rawBelleDst_w_func": "RawBelleDstW",
    "rawBelleDst_cosL_func": "RawBelleDstCosL",
    "rawBelleDst_cosV_func": "RawBelleDstCosV",
    "rawBelleDst_chi_func": "RawBelleDstChi",
}


@dataclass(frozen=True)
class FitFunctionSet:
    ff_model: str
    ff_nf: int
    np_model: str
    hqet_model: str
    cut_control: float
    param_list: list[str]
    function_bodies: dict[str, list[str]]
    output_bodies: dict[str, list[str]]
    ub_body: list[str]
    source: str

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> "FitFunctionSet":
        return cls(
            ff_model=data["ff_model"],
            ff_nf=int(data["ff_nf"]),
            np_model=data["np_model"],
            hqet_model=data.get("hqet_model", "3/2/1"),
            cut_control=float(data.get("cut_control", 1.0e-4)),
            param_list=list(data["param_list"]),
            function_bodies={key: list(value) for key, value in data["function_bodies"].items()},
            output_bodies={key: list(value) for key, value in data["output_bodies"].items()},
            ub_body=list(data.get("ub_body", [])),
            source=data.get("source", ""),
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "ff_model": self.ff_model,
            "ff_nf": self.ff_nf,
            "np_model": self.np_model,
            "hqet_model": self.hqet_model,
            "cut_control": self.cut_control,
            "param_list": self.param_list,
            "function_bodies": self.function_bodies,
            "output_bodies": self.output_bodies,
            "ub_body": self.ub_body,
            "source": self.source,
        }


def format_cut_control(cut_control: float) -> str:
    return f"{cut_control:.12g}"


def cache_name(ff_model: str, ff_nf: int, np_model: str, hqet_model: str, cut_control: float) -> str:
    hqet_tag = hqet_model.replace("/", "") if ff_model == "HQET" else ""
    parts = [ff_model, f"nf{ff_nf}", np_model]
    if hqet_tag:
        parts.append(f"HQET{hqet_tag}")
    parts.append(f"cut{format_cut_control(cut_control)}")
    return "_".join(parts) + ".json"


def load_fit_functions(
    *,
    cache_directory: Path,
    ff_model: str,
    ff_nf: int,
    np_model: str,
    hqet_model: str,
    cut_control: float,
    input_directory: Path | None = None,
) -> FitFunctionSet:
    np_model = canonical_np_model(np_model)
    path = cache_directory / cache_name(ff_model, ff_nf, np_model, hqet_model, cut_control)
    if ff_model == "HQET" and np_model in {"VR", "SL"} and not path.exists():
        path = cache_directory / cache_name(ff_model, ff_nf, "SM", hqet_model, cut_control)
    if not path.exists():
        if ff_model == "HQET" and np_model == "T":
            raise FileNotFoundError(
                "HQET CT integrated-function cache is supported, but the tensor "
                "fit-function structural cache is not present. Please create "
                f"{path} before running a full HQET CT Stan fit."
            )
        raise FileNotFoundError(
            "Python fit-function cache was not found for "
            f"FFModel={ff_model}, FFnf={ff_nf}, NPModel={np_model}, "
            f"HQETModel={hqet_model}, CutControl={cut_control}: {path}"
        )
    with path.open("r", encoding="utf-8") as handle:
        cache = FitFunctionSet.from_json(json.load(handle))
    if cache.np_model != np_model:
        cache = replace(cache, np_model=np_model)
    if input_directory is None:
        return cache
    return apply_intfunctions(cache, input_directory=input_directory)


def _obs_key(dataset: str, dist: str, bin_index: int | None = None) -> str:
    dist = "chi" if dist in {"χ", "\\[Chi]"} else dist
    parts = [dataset, dist]
    if bin_index is not None:
        parts.append(str(bin_index))
    return "|".join(parts)


def _required_intfunction_keys() -> list[str]:
    keys = [
        _obs_key("BrRatio2025", "D"),
        _obs_key("BrRatio2025", "Dst"),
        _obs_key("BrRatioTau", "D"),
        _obs_key("BrRatioTau", "Dst"),
    ]
    keys.extend(_obs_key("Belle15", "w", bin_index) for bin_index in range(1, 11))
    for dist in ["w", "cosL", "cosV", "chi"]:
        keys.extend(_obs_key("Belle17", dist, bin_index) for bin_index in range(1, 11))
    return keys


def _load_intfunction_polys(cache: FitFunctionSet, input_directory: Path, required_keys: list[str]) -> dict[str, Polynomial]:
    if cache.ff_model == "HQET":
        return load_simplified_intfunction_polys(
            input_directory,
            ff_model=cache.ff_model,
            ff_nf=cache.ff_nf,
            np_model=cache.np_model,
            required_keys=required_keys,
        )

    store = load_intfunctions(input_directory, cache.ff_model, required_keys=required_keys)
    return {
        key: simplify_int_poly(
            store[key],
            ff_model=cache.ff_model,
            ff_nf=cache.ff_nf,
            np_model=cache.np_model,
        )
        for key in required_keys
    }


def _simplify_obs(poly_store: dict[str, Polynomial], cache: FitFunctionSet, dataset: str, dist: str, bin_index: int | None = None, scale: float = 1.0) -> str:
    poly = poly_store[_obs_key(dataset, dist, bin_index)]
    if scale != 1.0:
        poly = scale_poly(poly, scale)
    return polynomial_to_stan(poly, cut_control=cache.cut_control)


def _sum_obs(poly_store: dict[str, Polynomial], cache: FitFunctionSet, items: list[tuple[str, str, int | None, float]]) -> str:
    total = {}
    for dataset, dist, bin_index, scale in items:
        poly = poly_store[_obs_key(dataset, dist, bin_index)]
        if scale != 1.0:
            poly = scale_poly(poly, scale)
        total = add_poly(total, poly)
    return polynomial_to_stan(total, cut_control=cache.cut_control)


def apply_intfunctions(cache: FitFunctionSet, *, input_directory: Path) -> FitFunctionSet:
    if cache.ff_model not in {"BGL", "BSZ", "HQET"}:
        return cache

    required_keys = _required_intfunction_keys()
    poly_store = _load_intfunction_polys(cache, input_directory, required_keys)
    d_raw = [_simplify_obs(poly_store, cache, "Belle15", "w", bin_index, scale=0.06) for bin_index in range(1, 11)]
    d_total = _sum_obs(poly_store, cache, [("Belle15", "w", bin_index, 0.06) for bin_index in range(1, 11)])

    dst_raw = {
        dist: [_simplify_obs(poly_store, cache, "Belle17", dist, bin_index) for bin_index in range(1, 11)]
        for dist in ["w", "cosL", "cosV", "chi"]
    }
    dst_total = _sum_obs(poly_store, cache, [("Belle17", "w", bin_index, 1.0) for bin_index in range(1, 11)])

    br = [
        _simplify_obs(poly_store, cache, "BrRatio2025", "D"),
        _simplify_obs(poly_store, cache, "BrRatio2025", "Dst"),
    ]
    br_tau = [
        _simplify_obs(poly_store, cache, "BrRatioTau", "D"),
        _simplify_obs(poly_store, cache, "BrRatioTau", "Dst"),
    ]

    function_bodies = {key: list(value) for key, value in cache.function_bodies.items()}
    function_bodies["BelleFITD"] = [
        f"BelleD_func[{index}]=({expr})/({d_total});"
        for index, expr in enumerate(d_raw, start=1)
    ]
    function_bodies["BelleFITDstW"] = [
        f"BelleDst_w_func[{index}]=({expr})/({dst_total});"
        for index, expr in enumerate(dst_raw["w"], start=1)
    ]
    function_bodies["BelleFITDstCosL"] = [
        f"BelleDst_cosL_func[{index}]=({expr})/({dst_total});"
        for index, expr in enumerate(dst_raw["cosL"], start=1)
    ]
    function_bodies["BelleFITDstCosV"] = [
        f"BelleDst_cosV_func[{index}]=({expr})/({dst_total});"
        for index, expr in enumerate(dst_raw["cosV"], start=1)
    ]
    function_bodies["BelleFITDstChi"] = [
        f"BelleDst_chi_func[{index}]=({expr})/({dst_total});"
        for index, expr in enumerate(dst_raw["chi"], start=1)
    ]

    output_bodies = {key: list(value) for key, value in cache.output_bodies.items()}
    output_bodies["VcbFromBr"] = [
        f"rawBrRatio_func[{index}]={expr};"
        for index, expr in enumerate(br, start=1)
    ]
    output_bodies["RD"] = [
        f"RD_func[{index}]= ({tau}) / ({light});"
        for index, (tau, light) in enumerate(zip(br_tau, br), start=1)
    ]
    output_bodies["RawBelleD"] = [
        f"rawBelleD_func[{index}]= {expr};"
        for index, expr in enumerate(d_raw, start=1)
    ]
    output_bodies["RawBelleDstW"] = [
        f"rawBelleDst_w_func[{index}]= {expr};"
        for index, expr in enumerate(dst_raw["w"], start=1)
    ]
    output_bodies["RawBelleDstCosL"] = [
        f"rawBelleDst_cosL_func[{index}]= {expr};"
        for index, expr in enumerate(dst_raw["cosL"], start=1)
    ]
    output_bodies["RawBelleDstCosV"] = [
        f"rawBelleDst_cosV_func[{index}]= {expr};"
        for index, expr in enumerate(dst_raw["cosV"], start=1)
    ]
    output_bodies["RawBelleDstChi"] = [
        f"rawBelleDst_chi_func[{index}]= {expr};"
        for index, expr in enumerate(dst_raw["chi"], start=1)
    ]

    return replace(
        cache,
        function_bodies=function_bodies,
        output_bodies=output_bodies,
        source=f"{cache.source}; IntFunctions_{cache.ff_model}",
    )


def _block(text: str, name: str) -> str:
    match = re.search(rf"{re.escape(name)}\s*\{{(.*?)\n\}}", text, flags=re.DOTALL)
    if not match:
        raise ValueError(f"Could not find Stan block: {name}")
    return match.group(1)


def _assignment_lines(block: str, prefixes: dict[str, str]) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {module: [] for module in prefixes.values()}
    for raw_line in block.splitlines():
        line = raw_line.strip()
        for prefix, module in prefixes.items():
            if line.startswith(prefix + "["):
                result[module].append(line)
                break
    return result


def extract_cache_from_stan(
    *,
    stan_file: Path,
    param_file: Path,
    ff_model: str,
    ff_nf: int,
    np_model: str,
    hqet_model: str,
    cut_control: float,
) -> FitFunctionSet:
    stan_text = stan_file.read_text(encoding="utf-8")
    transformed = _block(stan_text, "transformed parameters")
    generated = _block(stan_text, "generated quantities")

    with param_file.open("r", encoding="utf-8") as handle:
        param_list = json.load(handle)

    function_bodies = _assignment_lines(transformed, FUNCTION_PREFIXES)
    output_bodies = _assignment_lines(generated, OUTPUT_PREFIXES)
    ub_body = [line.strip() for line in transformed.splitlines() if line.strip().startswith("UBfunc[")]

    missing_functions = [key for key, value in function_bodies.items() if not value]
    missing_outputs = [key for key, value in output_bodies.items() if not value]
    if missing_functions:
        raise ValueError(f"Missing function assignment lines in {stan_file}: {', '.join(missing_functions)}")
    if missing_outputs:
        raise ValueError(f"Missing generated-output assignment lines in {stan_file}: {', '.join(missing_outputs)}")

    return FitFunctionSet(
        ff_model=ff_model,
        ff_nf=ff_nf,
        np_model=np_model,
        hqet_model=hqet_model,
        cut_control=cut_control,
        param_list=param_list,
        function_bodies=function_bodies,
        output_bodies=output_bodies,
        ub_body=ub_body,
        source=str(stan_file),
    )


def save_fit_function_cache(cache: FitFunctionSet, cache_directory: Path) -> Path:
    cache_directory.mkdir(parents=True, exist_ok=True)
    path = cache_directory / cache_name(
        cache.ff_model,
        cache.ff_nf,
        cache.np_model,
        cache.hqet_model,
        cache.cut_control,
    )
    with path.open("w", encoding="utf-8") as handle:
        json.dump(cache.to_json(), handle, indent=2)
        handle.write("\n")
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description="Create Python fit-function cache from an existing Stan file.")
    parser.add_argument("--stan-file", required=True, type=Path)
    parser.add_argument("--param-file", required=True, type=Path)
    parser.add_argument("--cache-directory", required=True, type=Path)
    parser.add_argument("--ff-model", required=True)
    parser.add_argument("--ff-nf", required=True, type=int)
    parser.add_argument("--np-model", default="SM")
    parser.add_argument("--hqet-model", default="3/2/1")
    parser.add_argument("--cut-control", default=1.0e-4, type=float)
    args = parser.parse_args()

    cache = extract_cache_from_stan(
        stan_file=args.stan_file,
        param_file=args.param_file,
        ff_model=args.ff_model,
        ff_nf=args.ff_nf,
        np_model=args.np_model,
        hqet_model=args.hqet_model,
        cut_control=args.cut_control,
    )
    path = save_fit_function_cache(cache, args.cache_directory)
    print(path)


if __name__ == "__main__":
    main()
