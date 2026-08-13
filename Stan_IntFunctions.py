from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


Monomial = tuple[tuple[str, int], ...]
Polynomial = dict[Monomial, float]

ZERO_MONOMIAL: Monomial = ()

NP_ZERO_VARS = {
    "SM": {"CT", "CVR", "CSR", "CSL"},
    "VR": {"CT", "CSR", "CSL"},
    "SL": {"CT", "CVR", "CSR"},
    "SR": {"CT", "CVR", "CSL"},
    "T": {"CVR", "CSR", "CSL"},
}

NP_MODEL_ALIASES = {
    "SM": "SM",
    "VR": "VR",
    "CVR": "VR",
    "SM+VR": "VR",
    "SL": "SL",
    "CSL": "SL",
    "SM+CSL": "SL",
    "SR": "SR",
    "CSR": "SR",
    "SM+CSR": "SR",
    "T": "T",
    "CT": "T",
    "SM+CT": "T",
}

NP_MODEL_CACHE_TAGS = {
    "SM": "SM",
    "VR": "VR",
    "SL": "CSL",
    "SR": "CSR",
    "T": "CT",
}

HQET_INTFUNCTION_NP_MODELS = {"SM", "VR", "SL", "T"}
SIMPLIFIED_CACHE_VERSION = 2

NF2_ZERO_VARS = {
    "bf02", "bfp2", "bfT2", "bf2", "bF12", "bF22", "bg2", "bT12", "bT22", "bT232",
    "af02", "afp2", "afT2", "aA02", "aA12", "aA122", "aV2", "aT12", "aT22", "aT232",
    "cxi3", "cet2", "cci22", "cci32", "cel11", "cel21", "cel31", "cel41", "cel51", "cel61",
}

BGL_KINETIC_REPLACEMENTS: dict[str, str] = {
    "bf00": "4.937920696927223*bfp0 + 0.3182047238876092*bfp1 + 0.020505441970220025*bfp2 - 0.06444103569456312*bf01 - 0.004152647081387959*bf02",
    "bF10": "0.16745498661163563*bf0",
    "bF20": "21.274993246047366*bF10 + 1.1931515128282635*bF11 + 0.06691473487677196*bF12 - 0.05608234508135209*bF21 - 0.003145229429823857*bF22",
    "bT20": "bT10",
    "bT230": "1.2517780921566837*bT20 - 0.07025793877125706*bT21 + 0.003943333080611104*bT22 + 0.05612651252763973*bT231 - 0.003150185408515299*bT232",
}

BSZ_KINETIC_REPLACEMENTS: dict[str, str] = {
    "af00": "afp0",
    "aA00": "3.562520345*aA120",
    "aT20": "aT10",
    "aA120": "0.2807001513389077*aA10 - 0.01575472056063357*aA11 + 0.0008842575209252033*aA12 + 0.056126512527639734*aA121 - 0.0031501854085153003*aA122",
    "aT230": "1.251787*aT20 - 0.07016627635063376*aT21 + 0.00393803379822834*aT22 + 0.05612651252763974*aT231 - 0.0031501854085153*aT232",
}

HQET_FIXED_REPLACEMENTS: dict[str, str] = {
    "als": "0.0716197243913529",
    "epc": "0.1807",
    "epb": "0.0522",
}


OBS_START_RE = re.compile(r'\s*Obs\["([^"]+)",\s*"([^"]+)"(?:,\s*(\d+))?\]\s*=\s*(.*)')


def canonical_np_model(np_model: str) -> str:
    key = np_model.strip().upper().replace(" ", "")
    try:
        return NP_MODEL_ALIASES[key]
    except KeyError as exc:
        raise ValueError(f"Unknown NPModel: {np_model}") from exc


def np_model_cache_tag(np_model: str) -> str:
    return NP_MODEL_CACHE_TAGS[canonical_np_model(np_model)]


def parse_intfunctions_file(path: Path, required_keys: Iterable[str] | None = None) -> dict[str, str]:
    required = set(required_keys) if required_keys is not None else None
    result: dict[str, str] = {}
    current_key: str | None = None
    current_parts: list[str] = []
    collect_current = False

    def flush_current() -> None:
        nonlocal current_key, current_parts, collect_current
        if current_key is not None and collect_current and current_key not in result:
            result[current_key] = clean_mathematica_expr("".join(current_parts))
        current_key = None
        current_parts = []
        collect_current = False

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            match = OBS_START_RE.match(line)
            if match:
                flush_current()
                if required is not None and required.issubset(result.keys()):
                    break
                dataset, dist, bin_index, expr = match.groups()
                current_key = intfunction_key(dataset, dist, bin_index)
                collect_current = required is None or current_key in required
                current_parts = [expr] if collect_current else []
            elif collect_current:
                current_parts.append(line)
        else:
            flush_current()

    if required is not None:
        missing = sorted(required.difference(result.keys()))
        if missing:
            raise KeyError(f"Missing Obs entries in {path}: {', '.join(missing)}")
    return result


def intfunction_key(dataset: str, dist: str, bin_index: str | int | None = None) -> str:
    dist = "chi" if dist in {"χ", "\\[Chi]"} else dist
    return "|".join(str(part) for part in [dataset, dist, bin_index] if part is not None)


def clean_mathematica_expr(expr: str) -> str:
    expr = re.sub(r"\(\*.*?\*\)", "", expr, flags=re.DOTALL)
    expr = re.sub(r"\s+", "", expr)
    expr = expr.rstrip(";")
    return re.sub(r"(?<=\d)\*\^([+-]?\d+)", r"e\1", expr)


def convert_intfunctions_to_json(
    ma_file: Path,
    json_file: Path,
    required_keys: Iterable[str] | None = None,
    existing_data: dict[str, str] | None = None,
) -> dict[str, str]:
    data = dict(existing_data or {})
    data.update(parse_intfunctions_file(ma_file, required_keys=required_keys))
    json_file.parent.mkdir(parents=True, exist_ok=True)
    with json_file.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)
        handle.write("\n")
    return data


def load_intfunctions(
    input_directory: Path,
    ff_model: str,
    required_keys: Iterable[str] | None = None,
) -> dict[str, str]:
    required = set(required_keys) if required_keys is not None else None
    cache_dir = input_directory / "IntFunctionsCache"
    json_file = cache_dir / f"IntFunctions_{ff_model}.json"
    cached: dict[str, str] = {}
    if json_file.exists():
        with json_file.open("r", encoding="utf-8") as handle:
            cached = json.load(handle)
        if required is None or required.issubset(cached.keys()):
            return cached
    ma_file = input_directory / f"IntFunctions_{ff_model}.ma"
    if not ma_file.exists():
        raise FileNotFoundError(f"Integrated function file was not found: {ma_file}")
    missing = required.difference(cached.keys()) if required is not None else None
    return convert_intfunctions_to_json(
        ma_file,
        json_file,
        required_keys=missing,
        existing_data=cached,
    )


def simplified_intfunctions_cache_path(
    input_directory: Path,
    *,
    ff_model: str,
    ff_nf: int,
    np_model: str,
) -> Path:
    cache_dir = input_directory / "IntFunctionsCache"
    tag = np_model_cache_tag(np_model)
    return cache_dir / f"IntFunctions_{ff_model}_{tag}_nf{ff_nf}.json"


def poly_to_json(poly: Polynomial) -> list[list[Any]]:
    return [
        [coeff, [[var, power] for var, power in monomial]]
        for monomial, coeff in sorted(poly.items())
    ]


def poly_from_json(data: list[list[Any]]) -> Polynomial:
    return {
        tuple((str(var), int(power)) for var, power in monomial): float(coeff)
        for coeff, monomial in data
    }


def load_simplified_intfunction_polys(
    input_directory: Path,
    *,
    ff_model: str,
    ff_nf: int,
    np_model: str,
    required_keys: Iterable[str],
) -> dict[str, Polynomial]:
    canonical_np = canonical_np_model(np_model)
    required = set(required_keys)
    cache_file = simplified_intfunctions_cache_path(
        input_directory,
        ff_model=ff_model,
        ff_nf=ff_nf,
        np_model=canonical_np,
    )

    cached_obs: dict[str, list[list[Any]]] = {}
    if cache_file.exists():
        with cache_file.open("r", encoding="utf-8") as handle:
            cached = json.load(handle)
        if cached.get("cache_version") == SIMPLIFIED_CACHE_VERSION:
            cached_obs = dict(cached.get("obs", {}))
            if required.issubset(cached_obs.keys()):
                return {key: poly_from_json(cached_obs[key]) for key in required}

    if ff_model == "HQET" and canonical_np not in HQET_INTFUNCTION_NP_MODELS:
        allowed = ", ".join(NP_MODEL_CACHE_TAGS[key] for key in sorted(HQET_INTFUNCTION_NP_MODELS))
        raise ValueError(f"HQET IntFunctions case cache is supported for {allowed}; got {np_model}")

    missing = required.difference(cached_obs.keys())
    ma_file = input_directory / f"IntFunctions_{ff_model}.ma"
    if not ma_file.exists():
        raise FileNotFoundError(f"Integrated function file was not found: {ma_file}")

    raw = parse_intfunctions_file(ma_file, required_keys=missing)
    for key, expr in raw.items():
        cached_obs[key] = poly_to_json(
            simplify_int_poly(expr, ff_model=ff_model, ff_nf=ff_nf, np_model=canonical_np)
        )

    cache_file.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "kind": "simplified_polynomial",
        "cache_version": SIMPLIFIED_CACHE_VERSION,
        "source": str(ma_file),
        "ff_model": ff_model,
        "ff_nf": ff_nf,
        "np_model": canonical_np,
        "np_tag": np_model_cache_tag(canonical_np),
        "obs": cached_obs,
    }
    with cache_file.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    return {key: poly_from_json(cached_obs[key]) for key in required}


def parse_polynomial(expr: str) -> Polynomial:
    expr = clean_mathematica_expr(expr)
    if not expr:
        return {}
    if expr[0] not in "+-":
        expr = "+" + expr
    poly: Polynomial = defaultdict(float)
    for sign, body in split_polynomial_terms(expr):
        coeff = -1.0 if sign == "-" else 1.0
        powers: dict[str, int] = defaultdict(int)
        for factor in body.split("*"):
            if not factor:
                continue
            if re.fullmatch(r"(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?", factor):
                coeff *= float(factor)
                continue
            match = re.fullmatch(r"([A-Za-z][A-Za-z0-9]*)(?:\^(\d+))?", factor)
            if not match:
                raise ValueError(f"Unsupported factor in IntFunctions expression: {factor!r}")
            var, power = match.groups()
            powers[var] += int(power or 1)
        monomial = tuple(sorted((var, power) for var, power in powers.items() if power))
        poly[monomial] += coeff
    return {monomial: coeff for monomial, coeff in poly.items() if coeff != 0}


def split_polynomial_terms(expr: str) -> list[tuple[str, str]]:
    terms: list[tuple[str, str]] = []
    sign = expr[0]
    start = 1
    for index in range(1, len(expr)):
        char = expr[index]
        previous = expr[index - 1]
        if char in "+-" and previous not in {"e", "E"}:
            terms.append((sign, expr[start:index]))
            sign = char
            start = index + 1
    terms.append((sign, expr[start:]))
    return terms


def multiply_poly(left: Polynomial, right: Polynomial) -> Polynomial:
    out: Polynomial = defaultdict(float)
    for left_monomial, left_coeff in left.items():
        for right_monomial, right_coeff in right.items():
            powers: dict[str, int] = defaultdict(int)
            for var, power in (*left_monomial, *right_monomial):
                powers[var] += power
            out[tuple(sorted(powers.items()))] += left_coeff * right_coeff
    return {monomial: coeff for monomial, coeff in out.items() if coeff != 0}


def add_poly(left: Polynomial, right: Polynomial) -> Polynomial:
    out: Polynomial = defaultdict(float, left)
    for monomial, coeff in right.items():
        out[monomial] += coeff
    return {monomial: coeff for monomial, coeff in out.items() if abs(coeff) > 0}


def pow_poly(poly: Polynomial, power: int) -> Polynomial:
    out: Polynomial = {ZERO_MONOMIAL: 1.0}
    for _ in range(power):
        out = multiply_poly(out, poly)
    return out


def zero_variables(poly: Polynomial, variables: set[str]) -> Polynomial:
    return {
        monomial: coeff
        for monomial, coeff in poly.items()
        if not any(var in variables for var, _ in monomial)
    }


def replacement_polys(ff_model: str) -> dict[str, Polynomial]:
    replacements = (
        BGL_KINETIC_REPLACEMENTS
        if ff_model == "BGL"
        else BSZ_KINETIC_REPLACEMENTS
        if ff_model == "BSZ"
        else HQET_FIXED_REPLACEMENTS
        if ff_model == "HQET"
        else {}
    )
    return {var: parse_polynomial(expr) for var, expr in replacements.items()}


def apply_replacements(poly: Polynomial, replacements: dict[str, Polynomial]) -> Polynomial:
    changed = True
    out = poly
    while changed:
        changed = False
        next_out: Polynomial = {}
        for monomial, coeff in out.items():
            term: Polynomial = {ZERO_MONOMIAL: coeff}
            replaced_any = False
            for var, power in monomial:
                if var in replacements:
                    term = multiply_poly(term, pow_poly(replacements[var], power))
                    replaced_any = True
                else:
                    term = multiply_poly(term, {((var, power),): 1.0})
            next_out = add_poly(next_out, term)
            changed = changed or replaced_any
        out = next_out
    return out


def simplify_int_expr(expr: str, *, ff_model: str, ff_nf: int, np_model: str, cut_control: float = 1.0e-4) -> str:
    poly = simplify_int_poly(expr, ff_model=ff_model, ff_nf=ff_nf, np_model=np_model)
    return polynomial_to_stan(poly, cut_control=cut_control)


def simplify_int_poly(expr: str, *, ff_model: str, ff_nf: int, np_model: str) -> Polynomial:
    np_model = canonical_np_model(np_model)
    poly = parse_polynomial(expr)
    zero_vars = set(NP_ZERO_VARS[np_model])
    if ff_nf == 1:
        zero_vars.update(NF2_ZERO_VARS)
    poly = zero_variables(poly, zero_vars)
    poly = apply_replacements(poly, replacement_polys(ff_model))
    poly = zero_variables(poly, zero_vars)
    return poly


def scale_poly(poly: Polynomial, factor: float) -> Polynomial:
    return {monomial: coeff * factor for monomial, coeff in poly.items()}


def monomial_to_stan(monomial: Monomial) -> str:
    factors: list[str] = []
    for var, power in monomial:
        factors.extend([var] * power)
    return "*".join(factors)


def balanced_sum(terms: list[str]) -> str:
    if not terms:
        return "0"
    if len(terms) == 1:
        return terms[0]
    midpoint = len(terms) // 2
    return f"({balanced_sum(terms[:midpoint])} + {balanced_sum(terms[midpoint:])})"


def polynomial_to_stan(poly: Polynomial, *, cut_control: float, balanced: bool = False) -> str:
    if not poly:
        return "0"
    max_coeff = max(abs(coeff) for coeff in poly.values())
    threshold = max_coeff * cut_control if max_coeff else 0.0
    terms = [(monomial, coeff) for monomial, coeff in sorted(poly.items()) if abs(coeff) >= threshold]
    if not terms:
        return "0"

    if balanced:
        signed_terms: list[str] = []
        for monomial, coeff in terms:
            body = f"{coeff:.17g}"
            mono = monomial_to_stan(monomial)
            if mono:
                body += "*" + mono
            signed_terms.append(f"({body})" if coeff < 0 else body)
        return balanced_sum(signed_terms)

    pieces: list[str] = []
    for monomial, coeff in terms:
        sign = "-" if coeff < 0 else "+"
        abs_coeff = abs(coeff)
        body = f"{abs_coeff:.17g}"
        mono = monomial_to_stan(monomial)
        if mono:
            body += "*" + mono
        if not pieces:
            pieces.append(("-" if sign == "-" else "") + body)
        else:
            pieces.append(f" {sign} {body}")
    return "".join(pieces)


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert Mathematica IntFunctions_*.ma to Python JSON.")
    parser.add_argument("ma_file", type=Path)
    parser.add_argument("json_file", type=Path)
    args = parser.parse_args()
    data = convert_intfunctions_to_json(args.ma_file, args.json_file)
    print(f"Converted {len(data)} integrated functions to {args.json_file}")


if __name__ == "__main__":
    main()
