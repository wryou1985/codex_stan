"""Generate and compile the HQET 2/1/0 SM soft-UB Stan model."""

from __future__ import annotations

import sys
from pathlib import Path

from cmdstanpy import CmdStanModel


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import Stan_FitAnalysis as analysis  # noqa: E402


def main() -> None:
    analysis.ANALYSIS_NAME = "CompileCheck-softUB_SM_HQET210"
    analysis.FF_MODEL = "HQET"
    analysis.FF_NF = 1
    analysis.NP_MODEL = "SM"
    analysis.HQET_MODEL = "2/1/0"
    analysis.UB_MODEL = "soft-UB"
    analysis.BELLE_FIT_OPTION = "rms"

    generated = analysis.run_python_generator()
    stan_file = Path(generated["StanFile"])
    data_file = Path(generated["DataFile"])
    param_file = Path(generated["ParamFile"])

    stan_data = analysis.load_stan_data(data_file)
    for key in ("UBmax", "UB_sigma"):
        if key not in stan_data:
            raise RuntimeError(f"soft-UB data is missing required key: {key}")

    model = CmdStanModel(stan_file=str(stan_file))
    executable = Path(model.exe_file)
    if not executable.is_file():
        raise RuntimeError(f"CmdStan did not create the executable: {executable}")

    print(f"Stan source: {stan_file}")
    print(f"Stan data: {data_file}")
    print(f"Parameter list: {param_file}")
    print(f"Compiled executable: {executable}")
    print("HQET 2/1/0 SM soft-UB compile check passed.")


if __name__ == "__main__":
    main()

