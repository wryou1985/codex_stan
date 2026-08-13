"""Run the production HQET 2/1/0 SM soft-UB analysis."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import Stan_FitAnalysis as analysis  # noqa: E402


def main() -> None:
    analysis.ANALYSIS_NAME = "CasePDG-softUB_SM_HQET210"
    analysis.FF_MODEL = "HQET"
    analysis.FF_NF = 1
    analysis.NP_MODEL = "SM"
    analysis.HQET_MODEL = "2/1/0"
    analysis.UB_MODEL = "soft-UB"
    analysis.BELLE_FIT_OPTION = "rms"

    output_file = analysis.run_fit()
    if not output_file.is_file() or output_file.stat().st_size == 0:
        raise RuntimeError(f"Full analysis did not create a non-empty CSV: {output_file}")
    print(f"Full HQET 2/1/0 SM soft-UB output: {output_file}")


if __name__ == "__main__":
    main()

