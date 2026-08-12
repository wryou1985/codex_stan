"""Run the production BGL nf=2 SM analysis on GitHub Actions."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import Stan_FitAnalysis as analysis  # noqa: E402


def main() -> None:
    analysis.ANALYSIS_NAME = "CasePDG-noUB_SM_BGLnf2"
    analysis.FF_MODEL = "BGL"
    analysis.FF_NF = 2
    analysis.NP_MODEL = "SM"
    analysis.HQET_MODEL = "3/2/1"
    analysis.UB_MODEL = "no-UB"
    analysis.BELLE_FIT_OPTION = "rms"

    output_file = analysis.run_fit()
    if not output_file.is_file() or output_file.stat().st_size == 0:
        raise RuntimeError(f"Full analysis did not create a non-empty CSV: {output_file}")
    print(f"Full BGL nf=2 SM output: {output_file}")


if __name__ == "__main__":
    main()

