"""Run the production BSZ nf=2 SM extended all-dataset LODO analysis."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import Stan_FitAnalysis_extended as analysis  # noqa: E402


def main() -> None:
    analysis.ANALYSIS_NAME = "BSZ_nf2_SM_extended"
    analysis.FF_MODEL = "BSZ"
    analysis.FF_NF = 2
    analysis.NP_MODEL = "SM"
    analysis.HQET_MODEL = "3/2/1"
    analysis.UB_MODEL = "no-UB"
    analysis.BELLE_FIT_OPTION = "rms"

    analysis.RUN_ALL_DATA = True
    analysis.RUN_LODO = True
    analysis.LODO_DATASETS = list(analysis.DATASET_ORDER)
    analysis.RUN_NODE_SPLIT = False
    analysis.NODE_SPLIT_DATASETS = []
    analysis.RUN_COVARIANCE_SCAN = False
    analysis.COVARIANCE_SCAN_DATASETS = []

    analysis.PRE_FIT_CHAIN_TIMEOUT = 200.0

    outputs = analysis.run_extended_analysis()
    metadata = outputs.get("metadata")
    if metadata is None or not metadata.is_file() or metadata.stat().st_size == 0:
        raise RuntimeError("Extended LODO analysis did not create run_metadata.json.")
    print(f"Full BSZ nf=2 SM extended LODO metadata: {metadata}")


if __name__ == "__main__":
    main()
