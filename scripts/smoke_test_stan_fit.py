"""Run the real Vcb fit pipeline with deliberately tiny CI settings.

This verifies data/model generation, Stan compilation, pre-fit, fit, and CSV
export.  The resulting samples are only a smoke test and have no physics use.
"""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import Stan_FitAnalysis as analysis  # noqa: E402


def main() -> None:
    analysis.ANALYSIS_NAME = "ci-smoke_SM_BGLnf1"
    # The normalized "sc" covariance is intentionally rank deficient.  Use a
    # positive-definite input option in CI so this test exercises the pipeline
    # without altering the production physics inputs.
    analysis.BELLE_FIT_OPTION = "rms"
    analysis.PRE_FIT_SETTINGS = {
        "chains": 1,
        "iter_warmup": 10,
        "iter_sampling": 10,
        "adapt_delta": 0.8,
        "max_treedepth": 8,
        "seed": 12345,
        "show_progress": False,
    }
    analysis.PRE_FIT_CHAIN_TIMEOUT = 30 * 60
    analysis.PRE_FIT_PARALLEL_CHAINS = 1
    analysis.FIT_SETTINGS = {
        "chains": 1,
        "parallel_chains": 1,
        "iter_warmup": 10,
        "iter_sampling": 10,
        "adapt_delta": 0.8,
        "max_treedepth": 8,
        "seed": 12345,
        "show_progress": False,
    }

    output_file = analysis.run_fit()
    if not output_file.is_file() or output_file.stat().st_size == 0:
        raise RuntimeError(f"Smoke test did not create a non-empty CSV: {output_file}")
    print(f"Smoke test output: {output_file}")


if __name__ == "__main__":
    main()
