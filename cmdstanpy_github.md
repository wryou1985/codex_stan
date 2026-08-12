# CmdStanPy GitHub Actions Setup

This note summarizes the CmdStanPy setup used in this repository so it can be reused in another GitHub project.

## Goal

Run a small CmdStanPy smoke test on GitHub Actions, cache CmdStan between runs, and save test CSV outputs as downloadable GitHub Actions artifacts.

## Files Used

```text
.
├── .github/workflows/cmdstanpy-smoke-test.yml
├── .gitignore
├── models/
│   └── normal_mean.stan
├── requirements.txt
└── scripts/
    └── smoke_test_cmdstanpy.py
```

## Python Dependencies

`requirements.txt` pins the Python-side environment:

```text
cmdstanpy==1.2.5
numpy==1.26.4
pandas==2.2.3
scipy==1.10.1
arviz==0.17.1
matplotlib==3.9.4
xarray==2024.7.0
```

For another project, start with these versions, then relax or update them only after confirming the workflow still passes.

## GitHub Actions Workflow

Place this file at:

```text
.github/workflows/cmdstanpy-smoke-test.yml
```

Current workflow:

```yaml
name: CmdStanPy smoke test

"on":
  workflow_dispatch:
  push:
    branches:
      - main

jobs:
  smoke-test:
    runs-on: ubuntu-latest
    timeout-minutes: 120

    env:
      CMDSTAN_VERSION: "2.36.0"
      MPLCONFIGDIR: "/tmp/matplotlib"

    steps:
      - uses: actions/checkout@v4

      - uses: actions/setup-python@v5
        with:
          python-version: "3.9"
          cache: "pip"

      - name: Install system dependencies
        run: |
          sudo apt-get update
          sudo apt-get install -y build-essential make

      - name: Install Python dependencies
        run: |
          python -m pip install --upgrade pip
          pip install -r requirements.txt

      - name: Cache CmdStan
        id: cache-cmdstan
        uses: actions/cache@v4
        with:
          path: ~/.cmdstan
          key: ${{ runner.os }}-cmdstan-${{ env.CMDSTAN_VERSION }}

      - name: Install CmdStan if missing
        if: steps.cache-cmdstan.outputs.cache-hit != 'true'
        run: |
          python - <<'PY'
          import os
          import cmdstanpy

          cmdstanpy.install_cmdstan(
              version=os.environ["CMDSTAN_VERSION"],
              cores=2,
              progress=True,
          )
          print("CmdStan path:", cmdstanpy.cmdstan_path())
          PY

      - name: Check CmdStan
        run: |
          python - <<'PY'
          import cmdstanpy

          print("CmdStan path:", cmdstanpy.cmdstan_path())
          PY

      - name: Run smoke test
        run: |
          python scripts/smoke_test_cmdstanpy.py

      - name: Upload smoke test CSV outputs
        if: always()
        uses: actions/upload-artifact@v4
        with:
          name: cmdstanpy-smoke-test-csv
          path: results/*.csv
          if-no-files-found: warn
```

## Important Workflow Points

- `workflow_dispatch` lets you run the workflow manually from the GitHub Actions tab.
- `push` runs the workflow automatically when changes are pushed to `main`.
- `CMDSTAN_VERSION` controls the CmdStan version installed on GitHub Actions.
- `actions/cache@v4` caches `~/.cmdstan`, so CmdStan is not rebuilt every run.
- The cache key includes `CMDSTAN_VERSION`; changing the version creates a fresh cache.
- `MPLCONFIGDIR` is set to `/tmp/matplotlib` to avoid matplotlib cache warnings in CI.
- `actions/upload-artifact@v4` saves CSV outputs from `results/*.csv`.

The first run may take several minutes because CmdStan is downloaded and built. Later runs should be faster if the CmdStan cache is restored.

## Smoke Test Script

Place the Python test script at:

```text
scripts/smoke_test_cmdstanpy.py
```

Current script:

```python
from pathlib import Path

from cmdstanpy import CmdStanModel, cmdstan_path


ROOT = Path(__file__).resolve().parents[1]
MODEL_PATH = ROOT / "models" / "normal_mean.stan"
RESULTS_DIR = ROOT / "results"
SUMMARY_PATH = RESULTS_DIR / "smoke_test_summary.csv"
DRAWS_PATH = RESULTS_DIR / "smoke_test_draws.csv"


def main() -> None:
    print(f"CmdStan path: {cmdstan_path()}")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    model = CmdStanModel(stan_file=str(MODEL_PATH))
    data = {"N": 8, "y": [0.1, 0.0, -0.2, 0.3, 0.2, -0.1, 0.1, 0.0]}
    fit = model.sample(
        data=data,
        chains=1,
        parallel_chains=1,
        iter_warmup=50,
        iter_sampling=50,
        seed=12345,
        output_dir=str(RESULTS_DIR),
        show_progress=False,
    )
    summary = fit.summary()
    selected_summary = summary.loc[["mu", "sigma"], ["Mean", "StdDev"]]
    selected_summary.to_csv(SUMMARY_PATH)
    fit.draws_pd().to_csv(DRAWS_PATH, index=False)
    print(selected_summary)
    print(f"Wrote summary CSV: {SUMMARY_PATH}")
    print(f"Wrote draws CSV: {DRAWS_PATH}")


if __name__ == "__main__":
    main()
```

This creates:

```text
results/smoke_test_summary.csv
results/smoke_test_draws.csv
results/normal_mean-*.csv
```

`normal_mean-*.csv` is the raw CmdStan output. `smoke_test_draws.csv` is a pandas-friendly CSV created from `fit.draws_pd()`.

## Example Stan Model

The sample model is:

```text
models/normal_mean.stan
```

```stan
data {
  int<lower=1> N;
  vector[N] y;
}

parameters {
  real mu;
  real<lower=0> sigma;
}

model {
  mu ~ normal(0, 5);
  sigma ~ exponential(1);
  y ~ normal(mu, sigma);
}
```

For another project, replace this file and the `data` object in `smoke_test_cmdstanpy.py` with the actual model and test data.

## Git Ignore

Use a `.gitignore` that keeps source files but ignores generated Stan outputs:

```text
__pycache__/
*.py[cod]
.venv/
stan_env/
.cmdstan/
.pytest_cache/
.mypy_cache/
.DS_Store

results/
*.csv
*.draws
*.hpp
*.o
*.d
*.exe
*.log

models/*
!models/
!models/*.stan
```

This keeps `.stan` files under version control while avoiding compiled binaries, CSV draws, logs, and local virtual environments.

## How To Run Locally

If you already have a working local environment:

```bash
source stan_env/bin/activate
MPLCONFIGDIR=/tmp/matplotlib python scripts/smoke_test_cmdstanpy.py
```

For a new local environment:

```bash
python3.9 -m venv stan_env
source stan_env/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
python -c "import cmdstanpy; cmdstanpy.install_cmdstan(version='2.36.0', cores=2)"
MPLCONFIGDIR=/tmp/matplotlib python scripts/smoke_test_cmdstanpy.py
```

## How To Run On GitHub

1. Push the files to GitHub.
2. Open the repository on GitHub.
3. Open the Actions tab.
4. Select `CmdStanPy smoke test`.
5. Click `Run workflow`.

The workflow also runs automatically when pushing to `main`.

## How To Download CSV Outputs

After the workflow finishes:

1. Open the completed workflow run.
2. Scroll to the Artifacts section.
3. Download `cmdstanpy-smoke-test-csv`.

The artifact contains CSV files from `results/*.csv`. These files are not committed back into the repository.

## Porting Checklist For Another Project

1. Copy `requirements.txt`.
2. Copy `.github/workflows/cmdstanpy-smoke-test.yml`.
3. Copy or adapt `scripts/smoke_test_cmdstanpy.py`.
4. Put Stan models under `models/`.
5. Update `MODEL_PATH` and test data in the smoke test script.
6. Keep generated outputs in `results/`.
7. Add the `.gitignore` rules above.
8. Run once manually from GitHub Actions.
9. Confirm the CmdStan cache is used on the second run.
10. Download the CSV artifact and verify the output format.

## Notes

- GitHub Actions runs on Linux, not macOS, so avoid hard-coded local paths such as `/Users/...`.
- Use repository-relative paths through `Path(__file__).resolve().parents[1]`.
- Keep `results/` ignored by Git. Use artifacts for outputs you want to download.
- If CmdStan needs to be rebuilt, change `CMDSTAN_VERSION` or clear the GitHub Actions cache.
- If a real model takes a long time, increase `timeout-minutes`.
