# codex_stan

Minimal CmdStanPy setup for local and GitHub Actions runs.

## Local setup

```bash
python3.9 -m venv stan_env
source stan_env/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
python -c "import cmdstanpy; cmdstanpy.install_cmdstan(version='2.36.0', cores=2)"
python scripts/smoke_test_cmdstanpy.py
```

If CmdStan is already installed, CmdStanPy should find it under `~/.cmdstan/cmdstan-2.36.0`.

## GitHub Actions

The workflow in `.github/workflows/cmdstanpy-smoke-test.yml` installs Python dependencies,
restores or builds CmdStan 2.36.0, and runs a small Stan smoke test.

The workflow can be started manually from the Actions tab with **Run workflow**.
