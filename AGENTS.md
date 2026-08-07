# Codex instructions

This project uses Python and CmdStanPy for Stan model work.

- Use Python 3.9 unless instructed otherwise.
- Install Python dependencies from `requirements.txt`.
- Use CmdStanPy with CmdStan 2.36.0.
- In GitHub Actions, cache `~/.cmdstan` to avoid rebuilding CmdStan on every run.
- Validate setup changes with `python scripts/smoke_test_cmdstanpy.py`.
- Keep generated samples, compiled Stan artifacts, caches, logs, and large outputs out of git.
- Keep changes minimal and show diffs after edits.
