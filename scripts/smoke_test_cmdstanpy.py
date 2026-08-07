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
