from pathlib import Path

from cmdstanpy import CmdStanModel, cmdstan_path


ROOT = Path(__file__).resolve().parents[1]
MODEL_PATH = ROOT / "models" / "normal_mean.stan"


def main() -> None:
    print(f"CmdStan path: {cmdstan_path()}")
    model = CmdStanModel(stan_file=str(MODEL_PATH))
    data = {"N": 8, "y": [0.1, 0.0, -0.2, 0.3, 0.2, -0.1, 0.1, 0.0]}
    fit = model.sample(
        data=data,
        chains=1,
        parallel_chains=1,
        iter_warmup=50,
        iter_sampling=50,
        seed=12345,
        show_progress=False,
    )
    summary = fit.summary()
    print(summary.loc[["mu", "sigma"], ["Mean", "StdDev"]])


if __name__ == "__main__":
    main()
