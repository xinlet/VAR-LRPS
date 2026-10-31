"""Section 4.3: expanding-window cross-validation on the resting-state data (Figs 6-7).

For every participant the models are refitted at each forecast origin t = t0, ..., n - 1 (default t0 = 120) on all
earlier observations y[0:t], and forecast y[t] one step ahead from y[t-1] (40 forecasts per participant).

Tuning (first 10 participants): the one-step error is recorded at every origin for VAR (OLS), LRPS k = 1..20 and
sVAR on 25 log-spaced lambda in [1e-3, 1]. Each participant's CV RMSE for a parameter value is
sqrt(mean over origins of the MSE over ROIs); k-hat and lambda-hat minimise the average over tuning participants.
Evaluation (remaining 25 participants): the same scheme with k-hat and lambda-hat; participant RMSE and MAE and
paired LRPS - competitor differences with paired-t 95% CIs.

Outputs (results/resting): tuning_errors_by_origin.csv.gz, eval_errors_by_origin.csv.gz, eval_participant.csv,
table_performance.csv, table_paired.csv, summary.json (includes the best lambda at the largest window, used for
the connectivity analysis).

When the analysis finishes, Figs 6-7 are drawn to results/figures (skip with --no-figures).

Example:
    python 03_resting_expanding.py --jobs 8
"""

from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd

from lrps import RESULTS_DIR
from lrps.figures import make_figures
from lrps.data import load_dataset
from lrps.models import lasso_path_columns, lrps_path
from lrps.stats import paired
from lrps.utils import parallel_map

LABEL = {"VAR": "VAR (OLS)", "LRPS": "LRPS", "sVAR": "sVAR"}


def origin_errors(y: np.ndarray, t: int, k_grid, lam_grid) -> list[dict]:
    """One-step errors at origin t (model fitted to y[0:t]) for every method and parameter value."""
    X, Y = y[:t - 1], y[1:t]
    x_new, y_new = y[t - 1], y[t]
    B, V = lrps_path(X, Y)
    f = x_new @ B
    preds = [("VAR", np.nan, f)] + [("LRPS", k, f @ V[:, :k] @ V[:, :k].T) for k in k_grid]
    if len(lam_grid):
        preds += [("sVAR", lam, x_new @ Bl) for lam, Bl in zip(lam_grid, lasso_path_columns(X, Y, lam_grid))]
    return [dict(method=m, param=par, train_size=t, mse=float(np.mean((y_new - pr) ** 2)),
                 mae=float(np.mean(np.abs(y_new - pr)))) for m, par, pr in preds]


def _worker(args):
    sid, y, t, k_grid, lam_grid = args
    return [dict(subject=sid, **r) for r in origin_errors(y, t, k_grid, lam_grid)]


def errors_for(subjects, t0, k_grid, lam_grid, jobs) -> pd.DataFrame:
    tasks = [(s, y, t, k_grid, lam_grid) for s, y in subjects.items() for t in range(t0, len(y))]
    df = pd.DataFrame([r for rows in parallel_map(_worker, tasks, jobs) for r in rows])
    df["rmse"] = np.sqrt(df.mse)
    return df


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--t0", type=int, default=120, help="first forecast origin (initial training size)")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--no-figures", action="store_true", help="do not draw the figures at the end")
    args = ap.parse_args()
    args.out = RESULTS_DIR / "resting"
    args.out.mkdir(parents=True, exist_ok=True)
    k_grid = list(range(1, 21))
    lam_grid = [float(v) for v in 10 ** np.linspace(-3, 0, 25)]
    data = load_dataset("resting")

    tune = errors_for(data["tune"], args.t0, k_grid, lam_grid, args.jobs)
    tune.to_csv(args.out / "tuning_errors_by_origin.csv.gz", index=False)
    subj_rmse = tune.groupby(["method", "param", "subject"], dropna=False).mse.mean().pow(0.5)
    mean_rmse = subj_rmse.groupby(["method", "param"], dropna=False).mean()
    k_hat, lam_hat = int(mean_rmse["LRPS"].idxmin()), float(mean_rmse["sVAR"].idxmin())
    choice = {m: subj_rmse[m].unstack("subject").idxmin() for m in ("LRPS", "sVAR")}
    last = tune[(tune.method == "sVAR") & (tune.train_size == tune.train_size.max())]
    lam_last = float(last.groupby("param").rmse.mean().idxmin())  # best lambda at the largest window
    print(f"resting: k-hat = {k_hat} (per participant {sorted(choice['LRPS'].astype(int))}); lambda-hat = "
          f"{lam_hat:.4g} (per participant {np.round(sorted(choice['sVAR']), 4)}); best lambda at t = "
          f"{tune.train_size.max()}: {lam_last:.4g}", flush=True)

    ev = errors_for(data["eval"], args.t0, [k_hat], [lam_hat], args.jobs)
    ev.to_csv(args.out / "eval_errors_by_origin.csv.gz", index=False)
    part = ev.groupby(["subject", "method"]).agg(mse=("mse", "mean"), mae=("mae", "mean")).reset_index()
    part["rmse"] = np.sqrt(part.mse)
    part.to_csv(args.out / "eval_participant.csv", index=False)
    perf = part.groupby("method").agg(rmse_mean=("rmse", "mean"), rmse_sd=("rmse", "std"), mae_mean=("mae", "mean"),
                                      mae_sd=("mae", "std")).loc[["VAR", "LRPS", "sVAR"]]
    perf.index = [LABEL[m] for m in perf.index]
    comp = []
    for metric in ("rmse", "mae"):
        w = part.pivot(index="subject", columns="method", values=metric)
        comp += [dict(competitor=LABEL[b], metric=metric.upper(), **paired(w["LRPS"], w[b])) for b in ("VAR", "sVAR")]
    comp = pd.DataFrame(comp)
    perf.to_csv(args.out / "table_performance.csv")
    comp.to_csv(args.out / "table_paired.csv", index=False)
    (args.out / "summary.json").write_text(json.dumps(
        {"t0": args.t0, "k_grid": k_grid, "lambda_grid": lam_grid, "k_hat": k_hat, "lambda_hat": lam_hat,
         "lambda_best_last_window": lam_last, "lrps_subject_choices": choice["LRPS"].astype(int).to_dict(),
         "svar_subject_choices": choice["sVAR"].to_dict(),
         "mean_cv_rmse": {m: {str(k): v for k, v in mean_rmse[m].items()} for m in ("LRPS", "sVAR", "VAR")}},
        indent=2, default=float))
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(perf.round(4).to_string(), "\n")
        print(comp.round(4).to_string(index=False))
    print(f"wrote {args.out}")
    if not args.no_figures:
        make_figures("03")


if __name__ == "__main__":
    main()
