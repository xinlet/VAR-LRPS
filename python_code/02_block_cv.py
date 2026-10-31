"""Section 4.2: block cross-validation on the DMCC task runs (Tables 2-4, Figs 4-5).

Each DMCC run has three task blocks separated by rest. Splits respect temporal order:
  one training block:  train on block 1 -> test on block 2, train on block 2 -> test on block 3
  two training blocks: train on blocks 1 + 2 -> test on block 3
Lag pairs (y_{t-1}, y_t) are used only when both time points lie inside one block; training pairs are stacked
across the training blocks. Forecasts are one step ahead on the test block.

Tuning (first 10 participants): LRPS k = 1..20 and sVAR lambda on 25 log-spaced values in [1e-3, 1]; k-hat and
lambda-hat minimise the mean over tuning participants of each participant's mean MSE over the splits.
Evaluation (remaining participants): participant RMSE = sqrt(mean over splits of the MSE over test points and
ROIs), MAE likewise; paired LRPS - competitor differences with paired-t 95% CIs.

Outputs (results/block_cv/train<1|2>/<task>/ and results/block_cv/): tuning_split_errors.csv, summary.json,
eval_participant.csv, table_tuning.csv (Table 2), table_performance.csv (Table 3), table_paired.csv (Table 4).

When the analysis finishes, Figs 4-5 are drawn to results/figures (skip with --no-figures).

Example:
    python 02_block_cv.py --jobs 8
"""

from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd

from lrps import RESULTS_DIR
from lrps.figures import make_figures
from lrps.data import block_pairs, load_dataset, subject_blocks
from lrps.models import fit_ols, lasso_path_columns
from lrps.stats import paired
from lrps.utils import parallel_map

TASKS = ("Axcpt", "Cuedts")
K_GRID = list(range(1, 21))
LAM_GRID = [float(v) for v in 10 ** np.linspace(-3, 0, 25)]
LABEL = {"VAR": "VAR (OLS)", "LRPS": "LRPS", "sVAR": "sVAR"}


def forward_splits(n_blocks: int, n_train: int):
    """(training blocks, test block): consecutive training blocks, tested on the next block."""
    return [(tuple(range(s, s + n_train)), s + n_train) for s in range(n_blocks - n_train)]


def split_errors(y: np.ndarray, blocks, k_grid, lam_grid, n_train: int) -> pd.DataFrame:
    """MSE and MAE of every method / parameter value on every split of one participant."""
    rows = []
    for train, test in forward_splits(len(blocks), n_train):
        parts = [block_pairs(y, blocks[i]) for i in train]
        X, Y = np.vstack([a for a, _ in parts]), np.vstack([b for _, b in parts])
        Xt, Yt = block_pairs(y, blocks[test])
        B = fit_ols(X, Y)
        V = np.linalg.svd(Y, full_matrices=False)[2].T  # right singular vectors of the training response
        F = Xt @ B
        preds = {("VAR", np.nan): F}
        preds.update({("LRPS", k): F @ V[:, :k] @ V[:, :k].T for k in k_grid})
        for lam, Bl in zip(lam_grid, lasso_path_columns(X, Y, lam_grid)):
            preds[("sVAR", lam)] = Xt @ Bl
        for (m, par), Pm in preds.items():
            e = Yt - Pm
            rows.append(dict(train_blocks="+".join(str(i + 1) for i in train), test_block=test + 1, method=m,
                             param=par, mse=float(np.mean(e ** 2)), mae=float(np.mean(np.abs(e))),
                             n_train=len(Y), n_test=len(Yt)))
    return pd.DataFrame(rows)


def _worker(args):
    task, sid, y, k_grid, lam_grid, n_train = args
    df = split_errors(y, subject_blocks(task, sid, len(y)), k_grid, lam_grid, n_train)
    df.insert(0, "subject", sid)
    return df


def errors_for(task, subjects, k_grid, lam_grid, n_train, jobs):
    tasks = [(task, s, y, k_grid, lam_grid, n_train) for s, y in subjects.items()]
    return pd.concat(parallel_map(_worker, tasks, jobs), ignore_index=True)


def run(task, n_train, out, jobs):
    data = load_dataset(task)
    tune = errors_for(task, data["tune"], K_GRID, LAM_GRID, n_train, jobs)
    tune.to_csv(out / "tuning_split_errors.csv", index=False)
    sel = {}
    for m in ("LRPS", "sVAR"):
        per_subj = tune[tune.method == m].groupby(["subject", "param"]).mse.mean().unstack()  # subjects x param
        sel[m] = dict(hat=float(per_subj.mean(axis=0).idxmin()), subject=per_subj.idxmin(axis=1),
                      curve=per_subj.mean(axis=0))
    k_hat, lam_hat = int(sel["LRPS"]["hat"]), sel["sVAR"]["hat"]
    print(f"{task}, {n_train} training block(s): k-hat = {k_hat} (per participant "
          f"{int(sel['LRPS']['subject'].min())}-{int(sel['LRPS']['subject'].max())}); lambda-hat = {lam_hat:.4g} "
          f"({sel['sVAR']['subject'].min():.4g}-{sel['sVAR']['subject'].max():.4g})", flush=True)

    ev = errors_for(task, data["eval"], [k_hat], [lam_hat], n_train, jobs)
    part = ev.groupby(["subject", "method"]).agg(mse=("mse", "mean"), mae=("mae", "mean")).reset_index()
    part["rmse"] = np.sqrt(part.mse)
    part.to_csv(out / "eval_participant.csv", index=False)
    perf = part.groupby("method").agg(rmse_mean=("rmse", "mean"), rmse_sd=("rmse", "std"), mae_mean=("mae", "mean"),
                                      mae_sd=("mae", "std")).loc[["VAR", "LRPS", "sVAR"]].reset_index()
    perf.insert(0, "task", task)
    perf.insert(0, "training_blocks", n_train)
    perf["method"] = perf.method.map(LABEL)
    comp = []
    for metric in ("rmse", "mae"):
        w = part.pivot(index="subject", columns="method", values=metric)
        for b in ("VAR", "sVAR"):
            comp.append(dict(training_blocks=n_train, task=task, competitor=LABEL[b], metric=metric.upper(),
                             **paired(w["LRPS"], w[b])))
    (out / "summary.json").write_text(json.dumps(
        {"task": task, "training_blocks": n_train, "k_hat": k_hat, "lambda_hat": lam_hat, "k_grid": K_GRID,
         "lambda_grid": LAM_GRID, "lrps_subject_choices": sel["LRPS"]["subject"].astype(int).to_dict(),
         "svar_subject_choices": sel["sVAR"]["subject"].to_dict(),
         "mean_cv_mse_lrps": {str(int(k)): v for k, v in sel["LRPS"]["curve"].items()},
         "mean_cv_mse_svar": {str(k): v for k, v in sel["sVAR"]["curve"].items()}}, indent=2))
    tun = dict(training_blocks=n_train, task=task, k_hat=k_hat,
               k_range=f"{int(sel['LRPS']['subject'].min())}-{int(sel['LRPS']['subject'].max())}", lambda_hat=lam_hat,
               lambda_range=f"{sel['sVAR']['subject'].min():.4f}-{sel['sVAR']['subject'].max():.4f}")
    return tun, perf, comp


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--no-figures", action="store_true", help="do not draw the figures at the end")
    args = ap.parse_args()
    args.out = RESULTS_DIR / "block_cv"
    tun, perf, comp = [], [], []
    for n_train in (1, 2):
        for task in TASKS:
            out = args.out / f"train{n_train}" / task
            out.mkdir(parents=True, exist_ok=True)
            t, p, c = run(task, n_train, out, args.jobs)
            tun.append(t)
            perf.append(p)
            comp += c
    tun, perf, comp = pd.DataFrame(tun), pd.concat(perf), pd.DataFrame(comp)
    tun.to_csv(args.out / "table_tuning.csv", index=False)
    perf.to_csv(args.out / "table_performance.csv", index=False)
    comp.to_csv(args.out / "table_paired.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print("\nTable 2\n", tun.to_string(index=False))
        print("\nTable 3\n", perf.round(3).to_string(index=False))
        print("\nTable 4\n", comp[comp.metric == "RMSE"].round(3).to_string(index=False))
    print(f"wrote {args.out}")
    if not args.no_figures:
        make_figures("02")


if __name__ == "__main__":
    main()
