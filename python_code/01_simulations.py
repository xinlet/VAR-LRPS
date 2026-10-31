"""Section 3: synthetic experiments (Figs 1-3, S1, S2 and the numbers quoted in Sections 3.2-3.3).

Design
  * VAR(1), p = 30, n = 100, standard normal errors, no intercept; y_1 = 0 and no burn-in (as
    tsDyn::VAR.sim with zero starting values).
  * True A: Experiment 1A flat (A_ii = rho, A_ij = 0.2), 1B flat (A_ij = 0.01), 2 Toeplitz (A_ij = rho^|i-j|),
    rho in {0.1, 0.5, 0.9}; if the largest eigenvalue modulus gamma_max >= 1, A -> A / (gamma_max + 0.1).
  * Estimation error (Fig 2): each method fitted to the whole series; MSE = ||A - A_hat||_F^2 / p^2 for
    LRPS (k = 1..30), the sVAR (lambda = 10^seq(-4, 0, length = 20)) and the OLS VAR; 100 replicates.
  * Rank selection (Fig 3, S1): fixed rolling window of 50 observations; at each t = 51..100 the model is fitted
    to the previous 50 observations, standardised with the window's mean and SD, and predicts y_t on that scale.
    cvMSPE = squared error averaged over the 50 forecasts and 30 variables. Each of 50 tuning series selects the
    k (and lambda) minimising cvMSPE; the consensus k-hat is the most frequent choice (ties -> smaller k;
    lambda ties -> larger lambda) and is evaluated on 100 independent series.
  * S2: the alternative rule k-hat = argmin of the cvMSPE averaged over the tuning series.

Random numbers: SeedSequence(--seed) spawns one child per series, first the 9 x 100 estimation replicates
(experiment, rho, replicate), then for each (experiment, rho) the 50 tuning and then 100 evaluation series.

Outputs (results/simulations): estimation_mse.csv, rank_selection.csv, eval_mspe.csv, cv_curves_tuning.csv,
selection_comparison.csv, summary_estimation.csv, summary_prediction.csv, settings.json.

When the analysis finishes, Figs 1-3, S1 and S2 are drawn to results/figures (skip with --no-figures).

Example:
    python 01_simulations.py --jobs 8
"""

from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd

from lrps import RESULTS_DIR
from lrps.figures import make_figures
from lrps.models import lasso_path_columns, lrps_path, simulate_var
from lrps.simulation import (EXPERIMENTS, K_GRID, LAM_GRID, N, N_EST, N_EVAL, N_TUNE, P, RHOS, WINDOW,  # noqa: F401
                             true_A)
from lrps.utils import consensus, parallel_map


def estimation_task(args):
    exp, rho, rep, seed = args
    A = true_A(exp, rho)
    y = simulate_var(A, N, np.random.default_rng(seed))
    X, Y = y[:-1], y[1:]
    B, V = lrps_path(X, Y)  # Y = X B, so B estimates A'
    At = A.T
    Bl = lasso_path_columns(X, Y, LAM_GRID)
    return dict(exp=exp, rho=rho, rep=rep, var=float(np.mean((At - B) ** 2)),
                lrps=[float(np.mean((At - B @ V[:, :k] @ V[:, :k].T) ** 2)) for k in K_GRID],
                lasso=[float(np.mean((At - b) ** 2)) for b in Bl])


def rolling_errors(y: np.ndarray, window: int = WINDOW):
    """cvMSPE of one series: VAR (scalar), LRPS for every k in K_GRID, sVAR for every lambda in LAM_GRID."""
    n, p = y.shape
    e_var, e_lrps, e_lasso = 0.0, np.zeros(K_GRID.size), np.zeros(LAM_GRID.size)
    for t in range(window, n):  # predict y[t] from a model fitted to y[t-window:t]
        w = y[t - window:t]
        mu, sd = w.mean(axis=0), w.std(axis=0, ddof=1)
        z = (w - mu) / sd
        X, Y = z[:-1], z[1:]
        x_new, y_new = z[-1], (y[t] - mu) / sd
        B, V = lrps_path(X, Y)
        f = x_new @ B
        e_var += np.mean((y_new - f) ** 2)
        cum = np.cumsum((f @ V)[:, None] * V.T, axis=0)  # row r-1: LRPS(k = r) forecast f V_r V_r'
        e_lrps += np.mean((y_new - cum) ** 2, axis=1)
        e_lasso += np.mean((y_new - np.einsum("j,ljm->lm", x_new, lasso_path_columns(X, Y, LAM_GRID))) ** 2, axis=1)
    m = n - window
    return e_var / m, e_lrps / m, e_lasso / m


def prediction_task(args):
    exp, rho, role, rep, seed = args
    y = simulate_var(true_A(exp, rho), N, np.random.default_rng(seed))
    v, lr, la = rolling_errors(y)
    return dict(exp=exp, rho=rho, role=role, rep=rep, var=v, lrps=lr, lasso=la)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--no-figures", action="store_true", help="do not draw the figures at the end")
    args = ap.parse_args()
    args.out = RESULTS_DIR / "simulations"
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    ss = np.random.SeedSequence(args.seed)

    # ---- estimation error (Fig 2) ----
    tasks = [(exp, rho, rep, ss.spawn(1)[0]) for exp in EXPERIMENTS for rho in RHOS for rep in range(N_EST)]
    est = pd.DataFrame(parallel_map(estimation_task, tasks, args.jobs))
    rows = []
    for (exp, rho), g in est.groupby(["exp", "rho"], sort=False):
        lr, la = np.stack(g.lrps.to_numpy()).mean(0), np.stack(g.lasso.to_numpy()).mean(0)
        rows.append(dict(experiment=exp, rho=rho, gamma_max=float(np.max(np.abs(np.linalg.eigvals(true_A(exp, rho))))),
                         mse_var=g["var"].mean(), best_k=int(K_GRID[lr.argmin()]), mse_lrps_best=lr.min(),
                         mse_lrps_k1=lr[0], best_lambda=float(LAM_GRID[la.argmin()]), mse_svar_best=la.min()))
    summ_est = pd.DataFrame(rows)
    est_csv = est.copy()
    for col in ("lrps", "lasso"):
        est_csv[col] = est_csv[col].apply(lambda v: json.dumps([float(x) for x in v]))
    est_csv.to_csv(out / "estimation_mse.csv", index=False)
    summ_est.to_csv(out / "summary_estimation.csv", index=False)
    print(summ_est.to_string(index=False, float_format=lambda v: f"{v:.4g}"), flush=True)

    # ---- rolling-window rank selection and prediction (Figs 3, S1, S2) ----
    tasks = [(exp, rho, role, rep, ss.spawn(1)[0]) for exp in EXPERIMENTS for rho in RHOS
             for role, nrep in (("tune", N_TUNE), ("eval", N_EVAL)) for rep in range(nrep)]
    pred = parallel_map(prediction_task, tasks, args.jobs)
    sel_rows, ev_rows, curve_rows, cmp_rows, summ = [], [], [], [], []
    for exp in EXPERIMENTS:
        for rho in RHOS:
            tune = [r for r in pred if r["exp"] == exp and r["rho"] == rho and r["role"] == "tune"]
            ev = [r for r in pred if r["exp"] == exp and r["rho"] == rho and r["role"] == "eval"]
            ks = [int(K_GRID[np.argmin(r["lrps"])]) for r in tune]
            ls = [float(LAM_GRID[np.argmin(r["lasso"])]) for r in tune]
            kh, lh = consensus(ks, "min"), consensus(ls, "max")
            ki, li = kh - 1, int(np.flatnonzero(LAM_GRID == lh)[0])
            k_avg = int(K_GRID[np.argmin(np.mean([r["lrps"] for r in tune], axis=0))])
            for r, k, lam in zip(tune, ks, ls):
                sel_rows.append(dict(exp=exp, rho=rho, rep=r["rep"], k_selected=k, lam_selected=lam, k_hat=kh,
                                     lam_hat=lh))
                curve_rows += [dict(exp=exp, rho=rho, rep=r["rep"], k=int(k), cvMSPE=float(c))
                               for k, c in zip(K_GRID, r["lrps"])]
            e = pd.DataFrame([dict(exp=exp, rho=rho, rep=r["rep"], mspe_var=r["var"], mspe_lrps=r["lrps"][ki],
                                   mspe_svar=r["lasso"][li], mspe_lrps_average_rule=r["lrps"][k_avg - 1])
                              for r in ev])
            ev_rows.append(e)
            d = e.mspe_lrps - e.mspe_var
            summ.append(dict(experiment=exp, rho=rho, k_hat=kh, k_selected_range=f"{min(ks)}-{max(ks)}", lam_hat=lh,
                             mspe_var_mean=e.mspe_var.mean(), mspe_var_sd=e.mspe_var.std(),
                             mspe_lrps_mean=e.mspe_lrps.mean(), mspe_lrps_sd=e.mspe_lrps.std(),
                             mspe_svar_mean=e.mspe_svar.mean(), mspe_svar_sd=e.mspe_svar.std(),
                             lrps_minus_var_mean=d.mean(), lrps_below_var=int((d < 0).sum())))
            cmp_rows.append(dict(exp=exp, rho=rho, k_consensus=kh, k_average=k_avg,
                                 eval_mspe_consensus=e.mspe_lrps.mean(), eval_mspe_average=e.mspe_lrps_average_rule.mean(),
                                 eval_mspe_var=e.mspe_var.mean()))
    pd.DataFrame(sel_rows).to_csv(out / "rank_selection.csv", index=False)
    pd.concat(ev_rows).to_csv(out / "eval_mspe.csv", index=False)
    pd.DataFrame(curve_rows).to_csv(out / "cv_curves_tuning.csv", index=False)
    pd.DataFrame(cmp_rows).to_csv(out / "selection_comparison.csv", index=False)
    summ_pred = pd.DataFrame(summ)
    summ_pred.to_csv(out / "summary_prediction.csv", index=False)
    with pd.option_context("display.width", 220, "display.max_columns", 30):
        print(summ_pred.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
        print(pd.DataFrame(cmp_rows).to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    (out / "settings.json").write_text(json.dumps(
        {"p": P, "n": N, "window": WINDOW, "rhos": RHOS, "k_grid": [1, P], "lambda_grid": list(LAM_GRID),
         "n_tune": N_TUNE, "n_eval": N_EVAL, "n_estimation_reps": N_EST, "seed": args.seed,
         "lasso": "sklearn.linear_model.lasso_path, per response, fit_intercept=False"}, indent=2))
    print(f"wrote {out}")
    if not args.no_figures:
        make_figures("01")


if __name__ == "__main__":
    main()
