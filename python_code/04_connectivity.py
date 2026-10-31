"""Section 4.4: population connectivity from per-participant VAR(1) estimates (Fig 8, S3 Fig).

Every participant (35 resting, 33 Axcpt, 35 Cuedts) is fitted once on all of their data: the whole run for the
resting state; for the task runs every lag pair inside one of the three task blocks, pooled across blocks.
  VAR (OLS)   no tuning
  LRPS        k = 10 (Fig 8), and k = 5, 20 (S3 Fig)
  sVAR        resting lambda = 0.010 (best lambda at the largest expanding window, t = 159);
              Axcpt 0.010 and Cuedts 0.0075 (lambda-hat with two training blocks, Table 2)
Convention: A[j, i] is the directed edge region j (t - 1) -> region i (t).

Population tests: one-sample t-tests across participants of every entry of A; Benjamini-Hochberg at FDR 5% over the
p(p - 1) = 13,340 directed edges (and separately over the 116 self-connections), and Bonferroni at the 95% level.

Outputs (results/connectivity/<dataset>/): estimates_<method>.npz (per-participant A, A-bar, t, p, q, masks) for
VAR, LRPS (k = 10), sVAR, LRPSk5, LRPSk20; results/connectivity/summary_significance.csv.

When the analysis finishes, Fig 8 and S3 Fig are drawn to results/figures (skip with --no-figures).

Example:
    python 04_connectivity.py --jobs 8
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from lrps import RESULTS_DIR
from lrps.figures import make_figures
from lrps.data import DATASETS, lag_pairs, load_dataset
from lrps.models import fit_ols, fit_svar
from lrps.stats import population_stats
from lrps.utils import parallel_map

LAMBDA = {"resting": 0.01, "Axcpt": 0.01, "Cuedts": 0.0074989420933245585}  # grid values 10^-2, 10^-2.125
K_MAIN, K_EXTRA = 10, (5, 20)


def tuned_lambdas(results: Path) -> dict:
    """The sVAR penalties as selected by 02_block_cv.py (two blocks) and 03_resting_expanding.py, if run."""
    out = {}
    f = results / "resting" / "summary.json"
    if f.exists():
        out["resting"] = json.loads(f.read_text())["lambda_best_last_window"]
    for task in ("Axcpt", "Cuedts"):
        f = results / "block_cv" / "train2" / task / "summary.json"
        if f.exists():
            out[task] = json.loads(f.read_text())["lambda_hat"]
    return out


def _fit(args):
    ds, sid, y, lam = args
    X, Y = lag_pairs(ds, sid, y)
    B = fit_ols(X, Y)
    V = np.linalg.svd(Y, full_matrices=False)[2].T  # right singular vectors of the response
    fits = {"VAR": B, "sVAR": fit_svar(X, Y, lam), "n_pairs": len(Y)}
    for k in (K_MAIN, *K_EXTRA):
        fits["LRPS" if k == K_MAIN else f"LRPSk{k}"] = B @ V[:, :k] @ V[:, :k].T
    return sid, fits


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--lambda-from-results", action="store_true",
                    help="use the sVAR penalties selected by scripts 02 and 03 instead of the paper's values")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--no-figures", action="store_true", help="do not draw the figures at the end")
    args = ap.parse_args()
    args.out = RESULTS_DIR / "connectivity"
    lam = dict(LAMBDA)
    tuned = tuned_lambdas(RESULTS_DIR)
    for ds, v in tuned.items():
        if not np.isclose(v, lam[ds]):
            print(f"note: tuned lambda for {ds} is {v:.4g} (paper: {lam[ds]:.4g})")
    if args.lambda_from_results:
        lam.update(tuned)
    summary = []
    for ds in DATASETS:
        data = load_dataset(ds)
        fits = dict(parallel_map(_fit, [(ds, s, y, lam[ds]) for s, y in data["all"].items()], args.jobs))
        sids = list(fits)
        out = args.out / ds
        out.mkdir(parents=True, exist_ok=True)
        for m in ("VAR", "LRPS", "sVAR") + tuple(f"LRPSk{k}" for k in K_EXTRA):
            As = np.stack([fits[s][m] for s in sids])
            st = population_stats(As, args.alpha)
            np.savez_compressed(out / f"estimates_{m}.npz", subjects=np.array(sids), A=As,
                                roi_names=np.array(data["roi_names"]), **st)
            tb = st["t"][st["bh_edges"]]
            summary.append(dict(dataset=ds, method=m, n_subjects=len(sids), lam=lam[ds] if m == "sVAR" else np.nan,
                                n_pairs_median=int(np.median([fits[s]["n_pairs"] for s in sids])),
                                bh_edges=int(st["bh_edges"].sum()), bh_pos=int((tb > 0).sum()),
                                bh_neg=int((tb < 0).sum()), bonferroni_edges=int(st["bonf_edges"].sum()),
                                bh_self=int(st["bh_self"].sum()), bonferroni_self=int(st["bonf_self"].sum())))
        print(f"{ds}: {len(sids)} participants", flush=True)
    summ = pd.DataFrame(summary)
    args.out.mkdir(parents=True, exist_ok=True)
    summ.to_csv(args.out / "summary_significance.csv", index=False)
    (args.out / "settings.json").write_text(json.dumps({"k_lrps": K_MAIN, "k_extra": K_EXTRA, "lambda": lam,
                                                        "alpha": args.alpha}, indent=2))
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(summ.to_string(index=False))
    if not args.no_figures:
        make_figures("04")


if __name__ == "__main__":
    main()
