"""Section 4.4: LRPS factor clusters, k = 10 components and m = 10 clusters (Figs 9-10).

Same data as 04_connectivity.py (all participants; task runs: lag pairs inside the task blocks, pooled).
1. Per participant: right singular vectors V_k of the response Y (k = 10) and the factor-level OLS coefficients
   B = (X'X)^{-1} X' Y V_k, so that the LRPS estimate is A = sum_r b_r v_r'.
2. The pooled vectors are clustered by sign-invariant spherical k-means with one-to-one assignment within
   participant (k = m, so each participant contributes one component to every cluster); clusters are numbered by
   the mean component index of their members.
3. Each cluster's sign-aligned centroid v-bar is the population network; the members' b_r are sign-aligned and
   Procrustes-rotated onto v-bar, averaged (b-bar, the driving pattern) and tested entrywise across participants
   (one-sample t), with Benjamini-Hochberg at FDR 5% over the 116 ROIs within each cluster.
4. Variance shares: the in-sample variance of each network score explained within vs off network
   (lrps.clusters.variance_shares).
5. The clusters of each dataset are matched to the Axcpt (reference) clusters by Hungarian assignment on |cos|
   between centroids.

Outputs (results/clusters/<dataset>/): cluster_membership.csv, cluster_summary_allm.csv,
cluster_roi_coefficients_allm.csv; results/clusters/: drive_varshare_by_subject.csv, matching.csv.

When the analysis finishes, Figs 9-10 are drawn to results/figures (skip with --no-figures).

Example:
    python 05_clusters.py
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

from lrps import RESULTS_DIR
from lrps.figures import make_figures
from lrps.clusters import (cluster_vectors_matched, match_to_reference, renumber_by_component, summarise_clusters,
                           variance_shares)
from lrps.data import lag_pairs, load_dataset
from lrps.models import lrps_factors
from lrps.stats import bh_adjust

DATASETS = ("Axcpt", "Cuedts", "resting")  # first = reference for matching


def run_dataset(ds, args):
    data = load_dataset(ds)
    names = data["roi_names"]
    V, B, pairs = {}, {}, {}
    for sid, y in data["all"].items():
        pairs[sid] = lag_pairs(ds, sid, y)
        V[sid], B[sid] = lrps_factors(*pairs[sid], args.k)
    sids = list(V)
    owners = [(s, r) for s in sids for r in range(args.k)]
    labels, obj = cluster_vectors_matched(V, args.m, n_init=args.n_init, seed=args.seed)
    labels = renumber_by_component(labels, args.k)
    clusters = summarise_clusters(V, B, labels, owners)

    out = args.out / ds
    out.mkdir(parents=True, exist_ok=True)
    memb = pd.DataFrame([{"subject": s, "component": r + 1, "cluster": lab} for (s, r), lab in zip(owners, labels)])
    memb.to_csv(out / "cluster_membership.csv", index=False)
    rows, frames = [], []
    for c in clusters:
        q = bh_adjust(c.p)
        comps = np.array([r + 1 for _, r in c.members])
        rows.append(dict(dataset=ds, cluster=c.cluster, n_subjects=len(c.b_aligned), mean_component=comps.mean(),
                         modal_component=int(np.bincount(comps).argmax()),
                         mean_angle_deg=float(np.degrees(c.angles.mean())), beta_persistence=c.beta,
                         n_sig_BH=int((q < args.alpha).sum()),
                         top_rois_vbar="; ".join(names[j] for j in np.argsort(-np.abs(c.centroid))[:5])))
        frames.append(pd.DataFrame({"cluster": c.cluster, "roi": names, "v_bar": c.centroid, "b_bar": c.b_bar,
                                    "t": c.t, "p": c.p, "p_BH": q, "off_network": c.off_network}))
    pd.DataFrame(rows).to_csv(out / "cluster_summary_allm.csv", index=False)
    pd.concat(frames).to_csv(out / "cluster_roi_coefficients_allm.csv", index=False)

    shares = []
    lab = {(s, r): l for (s, r), l in zip(owners, labels)}
    for sid in sids:
        X, Y = pairs[sid]
        for r in range(args.k):
            sh = variance_shares(X, Y, V[sid], B[sid], r)
            shares.append(dict(dataset=ds, subject=sid, component=r + 1, cluster=int(lab[sid, r]), n_pairs=len(X), **sh))
    print(f"{ds}: {len(sids)} participants, matched clustering mean |cos|^2 = {obj:.3f}", flush=True)
    print(pd.DataFrame(rows).head(args.keep).drop(columns="top_rois_vbar").round(3).to_string(index=False), flush=True)
    return np.stack([c.centroid for c in clusters]), pd.DataFrame(shares)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--m", type=int, default=10)
    ap.add_argument("--n-init", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--keep", type=int, default=4, help="leading clusters reported in the paper")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--no-figures", action="store_true", help="do not draw the figures at the end")
    args = ap.parse_args()
    args.out = RESULTS_DIR / "clusters"
    cents, shares = {}, []
    for ds in DATASETS:
        cents[ds], sh = run_dataset(ds, args)
        shares.append(sh)
    shares = pd.concat(shares)
    shares["R2_total"] = shares.R2_within + shares.R2_off + shares.R2_cross
    shares.to_csv(args.out / "drive_varshare_by_subject.csv", index=False)
    match = []
    for ds in DATASETS:
        for j, (c, sgn) in enumerate(match_to_reference(cents[ds], cents[DATASETS[0]], args.m)):
            match.append(dict(reference_cluster=j + 1, dataset=ds, matched_cluster=c, sign=sgn,
                              abs_cos=float(abs(cents[ds][c - 1] @ cents[DATASETS[0]][j]))))
    match = pd.DataFrame(match)
    match.to_csv(args.out / "matching.csv", index=False)
    print(match[match.reference_cluster <= args.keep].pivot(index="reference_cluster", columns="dataset",
                                                             values="abs_cos").round(2).to_string())
    print(f"wrote {args.out}")
    if not args.no_figures:
        make_figures("05")


if __name__ == "__main__":
    main()
