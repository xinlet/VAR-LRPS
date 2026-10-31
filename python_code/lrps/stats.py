"""Paired comparisons, population tests and multiple-testing corrections."""

from __future__ import annotations

import numpy as np
from scipy import stats


def paired(a, b) -> dict:
    """Paired comparison of per-participant errors: mean of a - b with a paired-t 95% CI, the number of
    participants with a < b, and the two-sided paired-t p-value."""
    d = np.asarray(a, float) - np.asarray(b, float)
    d = d[~np.isnan(d)]
    m, se = d.mean(), d.std(ddof=1) / np.sqrt(d.size)
    q = stats.t.ppf(0.975, d.size - 1)
    t = m / se
    return dict(n=d.size, mean_diff=m, ci_low=m - q * se, ci_high=m + q * se, n_lower=int((d < 0).sum()),
                p=2 * stats.t.sf(abs(t), d.size - 1))


def mean_ci(x) -> tuple[float, float, float]:
    x = np.asarray(x, float)
    m = x.mean()
    h = stats.t.ppf(0.975, x.size - 1) * x.std(ddof=1) / np.sqrt(x.size)
    return m, m - h, m + h


def bh_adjust(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (q-values), any shape."""
    flat = np.asarray(p, dtype=float).ravel()
    n = flat.size
    order = np.argsort(flat)
    adj = np.empty(n)
    adj[order] = np.minimum.accumulate((flat[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(adj, 1.0).reshape(np.shape(p))


def population_stats(As: np.ndarray, alpha: float = 0.05) -> dict:
    """Across-participant tests of every entry of the VAR coefficient matrices As (m, p, p).

    A-bar (mean), one-sample t statistics, two-sided p-values (t_{m-1}); Bonferroni at level alpha and
    Benjamini-Hochberg at FDR alpha, applied separately to the p(p-1) directed edges and to the p
    self-connections. Entries that are zero for every participant (e.g. sparse VAR) have no t statistic and
    enter the BH step with p = 1."""
    m, p, _ = As.shape
    Abar = As.mean(axis=0)
    sd = As.std(axis=0, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = Abar / (sd / np.sqrt(m))
    t[sd == 0] = np.nan
    pv = 2 * stats.t.sf(np.abs(t), m - 1)
    off = ~np.eye(p, dtype=bool)
    p1 = np.where(np.isfinite(pv), pv, 1.0)
    q = np.full_like(pv, np.nan)
    q[off], q[~off] = bh_adjust(p1[off]), bh_adjust(p1[~off])
    return dict(Abar=Abar, t=t, p=pv, q=q,
                bonf_edges=off & (pv < alpha / off.sum()), bonf_self=~off & (pv < alpha / p),
                bh_edges=off & (q <= alpha), bh_self=~off & (q <= alpha))
