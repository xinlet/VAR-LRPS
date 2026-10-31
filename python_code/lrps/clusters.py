"""Population summaries of LRPS in the factor space (Section 4.4).

For participant l with the rank-k right singular vectors V_k of the response Y, the factor-level coefficients
are B = (X'X)^{-1} X' Y V_k, so that the LRPS estimate decomposes as A = sum_r b_r v_r'. The pooled singular
vectors are clustered (sign-invariant spherical k-means, one component per participant per cluster); each
cluster gets a sign-aligned centroid v-bar (a population network), and each member's b_r is sign-aligned and
Procrustes-rotated onto v-bar before averaging and testing across participants.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import stats
from scipy.optimize import linear_sum_assignment


def cluster_vectors_matched(V: dict, m: int, n_init: int = 50, max_iter: int = 200, seed: int = 0):
    """Sign-invariant spherical k-means with one-to-one assignment within participant.

    V: participant -> (p x k) orthonormal loadings. Each participant's k vectors are assigned to distinct
    clusters (Hungarian matching on |cos|^2; requires k <= m). Centroids are the leading eigenvectors of the
    within-cluster scatter. Returns (labels 1..m in the order [(s, r) for s in V for r in range(k)],
    objective = mean |cos|^2 to the centroid); labels are provisional (see renumber_by_component)."""
    sids = list(V)
    k = V[sids[0]].shape[1]
    if k > m:
        raise ValueError("matched clustering needs k <= m")
    rng = np.random.default_rng(seed)
    best = (-np.inf, None)
    for _ in range(n_init):
        U = V[sids[rng.integers(len(sids))]].copy()  # init: one participant's own vectors (+ random extras)
        if m > k:
            U = np.linalg.qr(np.c_[U, rng.standard_normal((U.shape[0], m - k))])[0]
        assign = None
        for _ in range(max_iter):
            new = {s: linear_sum_assignment(-(V[s].T @ U) ** 2)[1] for s in sids}  # cluster of each component
            if assign is not None and all(np.array_equal(new[s], assign[s]) for s in sids):
                break
            assign = new
            for c in range(m):
                vs = [V[s][:, r] for s in sids for r in range(k) if assign[s][r] == c]
                if vs:
                    vs = np.stack(vs)
                    U[:, c] = np.linalg.eigh(vs.T @ vs)[1][:, -1]
        obj = np.mean([(V[s][:, r] @ U[:, assign[s][r]]) ** 2 for s in sids for r in range(k)])
        if obj > best[0]:
            best = (obj, assign)
    obj, assign = best
    labels = np.array([assign[s][r] + 1 for s in sids for r in range(k)])
    # provisional order: decreasing size, ties by mean component index
    keys = sorted(np.unique(labels), key=lambda c: (-np.sum(labels == c),
                  np.mean([r for (i, r) in enumerate(np.tile(np.arange(k), len(sids))) if labels[i] == c])))
    relabel = {old: new for new, old in enumerate(keys, start=1)}
    return np.array([relabel[c] for c in labels]), float(obj)


def renumber_by_component(labels: np.ndarray, k: int) -> np.ndarray:
    """Clusters 1..m ordered by the mean component index of their members (cluster 1 = leading components)."""
    comp = np.tile(np.arange(k), len(labels) // k)
    keys = sorted(np.unique(labels), key=lambda c: comp[labels == c].mean())
    relabel = {old: new for new, old in enumerate(keys, start=1)}
    return np.array([relabel[c] for c in labels])


def procrustes_to(v: np.ndarray, target: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Apply to b the rotation (closest to the identity) that maps unit v onto unit target within
    span{v, target}; v must satisfy v'target >= 0."""
    c = float(np.clip(v @ target, -1.0, 1.0))
    w = v - c * target
    nw = np.linalg.norm(w)
    if nw < 1e-12:
        return b.copy()
    w /= nw
    s = np.sqrt(max(0.0, 1.0 - c * c))
    tb, wb = target @ b, w @ b
    return b - (1 - c) * (target * tb + w * wb) + s * (target * wb - w * tb)


@dataclass
class ClusterSummary:
    cluster: int
    members: list  # (participant, component r) pairs
    centroid: np.ndarray  # v-bar (p,)
    b_aligned: dict  # participant -> sign-aligned, Procrustes-rotated b (p,)
    angles: np.ndarray  # angle of each member to the centroid (rad)
    b_bar: np.ndarray = field(init=False)
    t: np.ndarray = field(init=False)
    p: np.ndarray = field(init=False)

    def __post_init__(self):
        Bs = np.stack(list(self.b_aligned.values()))
        self.b_bar = Bs.mean(axis=0)
        self.t, self.p = stats.ttest_1samp(Bs, 0.0, axis=0)

    @property
    def beta(self) -> float:
        """Network persistence v-bar' b-bar."""
        return float(self.centroid @ self.b_bar)

    @property
    def off_network(self) -> np.ndarray:
        """(I - v-bar v-bar') b-bar: drive from outside the network."""
        return self.b_bar - self.beta * self.centroid


def summarise_clusters(V: dict, B: dict, labels: np.ndarray, owners: list) -> list[ClusterSummary]:
    """V, B: participant -> (p x k) loadings / factor coefficients; labels and owners give the cluster and the
    (participant, r) of each pooled vector, in the same order."""
    out = []
    for c in np.unique(labels):
        mem = [owners[i] for i in np.flatnonzero(labels == c)]
        vs = np.stack([V[s][:, r] for s, r in mem])
        u = np.linalg.eigh(vs.T @ vs)[1][:, -1]  # sign-invariant preliminary direction
        sgn = np.sign(vs @ u)
        sgn[sgn == 0] = 1
        vbar = (sgn[:, None] * vs).sum(axis=0)
        vbar /= np.linalg.norm(vbar)
        if vbar[np.argmax(np.abs(vbar))] < 0:
            vbar = -vbar
        b_aligned, angles = {}, []
        for (s, r), v in zip(mem, vs):
            sl = 1.0 if v @ vbar >= 0 else -1.0
            vb, bb = sl * v, sl * B[s][:, r]
            angles.append(np.arccos(np.clip(vb @ vbar, -1, 1)))
            b_aligned[s] = procrustes_to(vb, vbar, bb)
        out.append(ClusterSummary(int(c), mem, vbar, b_aligned, np.array(angles)))
    return out


def variance_shares(X: np.ndarray, Y: np.ndarray, V: np.ndarray, B: np.ndarray, r: int) -> dict:
    """In-sample variance of the r-th network score z = Y v_r predicted by its driving vector b = b_r, split
    into the within-network part (beta v_r) and the off-network part o = b - beta v_r, beta = v_r' b:
      R2_within = beta^2 v'Sv / var(z),  R2_off = o'So / var(z),  R2_cross = 2 beta v'So / var(z),
    with S = X'X / N; noise_off = sigma^2 (p - 1) / N / var(z) is the off-network share expected from
    estimation noise alone."""
    N, p = X.shape
    v, b = V[:, r], B[:, r]
    z = Y @ v
    S = X.T @ X / N
    var_z = z @ z / N
    beta = v @ b
    o = b - beta * v
    sigma2 = np.sum((z - X @ b) ** 2) / (N - p)
    return dict(R2_within=beta ** 2 * (v @ S @ v) / var_z, R2_off=(o @ S @ o) / var_z,
                R2_cross=2 * beta * (v @ S @ o) / var_z, noise_off=sigma2 * (p - 1) / N / var_z)


def match_to_reference(V: np.ndarray, V_ref: np.ndarray, keep: int):
    """Match clusters (rows of V, centroids) to reference centroids by Hungarian assignment on |cos|;
    returns [(cluster label, sign)] for reference clusters 1..keep."""
    C = V @ V_ref.T
    r, c = linear_sum_assignment(-np.abs(C))
    perm = np.empty(V_ref.shape[0], int)
    perm[c] = r
    return [(int(perm[j]) + 1, 1.0 if C[perm[j], j] >= 0 else -1.0) for j in range(keep)]
