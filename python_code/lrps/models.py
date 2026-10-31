"""VAR(1) estimators: OLS, LRPS (low-rank pre-smoothing) and the sparse VAR (lasso, via scikit-learn).

Row convention: with lag pairs X = (y_1, ..., y_{n-1})' and Y = (y_2, ..., y_n)', the model is Y = X B + E,
so B = A' in the notation y_t = A y_{t-1} + e_t of the paper, and B[j, i] is the effect of region j at t-1
on region i at t. No intercept (the series are centred).
"""

from __future__ import annotations

import warnings

import numpy as np


def fit_ols(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """OLS VAR coefficients B = (X'X)^{-1} X'Y (least squares; pseudo-inverse if X'X is singular)."""
    return np.linalg.lstsq(X, Y, rcond=None)[0]


def lrps_path(X: np.ndarray, Y: np.ndarray):
    """B_ols and the right singular vectors V of the response Y; the rank-k LRPS estimate is B_ols V_k V_k'
    (the OLS regression of the rank-k SVD approximation Y V_k V_k' on X)."""
    return fit_ols(X, Y), np.linalg.svd(Y, full_matrices=False)[2].T


def fit_lrps(X: np.ndarray, Y: np.ndarray, k: int) -> np.ndarray:
    B, V = lrps_path(X, Y)
    return B @ V[:, :k] @ V[:, :k].T


def lrps_factors(X: np.ndarray, Y: np.ndarray, k: int):
    """V_k (p x k) and the factor-level coefficients B_fac = (X'X)^{-1} X' Y V_k (p x k), so that the LRPS
    estimate is B_fac V_k' = sum_r b_r v_r'."""
    V = np.linalg.svd(Y, full_matrices=False)[2][:k].T
    return V, fit_ols(X, Y @ V)


def lasso_path_columns(X: np.ndarray, Y: np.ndarray, lams, tol: float = 1e-8, max_iter: int = 100_000):
    """Sparse VAR over a grid of penalties with scikit-learn: each column of Y is regressed separately on X with
    the lasso objective (1 / 2n) ||y - X b||^2 + lam ||b||_1 (no intercept, no standardisation), the same
    scaling as glmnet. Returns an array (len(lams), p_x, p_y), in the order of `lams`.

    Each response is fitted separately (sklearn's lasso_path with a 2-D target would fit a multi-task lasso)."""
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import lasso_path

    lams = np.atleast_1d(np.asarray(lams, dtype=float))
    order = np.argsort(-lams)  # sklearn computes the path from the largest penalty down
    G = X.T @ X
    out = np.zeros((lams.size, X.shape[1], Y.shape[1]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        for m in range(Y.shape[1]):
            y = Y[:, m]
            _, coefs, _ = lasso_path(X, y, alphas=lams[order], precompute=G, Xy=X.T @ y, tol=tol,
                                     max_iter=max_iter)
            out[order, :, m] = coefs.T
    return out


def fit_svar(X: np.ndarray, Y: np.ndarray, lam: float) -> np.ndarray:
    return lasso_path_columns(X, Y, [lam])[0]


def lagged(y: np.ndarray, lag: int = 1):
    return y[:-lag], y[lag:]


# ---- simulation -----------------------------------------------------------------------------------------

def simulate_var(A: np.ndarray, n: int, rng: np.random.Generator, burn_in: int = 0) -> np.ndarray:
    """y_1 = 0, y_t = A y_{t-1} + e_t, e_t ~ N(0, I) (as tsDyn::VAR.sim with zero starting values);
    returns an (n, p) array (rows = time)."""
    p = A.shape[0]
    total = n + burn_in
    y = np.zeros((total, p))
    e = rng.standard_normal((total, p))
    for t in range(1, total):
        y[t] = A @ y[t - 1] + e[t]
    return y[burn_in:]
