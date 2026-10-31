"""Design of the synthetic experiments (Section 3)."""

from __future__ import annotations

import numpy as np

P, N, WINDOW = 30, 100, 50
RHOS = (0.1, 0.5, 0.9)
EXPERIMENTS = {"1A": "flat, off-diagonal 0.2", "1B": "flat, off-diagonal 0.01", "2": "decaying (Toeplitz)"}
K_GRID = np.arange(1, P + 1)
LAM_GRID = 10.0 ** np.linspace(-4, 0, 20)
N_TUNE, N_EVAL, N_EST = 50, 100, 100


def true_A(exp: str, rho: float, p: int = P) -> np.ndarray:
    """1A: A_ii = rho, A_ij = 0.2; 1B: A_ij = 0.01; 2: A_ij = rho^|i-j|; rescaled to A / (gamma_max + 0.1) when
    the largest eigenvalue modulus gamma_max >= 1."""
    if exp in ("1A", "1B"):
        A = np.full((p, p), 0.2 if exp == "1A" else 0.01)
        np.fill_diagonal(A, rho)
    else:
        idx = np.arange(p)
        A = rho ** np.abs(idx[:, None] - idx[None, :])
    gmax = np.max(np.abs(np.linalg.eigvals(A)))
    return A / (gmax + 0.1) if gmax >= 1 else A
