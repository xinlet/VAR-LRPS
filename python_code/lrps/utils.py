"""Parallel map and tuning-parameter consensus."""

from __future__ import annotations

import os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor


def parallel_map(fn, tasks, n_jobs: int = 1) -> list:
    """[fn(t) for t in tasks], optionally across processes (one BLAS thread per worker)."""
    tasks = list(tasks)
    if n_jobs == 1 or len(tasks) <= 1:
        return [fn(t) for t in tasks]
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(v, "1")  # inherited by the worker processes
    with ProcessPoolExecutor(max_workers=min(n_jobs, len(tasks))) as ex:
        return list(ex.map(fn, tasks, chunksize=max(1, len(tasks) // (8 * n_jobs))))


def consensus(choices, prefer: str = "min"):
    """Most frequent value; ties broken towards the smallest ('min') or largest ('max') candidate."""
    counts = Counter(choices)
    top = max(counts.values())
    tied = [v for v, c in counts.items() if c == top]
    return min(tied) if prefer == "min" else max(tied)
