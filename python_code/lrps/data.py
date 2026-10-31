"""Loading the processed AAL ROI time series and the DMCC task-block timing."""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .config import DATA_DIR, EVENTS_DIR

DATASETS = ("Axcpt", "Cuedts", "resting")
N_TUNE = 10  # the first 10 participants (sorted by file name) are used for tuning


def subject_files(dataset: str, data_dir: Path = DATA_DIR) -> list[Path]:
    files = sorted((Path(data_dir) / dataset).glob("*.csv"))
    if not files:
        raise FileNotFoundError(f"no CSV files for {dataset!r} under {Path(data_dir) / dataset}")
    return files


def subject_id(path: Path | str) -> str:
    m = re.match(r"(sub-[^_]+)", Path(path).name)
    return m.group(1) if m else Path(path).stem


def load_subject(path: Path | str) -> tuple[np.ndarray, list[str]]:
    """One participant's (n, p) ROI series, each ROI z-scored (missing values replaced by the ROI mean)."""
    df = pd.read_csv(path)
    y = df.to_numpy(dtype=float)
    if np.isnan(y).any():
        y = np.where(np.isnan(y), np.nanmean(y, axis=0), y)
    sd = y.std(axis=0, ddof=1)
    sd[sd == 0] = 1.0
    return (y - y.mean(axis=0)) / sd, list(df.columns)


def load_dataset(dataset: str, data_dir: Path = DATA_DIR, n_tune: int = N_TUNE) -> dict:
    """{'tune': {sid: y}, 'eval': {sid: y}, 'all': {sid: y}, 'roi_names': [...]}; tuning = first n_tune files."""
    out = {"tune": {}, "eval": {}, "roi_names": None}
    for i, f in enumerate(subject_files(dataset, data_dir)):
        y, names = load_subject(f)
        out["roi_names"] = out["roi_names"] or names
        out["tune" if i < n_tune else "eval"][subject_id(f)] = y
    out["all"] = {**out["tune"], **out["eval"]}
    return out


# ---- task blocks (DMCC runs: three task blocks separated by rest) ----------------------------------

def load_blocks(events_tsv: Path, bold_json: Path, n_tr: int, gap: float = 20.0) -> list[tuple[int, int]]:
    """Task-block TR spans [start, stop): trials are grouped into blocks wherever consecutive onsets are more
    than `gap` seconds apart; a block spans [first onset, last onset + duration], converted to TRs."""
    ev = pd.read_csv(events_tsv, sep="\t").sort_values("onset").reset_index(drop=True)
    tr = float(json.loads(Path(bold_json).read_text())["RepetitionTime"])
    breaks = np.flatnonzero(np.diff(ev.onset.to_numpy()) > gap)
    starts, ends = np.r_[0, breaks + 1], np.r_[breaks, len(ev) - 1]
    blocks = []
    for a, b in zip(starts, ends):
        t0 = int(np.floor(ev.onset[a] / tr))
        t1 = int(np.ceil((ev.onset[b] + ev.duration[b]) / tr))
        blocks.append((max(t0, 0), min(t1, n_tr)))
    return blocks


def subject_blocks(task: str, sid: str, n_tr: int, events_dir: Path = EVENTS_DIR) -> list[tuple[int, int]]:
    tsv = next((Path(events_dir) / task).glob(f"{sid}_*_events.tsv"))
    return load_blocks(tsv, Path(str(tsv).replace("_events.tsv", "_bold.json")), n_tr)


def block_pairs(y: np.ndarray, block: tuple[int, int], lag: int = 1):
    """Lag pairs (y_{t-lag}, y_t) with both time points inside the block."""
    a, b = block
    seg = y[a:b]
    return seg[:-lag], seg[lag:]


def lag_pairs(dataset: str, sid: str, y: np.ndarray):
    """All of a participant's data as (X, Y) lag pairs: the whole run for the resting state; for the task runs
    the lag pairs inside the three task blocks, stacked, never bridging a break between blocks."""
    if dataset == "resting":
        return y[:-1], y[1:]
    parts = [block_pairs(y, b) for b in subject_blocks(dataset, sid, len(y))]
    return np.vstack([x for x, _ in parts]), np.vstack([t for _, t in parts])
