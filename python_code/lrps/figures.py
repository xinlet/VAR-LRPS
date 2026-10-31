"""Manuscript figures (PLOS ONE: 5.2 in text column, S3 Fig 7.5 in; Arial, text 8-12 pt), drawn from the saved
results. Each analysis script calls make_figures(<step>) when it finishes; every figure is written as PDF, SVG and
300 dpi LZW TIFF to <results>/figures.

  step  figures
  01    fig1_eigenvalues, fig2_estimation_mse, fig3_rank_1B_2, figS1_rank_1A_1B, figS2_cv_curves
  02    fig4_cv_scheme, fig5_cv_curves
  03    fig6_resting_surfaces, fig7_resting_by_size
  04    fig8_edges, figS3_edges_lrps_k
  05    fig9_clusters, fig10_driving
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .clusters import match_to_reference
from .config import RESULTS_DIR
from .data import subject_blocks, subject_files, subject_id
from .glassbrain import DIVERGING, aal_coords, plot_edges, plot_nodes
from .simulation import K_GRID, P, RHOS, true_A
from .stats import mean_ci

W = 5.2  # PLOS text-column width (in)
FS, FS_HEAD, FS_PANEL = 8, 9, 11
INK, INK2, GRID, AXIS = "#0b0b0b", "#52514e", "#e4e3df", "#8a8984"
C = {"LRPS": "#2a78d6", "LRPS_light": "#9fc2ec", "VAR": "#eb6834", "sVAR": "#1baf7a"}
C_RESTING = {"VAR": "#2a78d6", "LRPS": "#eb6834", "sVAR": "#1baf7a"}
SEQ = {1: "#7fb0e8", 2: "#1f5fae"}  # one / two training blocks
TRAIN_C, TEST_C, IDLE_C = "#2a78d6", "#eb6834", "#e4e3df"
DATASET_COLOR = {"Axcpt": "#2a78d6", "Cuedts": "#eb6834", "resting": "#1baf7a"}
DATASET_TINT = {"Axcpt": "#b3cff1", "Cuedts": "#f8c8b2", "resting": "#a9e3cc"}
SEP = "#c9c7c2"
EXP_NAME = {"1A": "Experiment 1A", "1B": "Experiment 1B", "2": "Experiment 2"}
LABEL = {"VAR": "VAR (OLS)", "LRPS": "LRPS", "sVAR": "sVAR"}
DATASETS = ("Axcpt", "Cuedts", "resting")
S = W / 7.5  # graphics scale of the glass-brain figures relative to the full page width

R = RESULTS_DIR
SIM, BCV, REST, CONN, CLU = (R / "simulations", R / "block_cv", R / "resting", R / "connectivity", R / "clusters")


# ---- helpers --------------------------------------------------------------------------------------------------

def setup():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42, "font.family": "Arial", "font.size": FS,
                         "axes.labelsize": FS, "axes.titlesize": FS, "xtick.labelsize": FS, "ytick.labelsize": FS,
                         "legend.fontsize": FS, "mathtext.fontset": "custom", "mathtext.rm": "Arial",
                         "mathtext.it": "Arial:italic", "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                         "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
                         "xtick.minor.size": 0, "ytick.minor.size": 0})
    return plt


class Lay:
    """Axes and text placed in inches from the top-left corner of a width x H figure."""

    def __init__(self, plt, H, width=W):
        self.fig, self.H, self.W = plt.figure(figsize=(width, H)), H, width

    def rect(self, x, y, w, h):
        return [x / self.W, 1 - (y + h) / self.H, w / self.W, h / self.H]

    def ax(self, x, y, w, h, **kw):
        return self.fig.add_axes(self.rect(x, y, w, h), **kw)

    def text(self, x, y, s, **kw):
        self.fig.text(x / self.W, 1 - y / self.H, s, **kw)

    def panel(self, x, y, letter, heading=""):
        self.text(x, y, f"{letter})", fontsize=FS_PANEL, fontweight="bold", color=INK, va="baseline")
        if heading:
            self.text(x + 0.27, y, heading, fontsize=FS_HEAD, color=INK, va="baseline")


def style(ax, grid="y", scale=1.0, color=INK):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(color)
    ax.tick_params(colors=color, labelsize=FS)
    if scale != 1.0:
        ax.spines[["left", "bottom"]].set_linewidth(0.8 * scale)
        ax.tick_params(width=0.8 * scale, length=3 * scale)
    if grid:
        ax.grid(axis=grid, color=GRID, lw=0.5)
    ax.set_axisbelow(True)


def plain_log_ticks(axis, ticks):
    """Log axis with plain decimal tick labels (no superscripts)."""
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

    axis.set_major_locator(FixedLocator(ticks))
    axis.set_minor_locator(NullLocator())
    axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))


def save(fig, out, name):
    import matplotlib.pyplot as plt

    for ext in ("pdf", "svg"):
        fig.savefig(out / f"{name}.{ext}")
    fig.savefig(out / f"{name}.tif", dpi=300, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("wrote", name)


# ---- Section 3 -------------------------------------------------------------------------------------------------

def fig1(plt, out):
    from matplotlib.lines import Line2D

    sty = {"1A": ("#e34948", "s"), "1B": ("#2a78d6", "o"), "2": ("#1baf7a", "^")}
    lay = Lay(plt, 2.25)
    for i, rho in enumerate(RHOS):
        ax = lay.ax(0.5 + 1.6 * i, 0.25, 1.38, 1.35)
        for exp, (col, mk) in sty.items():
            ev = np.sort(np.abs(np.linalg.eigvals(true_A(exp, rho))))[::-1]
            ax.plot(np.arange(1, P + 1), ev, color=col, lw=0.9, marker=mk, ms=2.6, mec="white", mew=0.3)
        ax.set_title(f"ρ = {rho}", loc="left", fontsize=FS, color=INK, pad=3)
        ax.set_xticks([1, 10, 20, 30])
        ax.set_ylim(-0.03, 1.03)
        ax.set_xlabel("eigenvalue index", color=INK)
        style(ax, grid="both")
        if i:
            ax.tick_params(labelleft=False)
        else:
            ax.set_ylabel("|eigenvalue| of A", color=INK)
    handles = [Line2D([], [], color=col, marker=mk, ms=3.5, lw=0.9, mec="white", mew=0.3, label=EXP_NAME[e])
               for e, (col, mk) in sty.items()]
    lay.fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0),
                   handlelength=1.8, columnspacing=1.5)
    save(lay.fig, out, "fig1_eigenvalues")


def fig2(plt, out):
    from matplotlib.lines import Line2D

    est = pd.read_csv(SIM / "estimation_mse.csv", dtype={"exp": str})
    lam = np.asarray(json.loads((SIM / "settings.json").read_text())["lambda_grid"])
    k = K_GRID
    lay = Lay(plt, 6.7)
    names = {"1A": "Experiment 1A: flat, off-diagonal 0.2", "1B": "Experiment 1B: flat, off-diagonal 0.01",
             "2": "Experiment 2: decaying (Toeplitz)"}
    row_h, top0 = 1.92, 0.67
    for j, rho in enumerate(RHOS):
        lay.text(0.55 + 1.58 * j + 0.7, 0.2, f"ρ = {rho}", fontsize=FS_HEAD, color=INK, ha="center", va="baseline")
    for i, exp in enumerate(("1A", "1B", "2")):
        g_rows = est[est.exp == exp]
        curves = {}
        for rho in RHOS:
            g = g_rows[g_rows.rho == rho]
            curves[rho] = (np.log(np.stack(g.lrps.apply(json.loads).to_numpy())).mean(axis=0),
                           np.log(np.stack(g.lasso.apply(json.loads).to_numpy())).mean(axis=0),
                           np.log(g["var"].to_numpy()).mean())
        lo = min(min(c[0].min(), c[1].min(), c[2]) for c in curves.values())
        hi = max(max(c[0].max(), c[1].max(), c[2]) for c in curves.values())
        pad = 0.06 * (hi - lo)
        y0 = top0 + row_h * i
        lay.panel(0.0, y0, "ABC"[i], names[exp])
        for j, rho in enumerate(RHOS):
            lr, la, vr = curves[rho]
            ax = lay.ax(0.55 + 1.58 * j, y0 + 0.5, 1.4, 1.0)
            ax.axhline(vr, color=C["VAR"], lw=1.0, ls="--", zorder=2)
            ax.plot(k, lr, color=C["LRPS"], lw=0.9, marker="o", ms=2.2, mec="white", mew=0.3, zorder=3)
            ax.set_xlim(0, k[-1] + 1)
            ax.set_xticks([1, 10, 20, 30])
            ax.set_ylim(lo - pad, hi + pad)
            ax.tick_params(axis="both", colors=INK, labelsize=FS, pad=1.5)
            ax.grid(color=GRID, lw=0.5)
            ax.set_axisbelow(True)
            if j:
                ax.tick_params(axis="y", labelleft=False)
            else:
                ax.set_ylabel("mean log MSE", color=INK)
            top = ax.twiny()  # sVAR against log10(lambda) on the top axis
            top.plot(np.log10(lam), la, color=C["sVAR"], lw=0.9, marker="^", ms=2.5, mec="white", mew=0.3, zorder=3)
            top.set_xlim(np.log10(lam[0]) - 0.2, np.log10(lam[-1]) + 0.2)
            top.set_xticks([-4, -3, -2, -1, 0])
            top.tick_params(axis="x", colors=INK, labelsize=FS, pad=1)
            for a in (ax, top):
                a.spines["right"].set_visible(False)
                a.spines["left"].set_color(INK)
            ax.spines["top"].set_visible(False)
            top.spines["bottom"].set_visible(False)
            top.spines["left"].set_visible(False)
            top.spines["top"].set_color(INK)
            ax.spines["bottom"].set_color(INK)
            if i == 2:
                ax.set_xlabel("rank k (LRPS)", color=INK, labelpad=1)
            if i == 0:
                top.set_xlabel("log10(λ) (sVAR)", color=INK, labelpad=2)
    handles = [Line2D([], [], color=C["LRPS"], lw=0.9, marker="o", ms=3, mec="white", label="LRPS (bottom axis: k)"),
               Line2D([], [], color=C["sVAR"], lw=0.9, marker="^", ms=3.5, mec="white",
                      label="sVAR (top axis: log10(λ))"),
               Line2D([], [], color=C["VAR"], lw=1.0, ls="--", label="VAR (OLS)")]
    lay.fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0),
                   handlelength=1.8, columnspacing=1.2)
    save(lay.fig, out, "fig2_estimation_mse")


def fig_rank(plt, out, exps, name, kmax=15):
    sel = pd.read_csv(SIM / "rank_selection.csv", dtype={"exp": str})
    ev = pd.read_csv(SIM / "eval_mspe.csv", dtype={"exp": str})
    lay = Lay(plt, 5.35)
    for b, exp in enumerate(exps):
        x0 = 2.62 * b
        lay.panel(x0, 0.22, "AB"[b], EXP_NAME[exp])
        ymax = max(np.bincount(sel[(sel.exp == exp) & (sel.rho == r)].k_selected).max() for r in RHOS)
        for i, rho in enumerate(RHOS):
            y = 0.5 + 1.55 * i
            s = sel[(sel.exp == exp) & (sel.rho == rho)]
            kh = int(s.k_hat.iloc[0])
            cnt = np.bincount(s.k_selected, minlength=P + 1)[1:]
            if cnt[kmax:].sum():
                print(f"note: {exp}, rho = {rho}: some tuning series select k > {kmax} (not shown)")
            ax = lay.ax(x0 + 0.42, y, 1.3, 1.05)
            ax.bar(np.arange(1, kmax + 1), cnt[:kmax], width=0.8,
                   color=[C["LRPS"] if kk == kh else C["LRPS_light"] for kk in range(1, kmax + 1)])
            ax.set_xlim(0.3, kmax + 0.7)
            ax.set_ylim(0, ymax * 1.08)
            ax.set_xticks([1, 5, 10, 15])
            ax.set_title(rf"ρ = {rho}:  $\hat{{k}}$ = {kh}", loc="left", fontsize=FS, color=INK, pad=3)
            ax.set_ylabel("tuning series", color=INK)
            if i == 2:
                ax.set_xlabel("selected rank k", color=INK)
            style(ax)
            e = ev[(ev.exp == exp) & (ev.rho == rho)]
            ax = lay.ax(x0 + 2.08, y, 0.48, 1.05)
            bp = ax.boxplot([e.mspe_var, e.mspe_lrps], widths=0.6, patch_artist=True,
                            medianprops=dict(color=INK, lw=0.8), flierprops=dict(ms=1.8, mec=AXIS, mew=0.4),
                            whiskerprops=dict(color=AXIS, lw=0.6), capprops=dict(color=AXIS, lw=0.6),
                            boxprops=dict(lw=0.6))
            for patch, m in zip(bp["boxes"], ("VAR", "LRPS")):
                patch.set_facecolor(C[m])
                patch.set_edgecolor(AXIS)
            ax.set_xticks([1, 2], ["VAR", "LRPS"], rotation=90)
            ax.set_ylim(0.9 * min(e.mspe_var.min(), e.mspe_lrps.min()), None)
            ax.set_ylabel("MSPE", color=INK, labelpad=1)
            style(ax)
    save(lay.fig, out, name)


def figS2(plt, out):
    from matplotlib.lines import Line2D

    curves = pd.read_csv(SIM / "cv_curves_tuning.csv", dtype={"exp": str})
    sel = pd.read_csv(SIM / "selection_comparison.csv", dtype={"exp": str})
    lay = Lay(plt, 5.55)
    for b, exp in enumerate(("1B", "2")):
        x0 = 2.62 * b
        lay.panel(x0, 0.22, "AB"[b], EXP_NAME[exp])
        for i, rho in enumerate(RHOS):
            ax = lay.ax(x0 + 0.45, 0.5 + 1.55 * i, 2.08, 1.05)
            c = curves[(curves.exp == exp) & (curves.rho == rho)]
            mat = c.pivot(index="rep", columns="k", values="cvMSPE").to_numpy()
            for row in mat:
                ax.plot(K_GRID, row, color="#9fc2ec", lw=0.4, alpha=0.6, zorder=1)
            ax.plot(K_GRID, mat.mean(axis=0), color=C["LRPS"], lw=1.4, zorder=3)
            s = sel[(sel.exp == exp) & (sel.rho == rho)].iloc[0]
            ax.axvline(s.k_consensus, color=INK, lw=0.8, ls="--", zorder=2)
            ax.axvline(s.k_average, color=C["VAR"], lw=0.9, zorder=2)
            ax.set_xlim(0.5, P + 0.5)
            ax.set_xticks([1, 10, 20, 30])
            kk = (rf"$\hat k$ = {s.k_consensus} (both)" if s.k_consensus == s.k_average
                  else rf"$\hat k$ = {s.k_consensus} (consensus), {s.k_average} (average)")
            ax.set_title(f"ρ = {rho}: {kk}", loc="left", fontsize=FS, color=INK, pad=3)
            ax.set_ylabel("cvMSPE", color=INK)
            if i == 2:
                ax.set_xlabel("rank k", color=INK)
            style(ax, grid="both")
    handles = [Line2D([], [], color="#9fc2ec", lw=0.8, label="tuning series (50)"),
               Line2D([], [], color=C["LRPS"], lw=1.4, label="average cvMSPE"),
               Line2D([], [], color=INK, lw=0.8, ls="--", label=r"consensus $\hat k$"),
               Line2D([], [], color=C["VAR"], lw=0.9, label=r"average $\hat k$")]
    lay.fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0),
                   handlelength=1.6, columnspacing=1.0)
    save(lay.fig, out, "figS2_cv_curves")


# ---- Section 4.2: block cross-validation -------------------------------------------------------------------------

def block_timing(n_tr=610):
    """Median task-block spans (TRs) over the Axcpt participants."""
    spans = np.array([[x for b in subject_blocks("Axcpt", subject_id(f), n_tr) for x in b]
                      for f in subject_files("Axcpt")])
    return np.median(spans, axis=0).astype(int).reshape(3, 2).tolist()


def fig4(plt, out):
    from matplotlib.patches import FancyArrowPatch, Rectangle

    blocks, n_tr = block_timing(), 610
    fig = plt.figure(figsize=(W, 2.55))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(-125, n_tr + 6)
    ax.set_ylim(0.38, 4.62)
    ax.axis("off")

    def run_row(y, roles, label):
        ax.text(-121, y + 0.2, label, fontsize=FS, va="center", ha="left", color=INK)
        ax.add_patch(Rectangle((0, y), n_tr, 0.4, facecolor="white", edgecolor=AXIS, lw=0.5))
        for b, (a, z) in enumerate(blocks):
            role = roles[b]
            ax.add_patch(Rectangle((a, y), z - a, 0.4, facecolor={"train": TRAIN_C, "test": TEST_C, "": IDLE_C}[role],
                                   edgecolor="white", lw=0.8))
            ax.text((a + z) / 2, y + 0.2, f"block {b + 1}: {role or 'not used'}", fontsize=FS, ha="center",
                    va="center", color="white" if role else INK2)

    y0 = 3.75
    ax.add_patch(Rectangle((0, y0), n_tr, 0.4, facecolor="white", edgecolor=AXIS, lw=0.5))
    for b, (a, z) in enumerate(blocks):
        ax.add_patch(Rectangle((a, y0), z - a, 0.4, facecolor="#c9c7c1", edgecolor="white", lw=0.8))
        ax.text((a + z) / 2, y0 + 0.2, f"task block {b + 1}", fontsize=FS, ha="center", va="center", color=INK)
    rests = [(0, blocks[0][0]), (blocks[0][1], blocks[1][0]), (blocks[1][1], blocks[2][0]), (blocks[2][1], n_tr)]
    for a, z in rests:
        ax.text((a + z) / 2, y0 - 0.2, "rest", fontsize=FS, ha="center", va="center", color=INK2)
    ax.text(-121, y0 + 0.2, f"one run\n({n_tr} TRs)", fontsize=FS, va="center", ha="left", color=INK)
    for t in (0, *[x for b in blocks for x in b], n_tr):
        ax.plot([t, t], [y0 + 0.4, y0 + 0.47], color=AXIS, lw=0.5)
    ax.text(n_tr / 2, y0 + 0.55, "time (TR = 1.2 s)", fontsize=FS, ha="center", va="bottom", color=INK2)
    ax.text(-121, 3.1, "(a) train on one block, test on the next (2 evaluations per subject)", fontsize=FS_HEAD,
            color=INK, fontweight="bold")
    run_row(2.42, ["train", "test", ""], "split 1")
    run_row(1.82, ["", "train", "test"], "split 2")
    ax.text(-121, 1.22, "(b) train on two blocks, test on the next (1 evaluation per subject)", fontsize=FS_HEAD,
            color=INK, fontweight="bold")
    run_row(0.55, ["train", "train", "test"], "split 1")
    for yy, a, b in ((2.42, blocks[0][1], blocks[1][0]), (1.82, blocks[1][1], blocks[2][0]),
                     (0.55, blocks[1][1], blocks[2][0])):
        ax.add_patch(FancyArrowPatch((a - 4, yy + 0.47), (b + 14, yy + 0.47), arrowstyle="-|>", mutation_scale=6,
                                     lw=0.7, color=INK, connectionstyle="arc3,rad=-0.35"))
    save(fig, out, "fig4_cv_scheme")


def fig5(plt, out):
    lay = Lay(plt, 3.72)
    for i, task in enumerate(("Axcpt", "Cuedts")):
        for j, kind in enumerate(("LRPS", "sVAR")):
            ax = lay.ax(0.55 + 2.4 * j, 0.25 + 1.62 * i, 2.05, 1.1)
            for n in (1, 2):
                d = BCV / f"train{n}" / task
                t = pd.read_csv(d / "tuning_split_errors.csv")
                s = json.loads((d / "summary.json").read_text())
                c = t[t.method == kind].groupby(["subject", "param"]).mse.mean().unstack().mean()
                ax.plot(c.index.to_numpy(float), np.sqrt(c.values), color=SEQ[n], lw=1.3,
                        label=f"{n} training block{'s' if n > 1 else ''}")
                ax.axvline(s["k_hat"] if kind == "LRPS" else s["lambda_hat"], color=SEQ[n], lw=0.8, ls="--")
            if kind == "sVAR":
                ax.set_xscale("log")
                plain_log_ticks(ax.xaxis, [0.001, 0.01, 0.1, 1])
                ax.set_xlabel("sVAR penalty λ", color=INK)
            else:
                ax.set_xticks([1, 5, 10, 15, 20])
                ax.set_xlabel("rank k", color=INK)
            ax.set_ylabel("CV RMSE", color=INK)
            ax.set_title(f"{task}: {kind}", loc="left", fontsize=FS_HEAD, color=INK, pad=3)
            style(ax)
    h, lab = lay.fig.axes[0].get_legend_handles_labels()
    lay.fig.legend(h, lab, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.0), columnspacing=2)
    save(lay.fig, out, "fig5_cv_curves")


# ---- Section 4.3: expanding window, resting state ------------------------------------------------------------------

def _resting():
    s = json.loads((REST / "summary.json").read_text())
    return s, pd.read_csv(REST / "tuning_errors_by_origin.csv.gz"), pd.read_csv(REST / "eval_errors_by_origin.csv.gz")


def _nice_ticks(lo, hi, n=5):
    cand = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
    t = [c for c in cand if lo <= c <= hi]
    while len(t) > n:
        t = t[::2]
    return t


def fig6(plt, out):
    from matplotlib.colors import LogNorm
    from matplotlib.lines import Line2D
    from matplotlib.ticker import NullLocator

    s, tune, _ = _resting()
    lay = Lay(plt, 2.75)
    specs = [("LRPS", s["k_grid"], s["k_hat"], "rank k", False),
             ("sVAR", s["lambda_grid"], s["lambda_hat"], "sVAR penalty λ", True)]
    for i, (m, grid, sel, xlabel, logx) in enumerate(specs):
        Sf = tune[tune.method == m].groupby(["train_size", "param"]).rmse.mean().unstack("param")
        sizes = Sf.index.to_numpy()
        x = np.array(grid, float)
        if logx:
            lx = np.log10(x)
            step = lx[1] - lx[0]
            edges = 10 ** np.r_[lx - step / 2, lx[-1] + step / 2]
        else:
            edges = np.r_[x - 0.5, x[-1] + 0.5]
        yedges = np.r_[sizes - 0.5, sizes[-1] + 0.5]
        vals = Sf.values
        vmin, vmax = np.nanmin(vals), np.nanpercentile(vals, 95)
        ax = lay.ax(0.48 + 2.48 * i, 0.25, 1.62, 1.85)
        im = ax.pcolormesh(edges, yedges, vals, cmap="Blues", norm=LogNorm(vmin=vmin, vmax=vmax), shading="flat",
                           rasterized=True)
        ax.plot(x[np.nanargmin(vals, axis=1)], sizes, "x", ms=3.8, mew=1.0, color="#eb6834", ls="none")
        ax.axvline(sel, color=INK, lw=0.9, ls="--")
        if logx:
            ax.set_xscale("log")
            plain_log_ticks(ax.xaxis, [0.001, 0.01, 0.1, 1])
        else:
            ax.set_xticks([1, 5, 10, 15, 20])
        ax.set_xlabel(xlabel, color=INK)
        if i == 0:
            ax.set_ylabel("training window size t", color=INK)
        ax.set_title(f"{LABEL[m]} (selected {'k' if m == 'LRPS' else 'λ'} = {sel:.3g})", loc="left",
                     fontsize=FS, color=INK, pad=3)
        ax.tick_params(colors=INK, labelsize=FS)
        cb = lay.fig.colorbar(im, cax=lay.ax(0.48 + 2.48 * i + 1.68, 0.25, 0.07, 1.85))
        ticks = _nice_ticks(vmin, vmax)
        cb.set_ticks(ticks)
        cb.set_ticklabels([f"{t:g}" for t in ticks])
        cb.ax.yaxis.set_minor_locator(NullLocator())
        cb.ax.tick_params(labelsize=FS, colors=INK)
        if i == 1:
            cb.set_label("RMSE (mean over tuning subjects)", fontsize=FS, color=INK)
    handles = [Line2D([], [], color="#eb6834", marker="x", ls="none", ms=4.5, mew=1.0,
                      label="best value at each training size"),
               Line2D([], [], color=INK, lw=0.9, ls="--", label="selected (min average CV RMSE)")]
    lay.fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.0))
    save(lay.fig, out, "fig6_resting_surfaces")


def fig7(plt, out):
    s, tune, ev = _resting()
    k_hat, lam_hat = s["k_hat"], s["lambda_hat"]
    lay = Lay(plt, 2.55)
    axes = [lay.ax(0.5, 0.25, 2.2, 1.6), lay.ax(2.9, 0.25, 2.2, 1.6)]
    top = 0
    for ax, (df, panel) in zip(axes, ((tune, "tuning subjects"), (ev, "held-out subjects"))):
        for m, par in (("VAR", None), ("LRPS", k_hat), ("sVAR", lam_hat)):
            d = df[df.method == m] if par is None else df[(df.method == m) & np.isclose(df.param, par)]
            g = d.groupby("train_size").rmse
            med, lo, hi = g.median(), g.quantile(0.25), g.quantile(0.75)
            top = max(top, hi.max())
            ax.fill_between(med.index, lo, hi, color=C_RESTING[m], alpha=0.15, lw=0)
            ax.plot(med.index, med.values, color=C_RESTING[m], lw=1.2,
                    label=LABEL[m] + ("" if par is None else (f" (k = {par})" if m == "LRPS" else f" (λ = {par:.3g})")))
        ax.set_xlabel("training window size t", color=INK)
        ax.set_title(panel, loc="left", fontsize=FS, color=INK, pad=3)
        style(ax)
    for i, ax in enumerate(axes):
        ax.set_ylim(0, top * 1.05)
        if i:
            ax.tick_params(labelleft=False)
    axes[0].set_ylabel("one-step RMSE", color=INK)
    h, lab = axes[0].get_legend_handles_labels()
    lay.fig.legend(h, lab, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0))
    save(lay.fig, out, "fig7_resting_by_size")


# ---- Section 4.4: connectivity ---------------------------------------------------------------------------------------

def _edge_norm(res):
    import matplotlib.colors as mcolors

    vmax = max(np.nanmax(np.abs(z["t"][z["bh_edges"]])) for z in res.values() if z["bh_edges"].any())
    return (mcolors.LinearSegmentedColormap.from_list("div", DIVERGING),
            mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax))


def edge_grid(plt, cols, width, H, gs_kw, scale, fs_small, fs_head, name, out, legend_xy, cbar,
              bonf_word="Bonf."):
    """Rows = datasets, columns = estimators x (coronal, axial); BH-significant edges, Bonferroni ones thicker."""
    from matplotlib.cm import ScalarMappable
    from matplotlib.lines import Line2D

    coords, _ = aal_coords()
    res = {(ds, lab): np.load(CONN / ds / f"estimates_{f}.npz") for ds in DATASETS for lab, f in cols}
    cmap, norm = _edge_norm(res)
    fig = plt.figure(figsize=(width, H))
    gs = fig.add_gridspec(len(DATASETS), 2 * len(cols), **gs_kw)
    grid = {}
    for i, ds in enumerate(DATASETS):
        for j, (lab, _) in enumerate(cols):
            z = res[(ds, lab)]
            src, dst = np.nonzero(z["bh_edges"])
            bonf = z["bonf_edges"][src, dst]
            axes = [fig.add_subplot(gs[i, 2 * j]), fig.add_subplot(gs[i, 2 * j + 1])]
            grid[i, j] = axes
            plot_edges(axes, coords, src, dst, z["t"][src, dst], bonf, cmap, norm, scale=scale)
            axes[0].text(0.0, 1.03, f"{len(src)} BH, {int(bonf.sum())} {bonf_word}", transform=axes[0].transAxes,
                         fontsize=fs_small, color=INK2, ha="left", va="bottom")
            print(f"  {ds:8s} {lab:14s} BH {len(src):4d}  Bonferroni {int(bonf.sum()):3d}")
            if i == 0:
                axes[0].text(0.0, 1.24, lab, transform=axes[0].transAxes, fontsize=fs_head, fontweight="bold",
                             color=INK, ha="left", va="bottom")
            if j == 0:
                p = axes[0].get_position()
                fig.text(p.x0 - 0.01, p.y0 + p.height / 2, ds, rotation=90, fontsize=fs_head, fontweight="bold",
                         color=DATASET_COLOR[ds], ha="right", va="center")
    y_top = grid[0, 0][0].get_position().y1 + 0.32 / H
    y_bot = grid[len(DATASETS) - 1, 0][0].get_position().y0
    for j in range(1, len(cols)):
        x = (grid[0, j - 1][1].get_position().x1 + grid[0, j][0].get_position().x0) / 2
        fig.add_artist(Line2D([x, x], [y_bot, y_top], transform=fig.transFigure, color=SEP, lw=0.8 * scale))
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=fig.add_axes(cbar), orientation="horizontal")
    cb.set_label("t-statistic of the mean coefficient across participants", fontsize=fs_small, color=INK2)
    cb.ax.tick_params(labelsize=fs_small, colors=INK2)
    handles = [Line2D([], [], color=INK2, lw=0.7 * scale, label="significant at FDR 5% (BH)"),
               Line2D([], [], color=INK2, lw=1.5 * scale, label="also Bonferroni-significant (95%)")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=legend_xy, fontsize=fs_small, frameon=False,
               handlelength=2.0)
    save(fig, out, name)


def fig8(plt, out):
    row = 0.3 + (1.45 - 0.3) * S
    H = 0.55 + row * len(DATASETS) + 0.95
    edge_grid(plt, [("VAR (OLS)", "VAR"), ("LRPS (k = 10)", "LRPS"), ("sVAR", "sVAR")], W, H,
              dict(left=0.07, right=0.995, top=1 - 0.55 / H, bottom=0.95 / H, wspace=0.03, hspace=0.42), S, 8,
              max(8, 10 * S), "fig8_edges", out, (0.55, 0.38 / H), [0.07, 0.5 / H, 0.38, 0.1 / H],
              bonf_word="Bonferroni")


def figS3(plt, out):
    H = 4.75
    cols = [("VAR (OLS)", "VAR")] + [(f"LRPS (k = {k})", f) for k, f in ((5, "LRPSk5"), (10, "LRPS"), (20, "LRPSk20"))]
    edge_grid(plt, cols, 7.5, H, dict(left=0.045, right=0.997, top=1 - 0.5 / H, bottom=0.85 / H, wspace=0.03,
                                      hspace=0.38), S, FS, FS_HEAD, "figS3_edges_lrps_k", out, (0.55, 0.3 / H),
              [0.06, 0.48 / H, 0.32, 0.08 / H])


# ---- Section 4.4: factor clusters ------------------------------------------------------------------------------------

ROW, TOP = 0.92, 0.62  # one row of brains, space above the brain grid (in)


def _clusters(keep=4):
    data = {}
    for ds in DATASETS:
        roi = pd.read_csv(CLU / ds / "cluster_roi_coefficients_allm.csv")
        by = {c: g.reset_index(drop=True) for c, g in roi.groupby("cluster")}
        data[ds] = (by, np.stack([by[c].v_bar.to_numpy() for c in sorted(by)]))
    matched = {ds: match_to_reference(data[ds][1], data[DATASETS[0]][1], keep) for ds in DATASETS}
    return data, matched


def brain_grid(lay, keep, matched, data, coords, field, vmax, ring_fn=None, title_fn=None):
    """Rows = matched clusters, columns = dataset x (left sagittal, axial); returns (ScalarMappable, bottom in)."""
    from matplotlib.lines import Line2D

    fig = lay.fig
    gs = fig.add_gridspec(keep, 2 * len(DATASETS), left=0.17 / W, right=1 - 0.02 / W, top=1 - TOP / lay.H,
                          bottom=1 - (TOP + ROW * keep) / lay.H, wspace=0.0, hspace=0.24,
                          width_ratios=[1.3, 0.85] * len(DATASETS))
    sm, grid = None, {}
    for j in range(keep):
        for d, ds in enumerate(DATASETS):
            c, sgn = matched[ds][j]
            g = data[ds][0][c]
            vals = sgn * np.nan_to_num(g[field].to_numpy())
            axes = [fig.add_subplot(gs[j, 2 * d]), fig.add_subplot(gs[j, 2 * d + 1])]
            grid[j, d] = axes
            sm = plot_nodes(axes, coords, vals, (2 + 30 * np.abs(vals) / vmax) * S ** 2, vmax,
                            ring=ring_fn(g) if ring_fn else None, views="lz", outline_scale=0.55 * S, edge_scale=S)
            if title_fn:
                axes[0].text(0.0, 1.03, title_fn(ds, j, c), transform=axes[0].transAxes, fontsize=FS, color=INK2,
                             ha="left", va="bottom")
            if d == 0:
                p = axes[0].get_position()
                fig.text(p.x0 - 0.03 / W, p.y0 + p.height / 2, f"{j + 1}", fontsize=FS, color=INK, fontweight="bold",
                         ha="right", va="center")
    for d, ds in enumerate(DATASETS):
        p = grid[0, d][0].get_position()
        fig.text(p.x0, p.y1 + 0.19 / lay.H, ds, fontsize=8.5, fontweight="bold", color=DATASET_COLOR[ds],
                 ha="left", va="bottom")
    y_top = grid[0, 0][0].get_position().y1 + 0.38 / lay.H
    y_bot = grid[keep - 1, 0][0].get_position().y0
    for d in range(1, len(DATASETS)):
        x = (grid[0, d - 1][1].get_position().x1 + grid[0, d][0].get_position().x0) / 2
        fig.add_artist(Line2D([x, x], [y_bot, y_top], transform=fig.transFigure, color=SEP, lw=0.8 * S))
    return sm, TOP + ROW * keep


def colorbar(lay, sm, y_top, width, label, label_right=False):
    x0 = (W - width) / 2
    cb = lay.fig.colorbar(sm, cax=lay.ax(x0, y_top, width, 0.08), orientation="horizontal")
    if label_right:
        lay.text(x0 + width + 0.1, y_top + 0.04, label, fontsize=FS, color=INK2, ha="left", va="center")
    else:
        cb.set_label(label, fontsize=FS, color=INK2, labelpad=2)
    cb.ax.tick_params(labelsize=FS, colors=INK2, width=0.8 * S, length=3 * S)
    cb.outline.set_linewidth(0.8 * S)


def fig9(plt, out, keep=4):
    data, matched = _clusters(keep)
    coords, _ = aal_coords()
    memb = {ds: pd.read_csv(CLU / ds / "cluster_membership.csv") for ds in DATASETS}
    lay = Lay(plt, 7.88)
    fig, ref = lay.fig, DATASETS[0]
    vmax = max(np.abs(data[ds][1]).max() for ds in DATASETS)

    def title(ds, j, c):
        return f"c{c}" if ds == ref else f"c{c}, |cos| = {abs(data[ds][1][c - 1] @ data[ref][1][j]):.2f}"

    lay.panel(0.0, 0.2, "A", "Population networks of the leading clusters")
    sm, y = brain_grid(lay, keep, matched, data, coords, "v_bar", vmax, title_fn=title)
    colorbar(lay, sm, y + 0.16, 2.2, "centroid loading", label_right=True)

    # B) principal angles between the matched centroids
    labels, vecs, groups = [], [], []
    for ds in DATASETS:
        for j in range(keep):
            c, _ = matched[ds][j]
            labels.append(f"c{c}")
            vecs.append(data[ds][1][c - 1])
            groups.append(ds)
    Vs = np.stack(vecs)
    ang = np.degrees(np.arccos(np.clip(np.abs(Vs @ Vs.T), 0, 1)))
    yb = y + 0.68
    lay.panel(0.0, yb, "B", "Angles between centroids")
    side, x0, y0 = 2.25, 0.42, yb + 0.2
    axB = lay.ax(x0, y0, side, side)
    im = axB.imshow(ang, cmap="Blues_r", vmin=0, vmax=90)
    for a in range(len(labels)):
        for b in range(len(labels)):
            axB.text(b, a, f"{ang[a, b]:.0f}", ha="center", va="center", fontsize=FS,
                     color="white" if ang[a, b] < 40 else INK)
    axB.set_xticks(range(len(labels)), labels, rotation=90)
    axB.set_yticks(range(len(labels)), labels)
    axB.tick_params(length=0, pad=1.5, labelsize=FS)
    for tl, g in zip(axB.get_xticklabels(), groups):
        tl.set_color(DATASET_COLOR[g])
    for tl, g in zip(axB.get_yticklabels(), groups):
        tl.set_color(DATASET_COLOR[g])
    for b in range(1, len(DATASETS)):
        axB.axhline(b * keep - 0.5, color="white", lw=2.0 * S)
        axB.axvline(b * keep - 0.5, color="white", lw=2.0 * S)
    for d, ds in enumerate(DATASETS):
        axB.text(-1.75, d * keep + (keep - 1) / 2, ds, rotation=90, ha="center", va="center", fontsize=FS,
                 fontweight="bold", color=DATASET_COLOR[ds], clip_on=False)
    for sp in axB.spines.values():
        sp.set_visible(False)
    cb2 = fig.colorbar(im, cax=lay.ax(x0 + side + 0.07, y0, 0.08, side))
    cb2.set_label("principal angle (°)", fontsize=FS, color=INK2, labelpad=2)
    cb2.ax.tick_params(labelsize=FS, colors=INK2, width=0.8 * S, length=3 * S)
    cb2.outline.set_linewidth(0.8 * S)

    # C) rank of the components in each matched cluster
    lay.panel(3.3, yb, "C", "Component rank")
    k = int(max(m.component.max() for m in memb.values()))
    fr = {}
    for j in range(keep):
        for ds in DATASETS:
            r = memb[ds][memb[ds].cluster == matched[ds][j][0]].component.to_numpy()
            fr[j, ds] = np.bincount(r, minlength=k + 1)[1:] / r.size
    ymax = max(v.max() for v in fr.values()) * 1.08
    xs, ys, wc, hc = (3.72, 4.44), (y0 + 0.17, y0 + 1.27), 0.66, 0.72
    for j in range(keep):
        ax = lay.ax(xs[j % 2], ys[j // 2], wc, hc)
        for ds in DATASETS:
            ax.step(np.arange(1, k + 1), fr[j, ds], where="mid", color=DATASET_COLOR[ds], lw=1.1 * S)
        ax.set_xlim(0.5, k + 0.5)
        ax.set_ylim(0, ymax)
        ax.set_xticks([1, 5, 10])
        ax.set_yticks([0, 0.4, 0.8])
        ax.set_title(f"cluster {j + 1}", fontsize=FS, color=INK2, pad=2)
        style(ax, scale=S, color=AXIS)
        ax.tick_params(colors=INK2)
        if j % 2:
            ax.tick_params(axis="y", labelleft=False)
        if j < 2:
            ax.tick_params(axis="x", labelbottom=False)
    lay.text(3.44, ys[0] + (ys[1] + hc - ys[0]) / 2, "share of participants", fontsize=FS, color=INK2, rotation=90,
             ha="center", va="center")
    lay.text(xs[0] + (xs[1] + wc - xs[0]) / 2, ys[1] + hc + 0.27, "rank r of the component", fontsize=FS,
             color=INK2, ha="center", va="top")
    save(fig, out, "fig9_clusters")


def fig10(plt, out, keep=4, alpha=0.05):
    from matplotlib.patches import Patch

    data, matched = _clusters(keep)
    coords, _ = aal_coords()
    shares = pd.read_csv(CLU / "drive_varshare_by_subject.csv")
    lay = Lay(plt, 7.55)
    fig = lay.fig
    vmax = max(np.nanmax(np.abs(data[ds][0][matched[ds][j][0]].t)) for ds in DATASETS for j in range(keep))

    def title(ds, j, c):
        return f"c{c}: {int((data[ds][0][c].p_BH < alpha).sum())} significant"

    lay.panel(0.0, 0.2, "A", "Driving patterns: lagged ROIs predicting each network")
    sm, y = brain_grid(lay, keep, matched, data, coords, "t", vmax, ring_fn=lambda g: (g.p_BH < alpha).to_numpy(),
                       title_fn=title)
    colorbar(lay, sm, y + 0.22, 2.6, f"t of the driving coefficient (ring: BH-significant, FDR {alpha:.0%})")

    yb = y + 0.98
    lay.panel(0.0, yb, "B", "Variance predicted: within vs off network")
    ya, ha = yb + 0.2, 1.35
    axB = lay.ax(0.66, ya, 4.46, ha)
    handles = ([Patch(color=DATASET_COLOR[ds], label=ds) for ds in DATASETS]
               + [Patch(color="#6f6e69", label="within-network"), Patch(color="#cfcdc8", label="off-network")])
    fig.legend(handles=handles, loc="center", bbox_to_anchor=((0.66 + 4.46 / 2) / W, 1 - (ya + ha + 0.47) / lay.H),
               ncol=5, fontsize=FS, frameon=False, handlelength=1.1, handletextpad=0.4, columnspacing=0.9,
               borderaxespad=0)
    x = np.arange(1, keep + 1)
    w = 0.84 / (2 * len(DATASETS))
    for d, ds in enumerate(DATASETS):
        cis = [(mean_ci(s.R2_within), mean_ci(s.R2_off)) for s in
               (shares[(shares.dataset == ds) & (shares.cluster == matched[ds][j][0])] for j in range(keep))]
        xs = x - 0.42 + w * (2 * d) + w / 2
        ekw = dict(ecolor=INK2, elinewidth=0.8 * S, capsize=1.2, capthick=0.8 * S)
        for off, idx, col in ((0, 0, DATASET_COLOR[ds]), (w, 1, DATASET_TINT[ds])):
            m = np.array([c[idx][0] for c in cis])
            err = [m - [c[idx][1] for c in cis], [c[idx][2] for c in cis] - m]
            axB.bar(xs + off, m, w * 0.9, color=col, yerr=err, error_kw=ekw)
    axB.set_xticks(x, [f"{j}" for j in x])
    axB.set_xlabel("matched cluster", fontsize=FS, color=INK2, labelpad=2)
    style(axB, scale=S, color=AXIS)
    axB.tick_params(colors=INK2)
    axB.set_ylabel("variance of network\npredicted (in-sample)", fontsize=FS, color=INK2)
    save(fig, out, "fig10_driving")


FIGURES = {
    "01": {"fig1": fig1, "fig2": fig2, "fig3": lambda plt, out: fig_rank(plt, out, ("1B", "2"), "fig3_rank_1B_2"),
           "figS1": lambda plt, out: fig_rank(plt, out, ("1A", "1B"), "figS1_rank_1A_1B"), "figS2": figS2},
    "02": {"fig4": fig4, "fig5": fig5},
    "03": {"fig6": fig6, "fig7": fig7},
    "04": {"fig8": fig8, "figS3": figS3},
    "05": {"fig9": fig9, "fig10": fig10},
}


def make_figures(step: str, out: Path | None = None, only=None):
    """Draw the figures of one analysis step ('01'..'05'), or of all steps if step is None; `only` restricts to
    the named figures (e.g. ['fig8'])."""
    out = Path(out) if out else RESULTS_DIR / "figures"
    out.mkdir(parents=True, exist_ok=True)
    plt = setup()
    for s in ([step] if step else FIGURES):
        for name, fn in FIGURES[s].items():
            if only is None or name in only:
                fn(plt, out)
