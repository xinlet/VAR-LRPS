"""Glass-brain plots of AAL ROI vectors and directed edges (nilearn glass-brain outlines)."""

from __future__ import annotations

import numpy as np

from .config import ATLAS_DIR

DIVERGING = ["#2a78d6", "#f0efec", "#e34948"]  # blue <-> red, neutral grey midpoint

VIEWS = {  # direction: (coordinate indices (horizontal, vertical), hemisphere filter)
    "l": ((1, 2), lambda c: c[:, 0] <= 0),
    "y": ((0, 2), lambda c: np.ones(len(c), bool)),
    "r": ((1, 2), lambda c: c[:, 0] >= 0),
    "z": ((0, 1), lambda c: np.ones(len(c), bool)),
}


def aal_coords(atlas_dir=ATLAS_DIR) -> tuple[np.ndarray, list[str]]:
    """MNI centre of mass of each AAL (SPM12) parcel, in ROI_MNI_V4.txt order (= data columns)."""
    import nibabel as nib
    import pandas as pd

    img = nib.load(atlas_dir / "atlas" / "AAL.nii")
    lab = np.asarray(img.dataobj).astype(int)
    tab = pd.read_csv(atlas_dir / "ROI_MNI_V4.txt", sep="\t", header=None, names=["short", "name", "value"])
    coords = []
    for v in tab["value"]:
        ijk = np.argwhere(lab == v)
        xyz = img.affine @ np.c_[ijk, np.ones(len(ijk))].T
        coords.append(xyz[:3].mean(axis=1))
    return np.array(coords), tab["name"].tolist()


def draw_outline(ax, direction, scale: float = 1.0):
    """nilearn glass-brain outline for one view, line widths multiplied by `scale`; returns the axis bounds."""
    from nilearn.plotting.glass_brain import plot_brain_schematics

    n_before = len(ax.patches)
    bounds = plot_brain_schematics(ax, direction)
    for patch in ax.patches[n_before:]:
        patch.set_linewidth(patch.get_linewidth() * scale)
    return bounds


def _finish(ax, bounds, direction):
    ax.set_xlim(bounds[0], bounds[1])
    ax.set_ylim(bounds[2], bounds[3])
    ax.set_aspect("equal")
    if direction == "l":
        ax.invert_xaxis()
    ax.axis("off")


def plot_nodes(axes, coords, values, sizes, vmax, ring=None, views="lz", outline_scale: float = 1.0,
               edge_scale: float = 1.0):
    """One row of glass-brain views with nodes coloured by `values` (diverging, +-vmax); `ring` marks nodes
    with a dark outline. Returns the ScalarMappable."""
    import matplotlib.colors as mcolors
    from matplotlib.cm import ScalarMappable

    cmap = mcolors.LinearSegmentedColormap.from_list("div", DIVERGING)
    norm = mcolors.Normalize(-vmax, vmax)
    ring = np.zeros(len(values), bool) if ring is None else ring
    for ax, direction in zip(axes, views):
        (h, v), hemi = VIEWS[direction]
        bounds = draw_outline(ax, direction, outline_scale)
        sel = np.flatnonzero(hemi(coords))
        sel = sel[np.argsort(np.abs(values[sel]))]  # strongest on top
        ax.scatter(coords[sel, h], coords[sel, v], s=sizes[sel], c=values[sel], cmap=cmap, norm=norm,
                   edgecolors=np.where(ring[sel, None], [[0.04, 0.04, 0.04, 1.0]], [[1, 1, 1, 0.9]]),
                   linewidths=np.where(ring[sel], 1.4, 0.4) * edge_scale, zorder=3)
        _finish(ax, bounds, direction)
    return ScalarMappable(norm=norm, cmap=cmap)


def plot_edges(axes, coords, src, dst, tvals, thick, cmap, norm, views="yz", scale: float = 1.0,
               outline_scale: float = 0.55):
    """Directed edges src -> dst as arcs bending to the left of their direction with an arrowhead at the
    target, coloured by t; `thick` marks edges drawn thicker (e.g. Bonferroni-significant)."""
    from matplotlib.patches import FancyArrowPatch

    order = np.argsort(np.abs(tvals))  # strongest on top
    for ax, direction in zip(axes, views):
        (h, v), _ = VIEWS[direction]
        bounds = draw_outline(ax, direction, outline_scale * scale)
        ax.scatter(coords[:, h], coords[:, v], s=1.2 * scale, color="#8a8984", linewidths=0, zorder=2)
        for e in order:
            ax.add_patch(FancyArrowPatch(coords[src[e], [h, v]], coords[dst[e], [h, v]],
                                         connectionstyle="arc3,rad=0.25", arrowstyle="-|>",
                                         mutation_scale=(5 if thick[e] else 4) * scale,
                                         lw=(1.5 if thick[e] else 0.7) * scale, color=cmap(norm(tvals[e])),
                                         shrinkA=1, shrinkB=1, zorder=3))
        _finish(ax, bounds, direction)
