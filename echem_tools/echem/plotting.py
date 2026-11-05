"""Shared matplotlib style: fixed categorical palette (validated, CVD-safe order), thin marks."""
from __future__ import annotations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cycler import cycler

PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]
STATUS = {"good": "#0ca30c", "warning": "#fab219", "serious": "#ec835a", "critical": "#d03b3b"}
SURFACE = "#fcfcfb"


def apply_style():
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.prop_cycle": cycler(color=PALETTE),
        "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6,
        "axes.edgecolor": "#b5b4ae", "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 1.8, "lines.markersize": 5, "font.size": 10, "legend.frameon": False,
        "figure.dpi": 110, "savefig.dpi": 150, "savefig.bbox": "tight",
        "text.color": "#0b0b0b", "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e", "ytick.color": "#52514e",
    })


def color_for(i: int) -> str:
    return PALETTE[i % len(PALETTE)]


def sequential(n: int, cmap: str = "Blues"):
    """One-hue light->dark ramp for ordered series (cycle number, time)."""
    import numpy as np
    cm = plt.get_cmap(cmap)
    return [cm(v) for v in np.linspace(0.35, 0.95, max(n, 1))]


apply_style()
