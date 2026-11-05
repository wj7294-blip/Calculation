from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ..plotting import color_for, MARKERS, STATUS


def plot_cell_curves(results: list, meta: pd.DataFrame, out: Path, color_by: str = "electrolyte"):
    """I-t (linear) and log-log current decay for all cells, coloured by one factor."""
    levels = list(pd.unique(meta[color_by].dropna())) if color_by in meta else []
    cmap = {l: color_for(i) for i, l in enumerate(levels)}
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.3))
    seen = set()
    for r in results:
        row = meta.loc[meta.cell_id == r.cell_id]
        lv = row[color_by].iloc[0] if len(row) and color_by in row else None
        c = cmap.get(lv, "#52514e")
        mass_g = r.metrics["active_mass_mg"] / 1000 if np.isfinite(r.metrics.get("active_mass_mg", np.nan)) else np.nan
        y = r.cv["current_mA"].to_numpy() / mass_g
        t = r.cv["t_cv_s"].to_numpy()
        lab = str(lv) if lv not in seen else None
        seen.add(lv)
        axes[0].plot(t / 3600, y, "-", lw=0.9, color=c, alpha=0.7, label=lab)
        m = (t > 0) & (y > 0)
        axes[1].loglog(t[m] / 3600, y[m], "-", lw=0.9, color=c, alpha=0.7)
    axes[0].set_xlabel("time in CV hold (h)"); axes[0].set_ylabel("current (mA/g)")
    axes[0].set_yscale("log"); axes[0].set_title("Current during 4.6 V hold", fontsize=10)
    axes[1].set_xlabel("time in CV hold (h)"); axes[1].set_ylabel("current (mA/g)")
    axes[1].set_title("log-log decay (slope = decay_slope)", fontsize=10)
    if levels:
        axes[0].legend(title=color_by, fontsize=8)
    fig.tight_layout(); out.parent.mkdir(parents=True, exist_ok=True); fig.savefig(out); plt.close(fig)
    return out


def plot_single_cell(r, out: Path, cfg_f: dict):
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    t = r.cv["t_cv_s"].to_numpy() / 3600; I = r.cv["current_mA"].to_numpy()
    axes[0].plot(t, I, color=color_for(0), lw=1.2)
    axes[0].axvspan(0, float(cfg_f.get("transient_h", 1.0)), color="#e4e3df", alpha=0.6, label="transient (excluded)")
    axes[0].axvspan(t[-1] - float(cfg_f.get("leak_window_h", 1.0)), t[-1], color=color_for(1), alpha=0.15, label="leak window")
    axes[0].set_yscale("log"); axes[0].set_xlabel("time in CV hold (h)"); axes[0].set_ylabel("current (mA)")
    axes[0].legend(fontsize=8); axes[0].set_title(f"{r.cell_id}: leak = {r.metrics['leak_mA_g']:.3f} mA/g", fontsize=10)
    m = (t > 0) & (I > 0)
    axes[1].loglog(t[m], I[m], ".", ms=2.5, color=color_for(0))
    if np.isfinite(r.metrics["decay_slope"]):
        w0, w1 = cfg_f.get("decay_fit_window_h", [2, 24])
        tt = np.logspace(np.log10(max(w0, t[m].min())), np.log10(min(w1, t[-1])), 50)
        axes[1].loglog(tt, 10 ** (r.metrics["decay_intercept"] + r.metrics["decay_slope"] * np.log10(tt * 3600)),
                       "--", color=color_for(1), label=f"slope = {r.metrics['decay_slope']:.2f}")
        axes[1].legend(fontsize=8)
    axes[1].set_xlabel("time (h)"); axes[1].set_ylabel("current (mA)"); axes[1].set_title("log-log decay", fontsize=10)
    fig.tight_layout(); out.parent.mkdir(parents=True, exist_ok=True); fig.savefig(out); plt.close(fig)
    return out


def plot_group_bars(cells: pd.DataFrame, metric: str, fA: str, fB: str, outliers: pd.Series, out: Path, ylabel: str = ""):
    """mean ± SD bars per (A,B) with individual cell points; outliers ringed (flagged, not removed)."""
    A = list(pd.unique(cells[fA])); B = list(pd.unique(cells[fB]))
    nB = len(B); width = 0.8 / nB
    fig, ax = plt.subplots(figsize=(1.6 * len(A) + 3, 4.2))
    rng = np.random.default_rng(0)
    for j, b in enumerate(B):
        c = color_for(j)
        xs, means, sds = [], [], []
        for i, a in enumerate(A):
            g = cells[(cells[fA] == a) & (cells[fB] == b)]
            x0 = i - 0.4 + width * (j + 0.5)
            y = g[metric].astype(float)
            xs.append(x0); means.append(y.mean()); sds.append(y.std(ddof=1) if len(y) > 1 else 0)
            jit = (rng.random(len(y)) - 0.5) * width * 0.6
            ax.scatter(x0 + jit, y, s=16, color=c, edgecolor="white", linewidth=0.6, zorder=3)
            o = outliers.loc[g.index] if outliers is not None else pd.Series(False, index=g.index)
            if o.any():
                ax.scatter(x0 + jit[o.values], y[o.values], s=70, facecolor="none", edgecolor=STATUS["critical"], linewidth=1.2, zorder=4)
        ax.bar(xs, means, width * 0.9, color=c, alpha=0.28, label=str(b), zorder=1)
        ax.errorbar(xs, means, yerr=sds, fmt="none", ecolor=c, elinewidth=1.2, capsize=3, zorder=2)
    ax.set_xticks(range(len(A))); ax.set_xticklabels(A); ax.set_xlabel(fA)
    ax.set_ylabel(ylabel or metric); ax.legend(title=fB, fontsize=8)
    ax.set_title(f"{metric}: mean ± SD, points = cells, red ring = outlier (flagged only)", fontsize=9.5)
    fig.tight_layout(); out.parent.mkdir(parents=True, exist_ok=True); fig.savefig(out); plt.close(fig)
    return out


def plot_interaction(cells: pd.DataFrame, metric: str, fA: str, fB: str, out: Path, ylabel: str = "", pvals: dict | None = None):
    """x = fA (surface), one line per fB (electrolyte), mean ± SD. Non-parallel lines = interaction."""
    A = list(pd.unique(cells[fA])); B = list(pd.unique(cells[fB]))
    fig, ax = plt.subplots(figsize=(1.3 * len(A) + 3.5, 4.2))
    for j, b in enumerate(B):
        g = cells[cells[fB] == b].groupby(fA, observed=True)[metric].agg(["mean", "std"]).reindex(A)
        x = np.arange(len(A)) + (j - (len(B) - 1) / 2) * 0.06
        ax.errorbar(x, g["mean"], yerr=g["std"], fmt="-" + MARKERS[j % len(MARKERS)], color=color_for(j),
                    capsize=3, ms=6, lw=1.8, label=str(b))
    ax.set_xticks(range(len(A))); ax.set_xticklabels(A); ax.set_xlabel(fA); ax.set_ylabel(ylabel or metric)
    ttl = f"Interaction plot: {metric}"
    if pvals:
        ttl += "   p: " + ", ".join(f"{k}={v:.3g}" for k, v in pvals.items())
    ax.set_title(ttl, fontsize=9.5); ax.legend(title=fB, fontsize=8)
    fig.tight_layout(); out.parent.mkdir(parents=True, exist_ok=True); fig.savefig(out); plt.close(fig)
    return out
