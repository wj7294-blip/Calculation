from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from ..plotting import color_for, sequential, MARKERS, STATUS
from .drt import DRTResult


def _colors(n, ordered):
    return sequential(n) if ordered else [color_for(i) for i in range(n)]


def plot_nyquist_bode_drt(results: list[DRTResult], labels: list[str], out: Path, title: str = "",
                          ordered: bool = False, kk_notes: list[str] | None = None):
    """3-panel figure: Nyquist (data + fit), Bode (|Z|, -phase), DRT γ vs τ. Works for 1..N spectra."""
    cols = _colors(len(results), ordered)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    axN, axB, axD = axes
    axB2 = axB.twinx()  # phase on the same frequency axis (same x, not a dual *y* of one quantity)
    axB2.grid(False)
    for i, (r, lab) in enumerate(zip(results, labels)):
        c = cols[i]
        m = MARKERS[i % len(MARKERS)]
        axN.plot(r.Z.real, -r.Z.imag, m, ms=4, mfc="none", color=c, label=lab)
        axN.plot(r.Z_fit.real, -r.Z_fit.imag, "-", lw=1.2, color=c, alpha=0.8)
        axB.plot(r.freq, np.abs(r.Z), m, ms=3.5, mfc="none", color=c)
        axB.plot(r.freq, np.abs(r.Z_fit), "-", lw=1.0, color=c, alpha=0.8)
        ph = -np.degrees(np.angle(r.Z))
        axB2.plot(r.freq, ph, m, ms=3.5, color=c, alpha=0.5)
        axD.plot(r.tau, r.gamma, "-", color=c, label=f"{lab}  (λ={r.lam:.1e})")
        if r.peaks is not None and len(results) == 1:
            for _, p in r.peaks.iterrows():
                axD.axvline(p.tau_s, color=c, ls=":", lw=0.8)
                axD.annotate(f"τ={p.tau_s:.2e}s\nR={p.R_ohm:.2f}Ω", (p.tau_s, p.gamma_max_ohm),
                             textcoords="offset points", xytext=(4, 4), fontsize=7.5, color="#52514e")
    axN.set_xlabel("Z' (Ω)"); axN.set_ylabel("-Z'' (Ω)"); axN.set_aspect("equal", adjustable="datalim")
    axN.set_title("Nyquist (markers: data, lines: DRT fit)", fontsize=10)
    axB.set_xscale("log"); axB.set_yscale("log"); axB.set_xlabel("Frequency (Hz)"); axB.set_ylabel("|Z| (Ω)")
    axB2.set_ylabel("-phase (°)", color="#52514e"); axB.set_title("Bode", fontsize=10)
    axD.set_xscale("log"); axD.set_xlabel("τ (s)"); axD.set_ylabel("γ(ln τ) (Ω)"); axD.set_title("DRT", fontsize=10)
    secax = axD.secondary_xaxis("top", functions=(lambda t: 1 / (2 * np.pi * np.maximum(t, 1e-300)),
                                                   lambda f: 1 / (2 * np.pi * np.maximum(f, 1e-300))))
    secax.set_xlabel("f = 1/(2πτ) (Hz)", fontsize=8)
    if len(results) > 1 or True:
        axD.legend(fontsize=7.5, loc="best")
    if kk_notes:
        fig.text(0.01, -0.02, " | ".join(kk_notes)[:300], fontsize=7.5,
                 color=STATUS["serious"] if any("WARNING" in k for k in kk_notes) else "#52514e")
    if title:
        fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return out


def plot_gcv(r: DRTResult, out: Path, title: str = ""):
    if r.gcv_lams is None:
        return None
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    ax[0].loglog(r.gcv_lams, r.gcv_scores, "o-", ms=3)
    ax[0].axvline(r.lam, color=color_for(1), ls="--", label=f"λ* = {r.lam:.2e}")
    ax[0].set_xlabel("λ"); ax[0].set_ylabel("GCV score"); ax[0].legend(); ax[0].set_title("GCV", fontsize=10)
    if r.lcurve is not None:
        ax[1].loglog(r.lcurve[:, 1], r.lcurve[:, 2], "o-", ms=3)
        k = int(np.argmin(np.abs(np.log(r.lcurve[:, 0]) - np.log(r.lam))))
        ax[1].plot(r.lcurve[k, 1], r.lcurve[k, 2], "o", ms=9, mfc="none", color=color_for(1))
        ax[1].set_xlabel("||residual||"); ax[1].set_ylabel("||U x|| (roughness)"); ax[1].set_title("L-curve", fontsize=10)
    if title:
        fig.suptitle(title, fontsize=10)
    fig.tight_layout(); out.parent.mkdir(parents=True, exist_ok=True); fig.savefig(out); plt.close(fig)
    return out


def plot_lambda_sweep(results: list[DRTResult], lam_auto: float, truth, out: Path, title: str = ""):
    """γ(τ) for several λ, stacked, with true peak positions. Shows splitting / merging."""
    n = len(results)
    fig, axes = plt.subplots(n, 1, figsize=(7.5, 1.7 * n + 1), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, r in zip(axes, results):
        is_auto = np.isclose(r.lam, lam_auto, rtol=0.05)
        c = color_for(1) if is_auto else color_for(0)
        ax.plot(r.tau, r.gamma, color=c, lw=1.6)
        ax.fill_between(r.tau, 0, r.gamma, color=c, alpha=0.12)
        npk = 0 if r.peaks is None else len(r.peaks)
        ax.text(0.01, 0.80, f"λ = {r.lam:.0e}" + ("  (GCV)" if is_auto else "") + f"   peaks: {npk}",
                transform=ax.transAxes, fontsize=8.5)
        if truth is not None:
            for _, t in truth.iterrows():
                ax.axvline(t.tau_true_s, color="#52514e", ls=":", lw=0.9)
        ax.set_ylabel("γ (Ω)", fontsize=8)
        ax.set_xscale("log")
    axes[-1].set_xlabel("τ (s)   (dotted: true τ)")
    if title:
        fig.suptitle(title, fontsize=10)
    fig.tight_layout(); out.parent.mkdir(parents=True, exist_ok=True); fig.savefig(out); plt.close(fig)
    return out


def plot_kk_residuals(freq, kk, out: Path, title: str = ""):
    if kk.res_real is None:
        return None
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.semilogx(freq, kk.res_real * 100, "o", ms=3, label="Δ real")
    ax.semilogx(freq, kk.res_imag * 100, "s", ms=3, label="Δ imag")
    ax.axhline(0, color="#52514e", lw=0.6)
    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("lin-KK residual (%)"); ax.legend()
    ax.set_title(title or kk.message.split("**")[0], fontsize=9)
    fig.tight_layout(); out.parent.mkdir(parents=True, exist_ok=True); fig.savefig(out); plt.close(fig)
    return out
