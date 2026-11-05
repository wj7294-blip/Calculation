"""EIS preprocessing: frequency cropping, inductive-tail handling, lin-KK validation."""
from __future__ import annotations
import warnings
from dataclasses import dataclass, field
import numpy as np

from ..io.loaders import EISData


@dataclass
class KKResult:
    M: int = 0
    mu: float = float("nan")
    res_real: np.ndarray | None = None   # fractional residuals (Z_meas - Z_kk)/|Z|
    res_imag: np.ndarray | None = None
    max_res_pct: float = float("nan")
    rms_res_pct: float = float("nan")
    warn: bool = False
    message: str = ""
    Z_fit: np.ndarray | None = None


def preprocess_eis(d: EISData, cfg_eis: dict) -> tuple[EISData, dict]:
    """Crop frequency range and handle inductive tail. Returns (data, info)."""
    info = {"n_in": len(d.freq)}
    f, Z = d.freq.copy(), d.Z.copy()
    fmin = cfg_eis.get("freq_min_Hz")
    fmax = cfg_eis.get("freq_max_Hz")
    m = np.ones_like(f, bool)
    if fmin is not None:
        m &= f >= float(fmin)
    if fmax is not None:
        m &= f <= float(fmax)
    # finite only
    m &= np.isfinite(Z.real) & np.isfinite(Z.imag)
    mode = str(cfg_eis.get("inductive_tail", "fit")).lower()
    n_ind = int(np.sum(Z.imag > 0))
    info["n_inductive_points"] = n_ind
    if mode == "drop":
        m &= ~(Z.imag > 0)
    info["inductive_mode"] = mode
    # remove exact duplicate frequencies (keep first)
    _, first_idx = np.unique(np.round(np.log10(f[m]), 6), return_index=True)
    keep_idx = np.where(m)[0][np.sort(first_idx)]
    sel = np.zeros_like(m)
    sel[keep_idx] = True
    extra = d.extra.iloc[sel].reset_index(drop=True) if d.extra is not None else None
    out = EISData(f[sel], Z[sel], d.source, extra, dict(d.meta)).sort()
    info["n_out"] = len(out.freq)
    return out, info


def linkk_check(d: EISData, cfg_kk: dict) -> KKResult:
    """Kramers-Kronig validation via impedance.py lin-KK (Schönleber et al. 2014)."""
    if not cfg_kk.get("enabled", True):
        return KKResult(message="disabled")
    try:
        from impedance.validation import linKK
        # impedance.py <=1.7.1 evaluates circuit strings with a namespace lacking numpy on
        # recent Python versions (NameError: np). Inject it; harmless otherwise.
        from impedance.models.circuits import elements as _el
        _el.circuit_elements.setdefault("np", np)
    except Exception as e:  # pragma: no cover
        return KKResult(message=f"impedance.py not available: {e}")
    f, Z = d.freq, d.Z
    if len(f) < 8:
        return KKResult(message="too few points for lin-KK")
    c = float(cfg_kk.get("c", 0.85))
    max_M = int(min(int(cfg_kk.get("max_M", 100)), max(3, len(f) - 2)))
    try:
        from impedance.validation import get_tc_distribution, fit_linKK, eval_linKK, residuals_linKK
        # Schönleber loop (M = 3 .. max_M) using impedance.py building blocks. impedance.linKK
        # stops at the *first* M with mu < c, which can be a badly under-fitted M (e.g. M=4 with
        # 50 % residuals on two close arcs); we add a residual-plateau safeguard:
        # choose the smallest M with mu < c whose rms residual is within `plateau_tol` of the best.
        rows = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for M in range(3, max_M + 1):
                ts = get_tc_distribution(f, M)
                el, mu = fit_linKK(f, ts, M, Z, "complex", True)
                Zf = eval_linKK(el, ts, f)
                rms = float(np.sqrt(np.mean(np.abs(Z - Zf) ** 2 / np.abs(Z) ** 2)))
                rows.append((M, float(mu), rms, el, ts))
        rms_min = min(r[2] for r in rows)
        tol = float(cfg_kk.get("plateau_tol", 1.2))
        cand = [r for r in rows if r[1] < c and r[2] <= tol * rms_min]
        if not cand:
            cand = [r for r in rows if r[1] < c] or [min(rows, key=lambda r: r[2])]
        M, mu, _, el, ts = cand[0]
        Z_fit = eval_linKK(el, ts, f)
        res_re = residuals_linKK(el, ts, Z, f, residuals="real")
        res_im = residuals_linKK(el, ts, Z, f, residuals="imag")
    except Exception as e:
        return KKResult(message=f"lin-KK failed: {e}")
    res_re = np.asarray(res_re, float)
    res_im = np.asarray(res_im, float)
    allres = np.concatenate([res_re, res_im]) * 100.0
    max_pct = float(np.nanmax(np.abs(allres)))
    rms_pct = float(np.sqrt(np.nanmean(allres ** 2)))
    thr_rms = float(cfg_kk.get("residual_warn_rms_pct", 1.0))
    thr_max = float(cfg_kk.get("residual_warn_max_pct", 5.0))
    warn = (rms_pct > thr_rms) or (max_pct > thr_max)
    msg = (f"lin-KK: M={M}, mu={mu:.3f}, rms|res|={rms_pct:.2f}%, max|res|={max_pct:.2f}%"
           + (f"  ** WARNING: rms>{thr_rms}% or max>{thr_max}% -> data may violate KK "
              f"(drift / nonlinearity / noise); interpret DRT with care" if warn else ""))
    return KKResult(int(M), float(mu), res_re, res_im, max_pct, rms_pct, warn, msg, np.asarray(Z_fit))
