"""Distribution of Relaxation Times (DRT).

Model (γ defined per ln τ):
    Z(ω) = R_inf + jωL + ∫ γ(ln τ) / (1 + jωτ) d ln τ ,   γ(ln τ) = Σ_m x_m φ_m(ln τ)
φ_m are Gaussian radial basis functions centred on a log-spaced τ grid.
Solve  min_x ||W(A x - b)||² + λ xᵀ M x   s.t. x ≥ 0  (NNLS on an augmented system).
M is the Gram matrix of the k-th derivative of the basis (k = reg_order, Tikhonov).
λ: generalized cross-validation (GCV) on the active set, or a user value.

If pyDRTtools can be imported and backend == 'pydrttools' (or 'auto'), it is tried first;
on any failure the built-in solver is used (a note is stored in result.backend).
"""
from __future__ import annotations
import warnings
from dataclasses import dataclass, field
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.linalg import cholesky
from scipy.signal import find_peaks

from ..io.loaders import EISData


@dataclass
class DRTResult:
    tau: np.ndarray                 # display grid (fine)
    gamma: np.ndarray               # γ(ln τ) [ohm]
    tau_nodes: np.ndarray           # basis centres
    x: np.ndarray                   # basis coefficients
    R_inf: float
    L: float
    lam: float
    lam_auto: bool
    R_pol: float                    # ∫ γ d ln τ
    freq: np.ndarray
    Z: np.ndarray
    Z_fit: np.ndarray
    resid_pct: float                # rms of |Z-Zfit|/|Z| in %
    gcv_lams: np.ndarray | None = None
    gcv_scores: np.ndarray | None = None
    lcurve: np.ndarray | None = None   # columns: lam, ||r||, ||Ux||
    peaks: pd.DataFrame | None = None
    backend: str = "builtin"
    label: str = ""

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame({"tau_s": self.tau, "freq_Hz": 1.0 / (2 * np.pi * self.tau), "gamma_ohm": self.gamma})


# ----------------------------------------------------------------------------
# basis & matrices
# ----------------------------------------------------------------------------
def _gauss_basis(y: np.ndarray, centres: np.ndarray, eps: float, order: int = 0) -> np.ndarray:
    """Gaussian RBF φ_m(y) = exp(-(eps (y - y_m))²) and its derivatives. Returns [len(y), len(centres)]."""
    d = y[:, None] - centres[None, :]
    g = np.exp(-(eps * d) ** 2)
    if order == 0:
        return g
    if order == 1:
        return -2 * eps ** 2 * d * g
    if order == 2:
        return (4 * eps ** 4 * d ** 2 - 2 * eps ** 2) * g
    raise ValueError("reg_order must be 0, 1 or 2")


def build_system(freq: np.ndarray, tau_nodes: np.ndarray, eps: float, reg_order: int,
                 fine_per_node: int = 12):
    """Return A_re, A_im [n_f, M], U (M = UᵀU), fine grid (y_fine, w_fine, Phi_fine)."""
    y_nodes = np.log(tau_nodes)
    dy = np.mean(np.diff(y_nodes))
    pad = 4.0 / eps
    n_fine = int(np.ceil((y_nodes[-1] - y_nodes[0] + 2 * pad) / dy * fine_per_node)) + 1
    y_fine = np.linspace(y_nodes[0] - pad, y_nodes[-1] + pad, n_fine)
    w = np.full(n_fine, y_fine[1] - y_fine[0])
    w[0] *= 0.5
    w[-1] *= 0.5
    Phi = _gauss_basis(y_fine, y_nodes, eps, 0)
    omega = 2 * np.pi * freq
    K = 1.0 / (1.0 + 1j * omega[:, None] * np.exp(y_fine)[None, :])     # [n_f, n_fine]
    A = K @ (Phi * w[:, None])                                          # [n_f, M] complex
    Phid = _gauss_basis(y_fine, y_nodes, eps, reg_order)
    Mreg = (Phid * w[:, None]).T @ Phid
    Mreg = 0.5 * (Mreg + Mreg.T)
    jitter = 1e-10 * np.trace(Mreg) / len(Mreg)
    U = cholesky(Mreg + jitter * np.eye(len(Mreg)), lower=False)
    return A.real, A.imag, U, (y_fine, w, Phi)


def _normalise_penalty(U, A_re, A_im, wts):
    """Scale U so that λ is dimensionless: tr(UᵀU) = tr(AᵀA) (weighted). Then λ ~ 1 means
    data misfit and smoothness are weighted equally; GCV optimum typically lies in 1e-4..1e-1."""
    Aw = np.vstack([A_re * wts[:, None], A_im * wts[:, None]])
    sA = np.sum(Aw * Aw)
    sU = np.sum(U * U)
    return U * np.sqrt(sA / sU)


# ----------------------------------------------------------------------------
# solver
# ----------------------------------------------------------------------------
def _solve(A_re, A_im, b_re, b_im, U, lam, wts, fit_L, omega):
    """NNLS of [R_inf, L, x...] with Tikhonov penalty on x. Returns (theta, resid_vec, Afull)."""
    n_f, M = A_re.shape
    one = np.ones(n_f)
    zero = np.zeros(n_f)
    cols_re = [one, zero if fit_L else None, A_re]
    cols_im = [zero, omega if fit_L else None, A_im]
    top_re = np.column_stack([c for c in cols_re if c is not None])
    top_im = np.column_stack([c for c in cols_im if c is not None])
    Afull = np.vstack([top_re * wts[:, None], top_im * wts[:, None]])
    b = np.concatenate([b_re * wts, b_im * wts])
    n_extra = 2 if fit_L else 1
    pen = np.hstack([np.zeros((M, n_extra)), np.sqrt(lam) * U])
    Aaug = np.vstack([Afull, pen])
    baug = np.concatenate([b, np.zeros(M)])
    theta, _ = nnls(Aaug, baug, maxiter=50 * Aaug.shape[1])
    r = Afull @ theta - b
    return theta, r, Afull, n_extra


def _gcv_curve(Afull, b, U, lams, n_extra):
    """Classical GCV on the *unconstrained* Tikhonov problem (the NNLS active set is too small
    to give a meaningful effective number of parameters; the unconstrained ridge is the
    standard choice, as in pyDRTtools).
        GCV(λ) = n ||A x_λ - b||² / (n - tr H_λ)²,  H_λ = A (AᵀA + λ M̃)⁻¹ Aᵀ
    M̃ = blockdiag(0_{n_extra}, UᵀU). Uses the GSVD-free direct form (small matrices)."""
    n, p = Afull.shape
    AtA = Afull.T @ Afull
    Atb = Afull.T @ b
    Mt = np.zeros((p, p))
    Mt[n_extra:, n_extra:] = U.T @ U
    scores = np.empty(len(lams))
    for i, l in enumerate(lams):
        G = AtA + l * Mt
        try:
            Ginv = np.linalg.inv(G)
        except np.linalg.LinAlgError:
            scores[i] = np.inf
            continue
        x = Ginv @ Atb
        r = Afull @ x - b
        trH = float(np.einsum("ij,ji->", AtA, Ginv))
        denom = n - trH
        scores[i] = n * float(r @ r) / denom ** 2 if denom > 0 else np.inf
    return scores


def compute_drt(d: EISData, cfg_drt: dict, label: str = "") -> DRTResult:
    backend = str(cfg_drt.get("backend", "auto")).lower()
    if backend in ("auto", "pydrttools"):
        res = _try_pydrttools(d, cfg_drt, label)
        if res is not None:
            return res
        if backend == "pydrttools":
            warnings.warn("pyDRTtools requested but not usable; falling back to builtin DRT")
    return compute_drt_builtin(d, cfg_drt, label)


def compute_drt_builtin(d: EISData, cfg_drt: dict, label: str = "", lam: float | None = None) -> DRTResult:
    freq = np.asarray(d.freq, float)
    Z = np.asarray(d.Z, complex)
    omega = 2 * np.pi * freq
    ppd = int(cfg_drt.get("tau_per_decade", 10))
    pad = float(cfg_drt.get("tau_pad_decades", 1.0))
    lt_min = np.log10(1 / (2 * np.pi * freq.max())) - pad
    lt_max = np.log10(1 / (2 * np.pi * freq.min())) + pad
    n_nodes = int(np.round((lt_max - lt_min) * ppd)) + 1
    tau_nodes = np.logspace(lt_min, lt_max, n_nodes)
    dy = np.log(10) / ppd
    shape = float(cfg_drt.get("rbf_shape", 0.5))
    # FWHM of basis = shape * (grid spacing) ; eps from FWHM: exp(-(eps*FWHM/2)^2) = 1/2
    fwhm = shape * dy
    eps = 2 * np.sqrt(np.log(2)) / fwhm
    reg_order = int(cfg_drt.get("reg_order", 2))
    A_re, A_im, U, fine = build_system(freq, tau_nodes, eps, reg_order)

    wmode = str(cfg_drt.get("weighting", "modulus")).lower()
    wts = 1.0 / np.abs(Z) if wmode == "modulus" else np.ones_like(freq)
    fit_L = bool(cfg_drt.get("fit_inductance", True))
    U = _normalise_penalty(U, A_re, A_im, wts)

    lam_cfg = cfg_drt.get("lambda", "auto") if lam is None else lam
    gcv_l = gcv_s = None
    lcurve = []
    lam_auto = False
    if isinstance(lam_cfg, str) and lam_cfg.lower() == "auto":
        lam_auto = True
        g = cfg_drt.get("lambda_grid", [1e-6, 10.0, 36])
        lams = np.logspace(np.log10(float(g[0])), np.log10(float(g[1])), int(g[2]))
        theta0, r0, Afull, n_extra = _solve(A_re, A_im, Z.real, Z.imag, U, lams[0], wts, fit_L, omega)
        b = np.concatenate([Z.real * wts, Z.imag * wts])
        scores = _gcv_curve(Afull, b, U, lams, n_extra)
        for l in lams:
            th, r, _, _ = _solve(A_re, A_im, Z.real, Z.imag, U, l, wts, fit_L, omega)
            lcurve.append([l, np.linalg.norm(r), np.linalg.norm(U @ th[n_extra:])])
        gcv_l, gcv_s = lams, scores
        ok = np.isfinite(scores)
        lam_best = float(lams[ok][np.argmin(scores[ok])]) if ok.any() else 1e-3
    else:
        lam_best = float(lam_cfg)
    theta, r, Afull, n_extra = _solve(A_re, A_im, Z.real, Z.imag, U, lam_best, wts, fit_L, omega)
    R_inf = float(theta[0])
    L = float(theta[1]) if fit_L else 0.0
    x = theta[n_extra:]
    y_fine, w, Phi = fine
    gamma_fine = Phi @ x
    # display grid: restrict to node span (padding of basis beyond nodes is numerical)
    msk = (y_fine >= np.log(tau_nodes[0])) & (y_fine <= np.log(tau_nodes[-1]))
    tau_disp = np.exp(y_fine[msk])
    gamma = gamma_fine[msk]
    R_pol = float(np.sum(x) * np.sqrt(np.pi) / eps)
    Z_fit = R_inf + 1j * omega * L + (A_re @ x) + 1j * (A_im @ x)
    resid_pct = float(np.sqrt(np.mean(np.abs(Z - Z_fit) ** 2 / np.abs(Z) ** 2)) * 100)
    res = DRTResult(tau_disp, gamma, tau_nodes, x, R_inf, L, lam_best, lam_auto, R_pol, freq, Z, Z_fit,
                    resid_pct, gcv_l, gcv_s, np.array(lcurve) if lcurve else None, None, "builtin", label)
    res.peaks = detect_peaks(res, cfg_drt)
    return res


# ----------------------------------------------------------------------------
# peaks
# ----------------------------------------------------------------------------
def detect_peaks(res: DRTResult, cfg_drt: dict) -> pd.DataFrame:
    g = res.gamma
    y = np.log(res.tau)
    if g.max() <= 0:
        return pd.DataFrame(columns=["peak", "tau_s", "freq_Hz", "R_ohm", "gamma_max_ohm", "frac_of_Rpol", "tau_lo_s", "tau_hi_s"])
    prom = float(cfg_drt.get("peak_prominence", 0.03)) * g.max()
    hmin = float(cfg_drt.get("peak_min_rel_height", 0.02)) * g.max()
    idx, props = find_peaks(g, prominence=prom, height=hmin)
    rows = []
    bounds = []
    for k, i in enumerate(idx):
        lo = 0 if k == 0 else int(idx[k - 1] + np.argmin(g[idx[k - 1]:i]))
        hi = len(g) - 1 if k == len(idx) - 1 else int(i + np.argmin(g[i:idx[k + 1]]))
        bounds.append((lo, hi))
    for k, (i, (lo, hi)) in enumerate(zip(idx, bounds)):
        R = float(np.trapezoid(g[lo:hi + 1], y[lo:hi + 1]))
        rows.append({"peak": k + 1, "tau_s": float(res.tau[i]), "freq_Hz": float(1 / (2 * np.pi * res.tau[i])),
                     "R_ohm": R, "gamma_max_ohm": float(g[i]),
                     "frac_of_Rpol": R / res.R_pol if res.R_pol > 0 else np.nan,
                     "tau_lo_s": float(res.tau[lo]), "tau_hi_s": float(res.tau[hi])})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# optional pyDRTtools backend (best effort)
# ----------------------------------------------------------------------------
def _try_pydrttools(d: EISData, cfg_drt: dict, label: str):
    import contextlib, io, os
    try:
        with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            warnings.simplefilter("ignore")
            from pyDRTtools.basics import EIS_object      # type: ignore
            from pyDRTtools.runs import simple_run         # type: ignore
    except BaseException:   # ImportError, but also SystemExit from GUI-less installs
        return None
    try:
        lam_cfg = cfg_drt.get("lambda", "auto")
        auto = isinstance(lam_cfg, str)
        data = EIS_object(np.asarray(d.freq, float), np.asarray(d.Z.real, float), np.asarray(d.Z.imag, float))
        out = simple_run(data, rbf_type="Gaussian", data_used="Combined Re-Im Data",
                         induct_used=1 if cfg_drt.get("fit_inductance", True) else 0,
                         der_used=f"{int(cfg_drt.get('reg_order', 2))}{'nd' if int(cfg_drt.get('reg_order', 2)) == 2 else 'st'} order",
                         cv_type="GCV" if auto else "none", reg_param=1e-3 if auto else float(lam_cfg),
                         shape_control="FWHM Coefficient", coeff=float(cfg_drt.get("rbf_shape", 0.5)))
        tau = np.asarray(getattr(out, "out_tau_vec"), float)
        gamma = np.asarray(getattr(out, "gamma"), float)
        x = np.asarray(getattr(out, "x"), float)
        R_inf = float(getattr(out, "R_inf", x[1] if len(x) > 1 else np.nan))
        L = float(getattr(out, "L_0", x[0] if len(x) > 0 else 0.0))
        Zf = np.asarray(getattr(out, "mu_Z_re"), float) + 1j * np.asarray(getattr(out, "mu_Z_im"), float)
        lam = float(getattr(out, "lambda_value", np.nan))
        order = np.argsort(tau)
        tau, gamma = tau[order], gamma[order]
        R_pol = float(np.trapezoid(gamma, np.log(tau)))
        resid = float(np.sqrt(np.mean(np.abs(d.Z - Zf) ** 2 / np.abs(d.Z) ** 2)) * 100)
        res = DRTResult(tau, gamma, tau, x, R_inf, L, lam, auto, R_pol, d.freq, d.Z, Zf, resid,
                        backend="pydrttools", label=label)
        res.peaks = detect_peaks(res, cfg_drt)
        return res
    except Exception as e:
        warnings.warn(f"pyDRTtools backend failed ({e}); using builtin solver")
        return None
