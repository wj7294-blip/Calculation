"""Synthetic EIS spectra from equivalent circuits (for validation)."""
from __future__ import annotations
import numpy as np
import pandas as pd
from ..io.loaders import EISData


def circuit_Z(freq: np.ndarray, spec: dict) -> np.ndarray:
    w = 2 * np.pi * np.asarray(freq, float)
    Z = np.full_like(w, float(spec.get("R0", 0.0)), dtype=complex)
    Z += 1j * w * float(spec.get("L", 0.0))
    for R, tau in spec.get("RC", []):
        Z += R / (1 + 1j * w * tau)
    if spec.get("FLW"):
        RW, tW = spec["FLW"]
        s = np.sqrt(1j * w * tW)
        Z += RW * np.tanh(s) / s
    if spec.get("W"):   # semi-infinite Warburg sigma
        sig = float(spec["W"])
        Z += sig / np.sqrt(w) * (1 - 1j)
    return Z


def circuit_truth(spec: dict, n_flw_terms: int = 3) -> pd.DataFrame:
    """Expected DRT peaks. RC element -> delta at tau with area R.
    Finite-length Warburg -> series: R_k = 8 R_W / ((2k-1)² π²), tau_k = 4 tau_W / ((2k-1)² π²)."""
    rows = []
    for i, (R, tau) in enumerate(spec.get("RC", [])):
        rows.append({"element": f"RC{i+1}", "tau_true_s": tau, "R_true_ohm": R})
    if spec.get("FLW"):
        RW, tW = spec["FLW"]
        for k in range(1, n_flw_terms + 1):
            c = (2 * k - 1) ** 2 * np.pi ** 2
            rows.append({"element": f"FLW_k{k}", "tau_true_s": 4 * tW / c, "R_true_ohm": 8 * RW / c})
    df = pd.DataFrame(rows)
    return df.sort_values("tau_true_s").reset_index(drop=True)


def synth_spectrum(spec: dict, f_hi: float = 1e5, f_lo: float = 1e-2, ppd: int = 10,
                   noise_pct: float = 0.5, rng: np.random.Generator | None = None,
                   label: str = "") -> EISData:
    rng = rng or np.random.default_rng(0)
    n = int(np.round((np.log10(f_hi) - np.log10(f_lo)) * ppd)) + 1
    freq = np.logspace(np.log10(f_hi), np.log10(f_lo), n)
    Z = circuit_Z(freq, spec)
    sig = noise_pct / 100.0 * np.abs(Z)
    Zn = Z + sig * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    return EISData(freq, Zn, label or spec.get("name", "synthetic"),
                   meta={"circuit": spec, "noise_pct": noise_pct, "Z_clean": Z})
