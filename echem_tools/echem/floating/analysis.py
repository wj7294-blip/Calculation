"""Per-cell analysis of a high-voltage floating (CC -> CV hold) test."""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
import pandas as pd
from ..io.loaders import ChronoData


@dataclass
class FloatResult:
    cell_id: str
    cv: pd.DataFrame                 # time_s, t_cv_s (from CV start), voltage_V, current_mA
    metrics: dict
    notes: list = field(default_factory=list)


def find_cv_segment(df: pd.DataFrame, V_hold: float, V_tol: float, min_points: int = 20,
                    cv_step_labels: list[str] | None = None) -> tuple[int, int, str]:
    """Return (i_start, i_end_exclusive, how). Primary: voltage window |V-V_hold|<=tol after reaching
    V_hold while current > 0 (charging). Optional: step labels from cycler export."""
    V = df["voltage_V"].to_numpy(float)
    I = df["current_mA"].to_numpy(float)
    n = len(df)
    # 1) step label based (if given & present)
    if cv_step_labels and "step" in df.columns:
        lab = df["step"].astype(str).str.strip().str.lower()
        m = lab.isin([s.lower() for s in cv_step_labels]).to_numpy()
        if m.sum() >= min_points:
            idx = np.where(m)[0]
            # take the longest contiguous run
            runs = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
            best = max(runs, key=len)
            return int(best[0]), int(best[-1]) + 1, "step_label"
    # 2) voltage based
    reach = np.where(V >= V_hold - V_tol)[0]
    if len(reach) == 0:
        raise ValueError(f"voltage never reaches V_hold-{V_tol} = {V_hold - V_tol:.3f} V (max {V.max():.3f} V)")
    i0 = int(reach[0])
    inwin = np.abs(V - V_hold) <= V_tol
    # end = first index after i0 where voltage leaves the window for >3 consecutive points, or current <= 0 (discharge)
    i1 = n
    bad = 0
    for i in range(i0, n):
        if (not inwin[i]) or (I[i] <= 0):
            bad += 1
            if bad > 3:
                i1 = i - bad + 1
                break
        else:
            bad = 0
    if i1 - i0 < min_points:
        raise ValueError(f"CV segment too short ({i1 - i0} points)")
    return i0, i1, "voltage_window"


def analyze_cell(cd: ChronoData, cell_id: str, active_mass_mg: float, area_cm2: float, cfg_f: dict) -> FloatResult:
    df = cd.df.sort_values("time_s").reset_index(drop=True)
    notes = []
    V_hold = float(cfg_f.get("V_hold", 4.6)); V_tol = float(cfg_f.get("V_tol", 0.02))
    i0, i1, how = find_cv_segment(df, V_hold, V_tol, int(cfg_f.get("min_cv_points", 20)),
                                  cd.meta.get("cv_step_labels"))
    cv = df.iloc[i0:i1].copy()
    cv["t_cv_s"] = cv["time_s"] - cv["time_s"].iloc[0]
    t = cv["t_cv_s"].to_numpy(float)
    I = cv["current_mA"].to_numpy(float)
    dur_h = t[-1] / 3600.0
    mass_g = float(active_mass_mg) / 1000.0 if active_mass_mg and np.isfinite(active_mass_mg) else np.nan
    area = float(area_cm2) if area_cm2 and np.isfinite(area_cm2) else np.nan
    if not np.isfinite(mass_g):
        notes.append("active_mass_mg missing -> per-gram metrics NaN")
    if not np.isfinite(area):
        notes.append("area_cm2 missing -> per-area metrics NaN")
    hold_h = float(cfg_f.get("hold_hours", 24))
    if dur_h < 0.9 * hold_h:
        notes.append(f"CV segment lasted {dur_h:.1f} h (< {hold_h} h expected)")

    # leakage current = mean over last window
    win = float(cfg_f.get("leak_window_h", 1.0)) * 3600
    m_last = t >= t[-1] - win
    I_leak_mA = float(np.mean(I[m_last]))
    I_leak_sd = float(np.std(I[m_last], ddof=1)) if m_last.sum() > 1 else np.nan

    # cumulative leakage capacity (mAh) = ∫ I dt / 3600
    Q_total_mAh = float(np.trapezoid(I, t) / 3600.0)
    tr = float(cfg_f.get("transient_h", 1.0)) * 3600
    m_tr = t >= tr
    Q_ex_mAh = float(np.trapezoid(I[m_tr], t[m_tr]) / 3600.0) if m_tr.sum() > 2 else np.nan
    Q_tr_mAh = Q_total_mAh - Q_ex_mAh if np.isfinite(Q_ex_mAh) else np.nan

    # log-log decay slope: log10(I) = a + n log10(t) in fit window
    w0, w1 = cfg_f.get("decay_fit_window_h", [2.0, 24.0])
    mf = (t >= float(w0) * 3600) & (t <= float(w1) * 3600) & (I > 0)
    if mf.sum() >= 5:
        p = np.polyfit(np.log10(t[mf]), np.log10(I[mf]), 1)
        slope, intercept = float(p[0]), float(p[1])
        pred = np.polyval(p, np.log10(t[mf]))
        r2 = 1 - np.sum((np.log10(I[mf]) - pred) ** 2) / np.sum((np.log10(I[mf]) - np.mean(np.log10(I[mf]))) ** 2)
    else:
        slope = intercept = r2 = np.nan
        notes.append("too few points in decay fit window")
    # early transient slope (first hour) for reference
    me = (t > 60) & (t <= tr) & (I > 0)
    slope_early = float(np.polyfit(np.log10(t[me]), np.log10(I[me]), 1)[0]) if me.sum() >= 5 else np.nan

    metrics = {
        "cell_id": cell_id,
        "cv_detect": how,
        "cv_start_s": float(df["time_s"].iloc[i0]), "cv_duration_h": dur_h, "n_cv_points": int(len(cv)),
        "V_cv_mean": float(cv["voltage_V"].mean()),
        "I_cc_mA": float(df["current_mA"].iloc[:i0].mean()) if i0 > 5 else np.nan,
        "leak_mA": I_leak_mA, "leak_mA_sd": I_leak_sd,
        "leak_mA_g": I_leak_mA / mass_g, "leak_uA_cm2": I_leak_mA * 1000.0 / area,
        "Q_leak_mAh": Q_total_mAh, "Q_leak_mAh_g": Q_total_mAh / mass_g,
        "Q_leak_ex_transient_mAh": Q_ex_mAh, "Q_leak_ex_transient_mAh_g": Q_ex_mAh / mass_g,
        "Q_transient_mAh_g": Q_tr_mAh / mass_g,
        "decay_slope": slope, "decay_intercept": intercept, "decay_r2": float(r2) if np.isfinite(r2) else np.nan,
        "decay_slope_early": slope_early,
        "I_1h_mA_g": float(np.interp(3600, t, I)) / mass_g if t[-1] >= 3600 else np.nan,
        "I_end_mA_g": float(I[-1]) / mass_g,
        "active_mass_mg": active_mass_mg, "area_cm2": area_cm2,
        "notes": "; ".join(notes),
    }
    return FloatResult(cell_id, cv[["time_s", "t_cv_s", "voltage_V", "current_mA"]], metrics, notes)
