"""Synthetic floating-test data: surface x electrolyte x n cells, written as .mpt / .DTA / Neware CSV
so that all three loaders are exercised. Returns the metadata DataFrame and the hidden truth table."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd


def _cell_curve(rng, leak_mA, mass_g, cfg_f, cfg_s):
    """CC at C-rate to V_hold, then CV for hold_hours. Returns DataFrame(time_s, voltage_V, current_mA, step)."""
    dt = float(cfg_s.get("sample_interval_s", 30))
    V_hold = float(cfg_f.get("V_hold", 4.6))
    Qn = float(cfg_f.get("nominal_capacity_mAh_g", 200)) * mass_g      # mAh
    I_cc = float(cfg_f.get("c_rate", 0.1)) * Qn                          # mA
    # CC: from 3.6 V, simple OCV-ish curve reaching V_hold at ~75 % of 1/C-rate hours
    t_cc = 0.75 / float(cfg_f.get("c_rate", 0.1)) * 3600 * rng.uniform(0.9, 1.1)
    n_cc = int(t_cc / dt)
    tc = np.arange(n_cc) * dt
    x = tc / t_cc
    V_cc = 3.6 + (V_hold - 3.6) * (0.15 * x + 0.85 * x ** 3) + 0.01 * rng.standard_normal(n_cc) * 0.2
    V_cc[-1] = V_hold
    I_cc_arr = np.full(n_cc, I_cc) * (1 + 0.002 * rng.standard_normal(n_cc))
    # CV: I(t) = I_leak + A t^-n, A set so that current at t=dt is ~ I_cc (continuity)
    hold_s = float(cfg_f.get("hold_hours", 24)) * 3600
    n_cv = int(hold_s / dt) + 1
    tv = np.arange(1, n_cv + 1) * dt
    nexp = float(cfg_s.get("decay_exponent", 0.55)) * rng.uniform(0.9, 1.1)
    A = (I_cc - leak_mA) * dt ** nexp
    # slow additional relaxation term (electrolyte oxidation build-up) for realism
    I_model = leak_mA * (1 + 0.15 * np.exp(-tv / (4 * 3600))) + A * tv ** (-nexp)
    I_cv = I_model * (1 + float(cfg_s.get("noise_rel", 0.02)) * rng.standard_normal(n_cv))
    I_cv += 0.002 * I_cc * rng.standard_normal(n_cv)
    V_cv = V_hold + 0.002 * rng.standard_normal(n_cv)
    df = pd.DataFrame({
        "time_s": np.concatenate([tc, tc[-1] + tv]),
        "voltage_V": np.concatenate([V_cc, V_cv]),
        "current_mA": np.concatenate([I_cc_arr, I_cv]),
        "step": ["CC_Chg"] * n_cc + ["CV_Chg"] * n_cv,
    })
    last = tv >= tv[-1] - 3600
    df.attrs["I_model_last1h_mA"] = float(I_model[last].mean())      # noise-free value of what leak metric measures
    return df


def write_mpt(df: pd.DataFrame, path: Path):
    lines = ["EC-Lab ASCII FILE", "Nb header lines : 6", "", "Chronopotentiometry / Chronoamperometry (synthetic)", "",
             "mode\tox/red\terror\ttime/s\tcontrol/V/mA\tEwe/V\t<I>/mA\tcycle number\tNs"]
    mode = np.where(df["step"] == "CC_Chg", 1, 2)
    ns = np.where(df["step"] == "CC_Chg", 0, 1)
    with open(path, "w", encoding="latin-1") as f:
        f.write("\n".join(lines) + "\n")
        for i, r in enumerate(df.itertuples(index=False)):
            # BioLogic style: decimal comma to exercise the parser
            f.write(f"{mode[i]}\t1\t0\t{r.time_s:.3f}\t0\t{r.voltage_V:.6f}\t{r.current_mA:.6f}\t0\t{ns[i]}\n".replace(".", ","))


def write_dta(df: pd.DataFrame, path: Path):
    hdr = ["EXPLAIN", "TAG\tCHRONOA", "TITLE\tLABEL\tSynthetic floating test\tTest &Identifier", "NOTES\t0", ""]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(hdr) + "\n")
        for k, (name, g) in enumerate(df.groupby("step", sort=False)):
            f.write(f"CURVE{k+1}\tTABLE\t{len(g)}\n")
            f.write("\tPt\tT\tVf\tIm\tVu\tSig\tAch\tIERange\tOver\n")
            f.write("\t#\ts\tV vs. Ref.\tA\tV\tV\tV\t#\tbits\n")
            for j, r in enumerate(g.itertuples(index=False)):
                f.write(f"\t{j}\t{r.time_s:.4f}\t{r.voltage_V:.6E}\t{r.current_mA/1000:.6E}\t0\t0\t0\t10\t..........\n")


def write_neware_csv(df: pd.DataFrame, path: Path):
    out = pd.DataFrame({
        "Record Index": np.arange(1, len(df) + 1),
        "Step Type": df["step"].values,
        "Cycle Index": 1,
        "Total Time(s)": df["time_s"].round(2).values,
        "Voltage(V)": df["voltage_V"].round(5).values,
        "Current(A)": (df["current_mA"] / 1000).round(9).values,       # exercise A->mA scale
        "Capacity(mAh)": 0.0,
    })
    out.to_csv(path, index=False)


def make_synthetic_floating(cfg: dict, out_dir: Path, rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame]:
    cfg_s = cfg["synthetic"]["floating"]; cfg_f = cfg["floating"]
    out_dir = Path(out_dir); (out_dir / "float").mkdir(parents=True, exist_ok=True)
    fmts = cfg_s.get("formats", ["mpt", "dta", "neware_csv"])
    rows, truth = [], []
    k = 0
    for s in cfg_s["surfaces"]:
        for e in cfg_s["electrolytes"]:
            mult = cfg_s["surface_effect"][s] * cfg_s["electrolyte_effect"][e] * cfg_s.get("interaction", {}).get(f"{s}|{e}", 1.0)
            mu = cfg_s["base_leak_mA_g"] * mult
            for i in range(int(cfg_s["n_per_group"])):
                k += 1
                mass_mg = float(cfg_s["active_mass_mg"]) * rng.uniform(0.95, 1.05)
                area = float(cfg_s["area_cm2"])
                leak_g = mu * (1 + float(cfg_s["cell_cv"]) * rng.standard_normal())
                leak_g = max(leak_g, 0.05 * mu)
                cell_id = f"{s}_{e}_{i+1:02d}"
                fmt = fmts[(k - 1) % len(fmts)]
                df = _cell_curve(rng, leak_g * mass_mg / 1000, mass_mg / 1000, cfg_f, cfg_s)
                if fmt == "mpt":
                    fn = f"float/{cell_id}.mpt"; write_mpt(df, out_dir / fn); prof = ""
                elif fmt == "dta":
                    fn = f"float/{cell_id}.DTA"; write_dta(df, out_dir / fn); prof = ""
                else:
                    fn = f"float/{cell_id}.csv"; write_neware_csv(df, out_dir / fn); prof = "neware"
                rows.append({"cell_id": cell_id, "file": fn, "type": "float", "surface": s, "electrolyte": e,
                             "active_mass_mg": round(mass_mg, 3), "area_cm2": area, "profile": prof,
                             "note": f"synthetic ({fmt})"})
                truth.append({"cell_id": cell_id, "surface": s, "electrolyte": e, "true_leak_mA_g": leak_g,
                              "true_I_last1h_mA_g": df.attrs["I_model_last1h_mA"] / (mass_mg / 1000),
                              "group_mean_leak_mA_g": mu, "interaction_factor": cfg_s.get("interaction", {}).get(f"{s}|{e}", 1.0)})
    return pd.DataFrame(rows), pd.DataFrame(truth)
