"""EIS -> DRT pipeline: per-file analysis, comparison overlays, synthetic validation."""
from __future__ import annotations
import re
import warnings
from pathlib import Path
import numpy as np
import pandas as pd

from .io.loaders import load_eis, EISData
from .eis.preprocess import preprocess_eis, linkk_check
from .eis.drt import compute_drt, compute_drt_builtin
from .eis.synth import synth_spectrum, circuit_truth
from .eis import plots as P
from .report import Report


def _safe(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s))


def analyze_one(d: EISData, cfg: dict, label: str, out_dir: Path, make_plots: bool = True):
    cfg_e = cfg["eis"]
    d2, info = preprocess_eis(d, cfg_e)
    kk = linkk_check(d2, cfg_e.get("kk", {}))
    if kk.warn:
        warnings.warn(f"[{label}] {kk.message}")
    res = compute_drt(d2, cfg_e["drt"], label)
    files = {}
    if make_plots:
        files["fig"] = P.plot_nyquist_bode_drt([res], [label], out_dir / "eis" / f"{_safe(label)}_drt.png",
                                               title=f"{label}  R∞={res.R_inf:.3g} Ω, L={res.L:.2e} H, R_pol={res.R_pol:.3g} Ω, "
                                                     f"fit resid={res.resid_pct:.2f}%", kk_notes=[kk.message])
        files["gcv"] = P.plot_gcv(res, out_dir / "eis" / f"{_safe(label)}_gcv.png", label)
        files["kk"] = P.plot_kk_residuals(d2.freq, kk, out_dir / "eis" / f"{_safe(label)}_kk.png")
    (out_dir / "eis").mkdir(parents=True, exist_ok=True)
    res.to_frame().to_csv(out_dir / "eis" / f"{_safe(label)}_drt.csv", index=False)
    pk = res.peaks.copy(); pk.insert(0, "label", label)
    pk.to_csv(out_dir / "eis" / f"{_safe(label)}_peaks.csv", index=False)
    pd.DataFrame({"freq_Hz": res.freq, "Zre": res.Z.real, "Zim": res.Z.imag, "Zre_fit": res.Z_fit.real,
                  "Zim_fit": res.Z_fit.imag}).to_csv(out_dir / "eis" / f"{_safe(label)}_fit.csv", index=False)
    summary = {"label": label, "source": d.source, "n_points": info["n_out"], "n_inductive": info["n_inductive_points"],
               "kk_M": kk.M, "kk_mu": kk.mu, "kk_rms_pct": kk.rms_res_pct, "kk_max_pct": kk.max_res_pct, "kk_warn": kk.warn,
               "R_inf_ohm": res.R_inf, "L_H": res.L, "R_pol_ohm": res.R_pol, "lambda": res.lam, "lambda_auto": res.lam_auto,
               "fit_resid_pct": res.resid_pct, "n_peaks": len(res.peaks), "backend": res.backend}
    return res, kk, summary, files


def run_eis(entries: list[dict], cfg: dict, out_dir: Path, report: Report | None = None,
            compare_key: str | None = None, title: str = "EIS → DRT") -> pd.DataFrame:
    """entries: [{label, path, profile?, group?, order?}] ; writes CSV/PNG, returns summary table."""
    out_dir = Path(out_dir)
    results, summaries, peaks_all, kks = [], [], [], []
    for e in entries:
        try:
            d = load_eis(e["path"], cfg, e.get("profile") or None)
        except Exception as ex:
            warnings.warn(f"skip {e['path']}: {ex}")
            summaries.append({"label": e["label"], "source": str(e["path"]), "error": str(ex)})
            continue
        res, kk, summ, files = analyze_one(d, cfg, e["label"], out_dir)
        for k in ("group", "order"):
            if k in e:
                summ[k] = e[k]
        results.append((e, res, kk, files)); summaries.append(summ)
        pk = res.peaks.copy(); pk.insert(0, "label", e["label"]); peaks_all.append(pk)
    summ_df = pd.DataFrame(summaries)
    summ_df.to_csv(out_dir / "eis_summary.csv", index=False)
    if peaks_all:
        pd.concat(peaks_all, ignore_index=True).to_csv(out_dir / "eis_peaks.csv", index=False)
    # comparison overlays
    cmp_files = []
    if len(results) > 1:
        key = compare_key or cfg["eis"].get("compare", {}).get("group_by") or "all"
        groups: dict[str, list] = {}
        for e, res, kk, files in results:
            g = str(e.get("group", "all")) if key != "all" else "all"
            groups.setdefault(g, []).append((e, res))
        if len(groups) == 1 or all(len(v) == 1 for v in groups.values()):
            groups = {"all": [(e, r) for e, r, _, _ in results]}
        for g, items in groups.items():
            items = sorted(items, key=lambda t: t[0].get("order", 0))
            ordered = any("order" in e for e, _ in items)
            f = P.plot_nyquist_bode_drt([r for _, r in items], [e["label"] for e, _ in items],
                                        out_dir / "eis" / f"compare_{_safe(g)}.png", title=f"comparison: {g}", ordered=ordered)
            cmp_files.append((g, f))
    if report is not None:
        report.section(title)
        report.text("전처리 → lin-KK 검증 → DRT(γ vs τ) → 피크 검출. kk_warn=True 인 행은 KK 잔차가 커서 해석에 주의.")
        cols = [c for c in ["label", "n_points", "kk_M", "kk_rms_pct", "kk_max_pct", "kk_warn", "R_inf_ohm", "L_H", "R_pol_ohm",
                            "lambda", "fit_resid_pct", "n_peaks", "backend", "error"] if c in summ_df]
        report.table(summ_df[cols], "요약 (eis_summary.csv)")
        if peaks_all:
            report.table(pd.concat(peaks_all, ignore_index=True), "DRT 피크 (eis_peaks.csv)")
        for g, f in cmp_files:
            report.image(f, f"비교: {g}")
        for e, res, kk, files in results:
            report.image(files.get("fig"), e["label"])
            if kk.warn:
                report.html(f'<p class="warn">{kk.message}</p>')
            report.image(files.get("gcv"), f"{e['label']} — λ selection (GCV / L-curve)")
    return summ_df


# ----------------------------------------------------------------------------
# synthetic validation
# ----------------------------------------------------------------------------
def _match_peaks(truth: pd.DataFrame, peaks: pd.DataFrame, tol_dec: float = 0.5) -> pd.DataFrame:
    rows = []
    used = set()
    for _, t in truth.iterrows():
        best, bd = None, np.inf
        for i, p in peaks.iterrows():
            dd = abs(np.log10(p.tau_s) - np.log10(t.tau_true_s))
            if dd < bd and i not in used:
                best, bd = i, dd
        if best is not None and bd <= tol_dec:
            used.add(best); p = peaks.loc[best]
            rows.append({**t.to_dict(), "tau_found_s": p.tau_s, "R_found_ohm": p.R_ohm,
                         "tau_err_decades": float(np.log10(p.tau_s / t.tau_true_s)),
                         "R_err_pct": float((p.R_ohm / t.R_true_ohm - 1) * 100), "found": True})
        else:
            rows.append({**t.to_dict(), "tau_found_s": np.nan, "R_found_ohm": np.nan, "tau_err_decades": np.nan,
                         "R_err_pct": np.nan, "found": False})
    return pd.DataFrame(rows)


def validate_eis(cfg: dict, out_dir: Path, report: Report | None = None, write_files: Path | None = None) -> dict:
    """Synthetic circuits x noise levels -> recovery table + lambda sweep figure.
    If write_files is given, spectra are also written as .mpt/.DTA/.csv to exercise the loaders."""
    out_dir = Path(out_dir); (out_dir / "eis").mkdir(parents=True, exist_ok=True)
    cs = cfg["synthetic"]["eis"]
    rng = np.random.default_rng(int(cfg["synthetic"].get("seed", 0)))
    f_hi, f_lo = cs["freq_range"]
    rec_rows, summ_rows, figs = [], [], []
    meta_rows = []
    for ci, spec in enumerate(cs["circuits"]):
        truth = circuit_truth(spec)
        for noise in cs["noise_pct"]:
            label = f"{spec['name']}_noise{noise}"
            d = synth_spectrum(spec, f_hi, f_lo, int(cs.get("points_per_decade", 10)), float(noise), rng, label)
            if write_files is not None:
                fn = _write_synthetic_eis(d, write_files, label, ci)
                meta_rows.append({"cell_id": label, "file": fn, "type": "eis", "surface": spec["name"], "electrolyte": f"noise{noise}",
                                  "active_mass_mg": "", "area_cm2": "", "profile": "generic" if fn.endswith(".csv") else "",
                                  "note": "synthetic"})
                d = load_eis(write_files / fn, cfg, "generic" if fn.endswith(".csv") else None)   # round-trip through loader
            res, kk, summ, files = analyze_one(d, cfg, label, out_dir)
            summ.update({"circuit": spec["name"], "noise_pct": noise, "R0_true": spec.get("R0"), "L_true": spec.get("L", 0.0),
                         "R0_err_pct": (res.R_inf / spec["R0"] - 1) * 100 if spec.get("R0") else np.nan,
                         "R_pol_true": sum(r for r, _ in spec.get("RC", [])) + (spec["FLW"][0] if spec.get("FLW") else 0),
                         "R_pol_err_pct": (res.R_pol / (sum(r for r, _ in spec.get("RC", [])) + (spec["FLW"][0] if spec.get("FLW") else 0)) - 1) * 100})
            summ_rows.append(summ)
            m = _match_peaks(truth, res.peaks); m.insert(0, "noise_pct", noise); m.insert(0, "circuit", spec["name"])
            m["n_peaks_found_total"] = len(res.peaks); m["lambda"] = res.lam
            rec_rows.append(m); figs.append((label, files.get("fig")))
    rec = pd.concat(rec_rows, ignore_index=True)
    summ = pd.DataFrame(summ_rows)
    rec.to_csv(out_dir / "eis_validation_peaks.csv", index=False)
    summ.to_csv(out_dir / "eis_validation_summary.csv", index=False)
    # lambda sweep on the "close peaks" circuit
    spec = next((s for s in cs["circuits"] if "close" in s["name"]), cs["circuits"][-1])
    d = synth_spectrum(spec, f_hi, f_lo, int(cs.get("points_per_decade", 10)), float(cs["noise_pct"][-1]), rng, "lambda_sweep")
    d, _ = preprocess_eis(d, cfg["eis"])
    base = compute_drt_builtin(d, cfg["eis"]["drt"])
    lams = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]
    sweep = [compute_drt_builtin(d, cfg["eis"]["drt"], lam=l) for l in lams]
    if not any(np.isclose(base.lam, l, rtol=0.05) for l in lams):
        sweep.append(base); sweep.sort(key=lambda r: r.lam)
    f_sweep = P.plot_lambda_sweep(sweep, base.lam, circuit_truth(spec), out_dir / "eis" / "lambda_sweep.png",
                                  title=f"λ sweep — {spec['name']} (two RC with τ ratio {spec['RC'][1][1]/spec['RC'][0][1]:.0f}, noise {cs['noise_pct'][-1]}%)")
    sweep_tbl = pd.DataFrame([{"lambda": r.lam, "gcv_choice": bool(np.isclose(r.lam, base.lam, rtol=0.05)), "n_peaks": len(r.peaks),
                               "peak_taus": ", ".join(f"{t:.2e}" for t in r.peaks.tau_s), "peak_R": ", ".join(f"{t:.2f}" for t in r.peaks.R_ohm),
                               "fit_resid_pct": r.resid_pct} for r in sweep])
    sweep_tbl.to_csv(out_dir / "eis_lambda_sweep.csv", index=False)
    if report is not None:
        report.section("EIS→DRT 검증 (합성 스펙트럼)")
        report.text("등가회로(R0 + RC 2–3개 (+ 유한길이 Warburg) + L)로 합성 → 0.5 %, 1 % 노이즈 → 로더 왕복(.mpt/.DTA/.csv) → DRT. "
                    "FLW(유한길이 Warburg)는 DRT에서 τ_k = 4τ_W/((2k−1)²π²), R_k = 8R_W/((2k−1)²π²) 의 급수로 나타나므로 k=1..3 항을 '참값'으로 둠.")
        cols = ["circuit", "noise_pct", "kk_rms_pct", "kk_warn", "lambda", "R_inf_ohm", "R0_true", "R0_err_pct", "L_H", "L_true",
                "R_pol_ohm", "R_pol_true", "R_pol_err_pct", "fit_resid_pct", "n_peaks"]
        report.table(summ[[c for c in cols if c in summ]], "스펙트럼별 요약 (eis_validation_summary.csv)")
        report.table(rec[["circuit", "noise_pct", "element", "tau_true_s", "tau_found_s", "tau_err_decades", "R_true_ohm", "R_found_ohm", "R_err_pct", "found", "n_peaks_found_total"]],
                     "피크 복원 (eis_validation_peaks.csv)")
        report.image(f_sweep, "λ에 따른 피크 분리/합침 (점선 = 참 τ)")
        report.table(sweep_tbl, "λ sweep 표")
        for label, f in figs:
            report.image(f, label)
    return {"recovery": rec, "summary": summ, "sweep": sweep_tbl, "meta_rows": meta_rows}


def _write_synthetic_eis(d: EISData, root: Path, label: str, idx: int) -> str:
    (root / "eis").mkdir(parents=True, exist_ok=True)
    fmt = ["mpt", "dta", "csv"][idx % 3]
    if fmt == "mpt":
        fn = f"eis/{label}.mpt"
        with open(root / fn, "w", encoding="latin-1") as f:
            f.write("EC-Lab ASCII FILE\nNb header lines : 5\n\nPotentio Electrochemical Impedance Spectroscopy (synthetic)\n")
            f.write("freq/Hz\tRe(Z)/Ohm\t-Im(Z)/Ohm\t|Z|/Ohm\tPhase(Z)/deg\ttime/s\t<Ewe>/V\tcycle number\n")
            for i, (fr, z) in enumerate(zip(d.freq, d.Z)):
                f.write(f"{fr:.6E}\t{z.real:.6E}\t{-z.imag:.6E}\t{abs(z):.6E}\t{np.degrees(np.angle(z)):.4f}\t{i*2.0:.3f}\t3.7\t1\n".replace(".", ","))
    elif fmt == "dta":
        fn = f"eis/{label}.DTA"
        with open(root / fn, "w") as f:
            f.write("EXPLAIN\nTAG\tEISPOT\nTITLE\tLABEL\tsynthetic EIS\n\n")
            f.write(f"ZCURVE\tTABLE\t{len(d.freq)}\n\tPt\tTime\tFreq\tZreal\tZimag\tZsig\tZmod\tZphz\tIdc\tVdc\tIERange\n")
            f.write("\t#\ts\tHz\tohm\tohm\tV\tohm\t°\tA\tV\t#\n")
            for i, (fr, z) in enumerate(zip(d.freq, d.Z)):
                f.write(f"\t{i}\t{i*2.0:.2f}\t{fr:.6E}\t{z.real:.6E}\t{z.imag:.6E}\t0.01\t{abs(z):.6E}\t{np.degrees(np.angle(z)):.4f}\t0\t3.7\t8\n")
    else:
        fn = f"eis/{label}.csv"
        pd.DataFrame({"Frequency (Hz)": d.freq, "Z' (ohm)": d.Z.real, "Z'' (ohm)": d.Z.imag}).to_csv(root / fn, index=False)
    return fn
