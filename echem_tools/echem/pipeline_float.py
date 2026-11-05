"""Floating leakage-current pipeline: per-cell metrics, group statistics, synthetic validation."""
from __future__ import annotations
import re
import warnings
from pathlib import Path
import numpy as np
import pandas as pd

from .io.loaders import load_chrono
from .floating.analysis import analyze_cell
from .floating.stats import group_stats
from .floating.synth import make_synthetic_floating
from .floating import plots as P
from .report import Report

METRIC_LABELS = {
    "leak_mA_g": "leakage current (mA/g)", "leak_uA_cm2": "leakage current (µA/cm²)",
    "Q_leak_mAh_g": "cumulative leak capacity, full CV (mAh/g)",
    "Q_leak_ex_transient_mAh_g": "cumulative leak capacity, transient excluded (mAh/g)",
    "decay_slope": "log-log decay slope (−n)",
}


def _safe(s):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s))


def run_floating(meta: pd.DataFrame, cfg: dict, out_dir: Path, report: Report | None = None,
                 per_cell_plots: bool = True, title: str = "고전압 플로팅 누설전류") -> dict:
    out_dir = Path(out_dir); (out_dir / "float").mkdir(parents=True, exist_ok=True)
    cfg_f = cfg["floating"]
    results, rows = [], []
    for _, m in meta.iterrows():
        try:
            cd = load_chrono(m["path"], cfg, (m.get("profile") or None) if isinstance(m.get("profile"), str) and m.get("profile") else None)
            r = analyze_cell(cd, m["cell_id"], m.get("active_mass_mg", np.nan), m.get("area_cm2", np.nan), cfg_f)
        except Exception as ex:
            warnings.warn(f"skip {m['cell_id']} ({m['path']}): {ex}")
            rows.append({"cell_id": m["cell_id"], "error": str(ex)})
            continue
        results.append(r)
        row = {**r.metrics, "file": m["file"], "format": cd.meta.get("format")}
        for c in ("surface", "electrolyte", "note"):
            if c in m:
                row[c] = m[c]
        rows.append(row)
        r.cv.to_csv(out_dir / "float" / f"{_safe(r.cell_id)}_cv.csv", index=False)
        if per_cell_plots:
            P.plot_single_cell(r, out_dir / "float" / f"{_safe(r.cell_id)}_cv.png", cfg_f)
    cells = pd.DataFrame(rows)
    fA, fB = cfg_f["stats"].get("factors", ["surface", "electrolyte"])
    stats_out = None
    ok = cells[cells.get("error", pd.Series(index=cells.index, dtype=object)).isna()] if "error" in cells else cells
    if len(ok) and fA in ok and fB in ok and ok[fA].notna().any():
        for c in (fA, fB):
            ok[c] = pd.Categorical(ok[c], categories=list(pd.unique(meta[c].dropna())))
        stats_out = group_stats(ok, cfg_f["stats"])
        for c in stats_out["outliers"].columns:
            cells.loc[stats_out["outliers"].index, c] = stats_out["outliers"][c]
        stats_out["summary"].to_csv(out_dir / "float_group_summary.csv", index=False)
    cells.to_csv(out_dir / "float_cells.csv", index=False)

    figs = {}
    if results:
        figs["curves"] = P.plot_cell_curves(results, meta, out_dir / "float" / "all_cells_decay.png",
                                            color_by=fB if fB in meta else "cell_id")
    if stats_out is not None:
        for mname in stats_out["metrics"]:
            r = stats_out["metrics"][mname]
            pv = {}
            if "anova" in r:
                a = r["anova"]
                pv = {fA: a.loc[fA, "PR(>F)"], fB: a.loc[fB, "PR(>F)"], f"{fA}×{fB}": a.loc[f"{fA}:{fB}", "PR(>F)"]}
            figs[f"bars_{mname}"] = P.plot_group_bars(ok, mname, fA, fB, stats_out["outliers"][f"outlier_{mname}"],
                                                      out_dir / "float" / f"group_{mname}.png", METRIC_LABELS.get(mname, mname))
            figs[f"inter_{mname}"] = P.plot_interaction(ok, mname, fA, fB, out_dir / "float" / f"interaction_{mname}.png",
                                                        METRIC_LABELS.get(mname, mname), pv)
        _write_stats_files(stats_out, out_dir)
    if report is not None:
        _report_floating(report, title, cells, stats_out, figs, cfg_f, results, out_dir)
    return {"cells": cells, "stats": stats_out, "results": results, "figs": figs}


def _write_stats_files(st: dict, out_dir: Path):
    fA, fB = st["factors"]
    an, tk, npr = [], [], []
    for m, r in st["metrics"].items():
        if "anova" in r:
            a = r["anova"].copy(); a.insert(0, "metric", m); a.insert(1, "term", a.index)
            a["shapiro_p"] = r["shapiro_p"]; a["levene_p"] = r["levene_p"]; an.append(a)
            for key in ("tukey_all", "tukey_within_B", "tukey_A", "tukey_B"):
                t = r.get(key)
                if t is not None and len(t):
                    t = t.copy(); t.insert(0, "metric", m); t.insert(1, "comparison", key); tk.append(t)
            if r.get("nonparametric"):
                s = r["nonparametric"]["srh"].copy(); s.insert(0, "metric", m); s.insert(1, "term", s.index)
                s["kruskal_p_" + fA] = r["nonparametric"]["kruskal_p"].get(fA)
                s["kruskal_p_" + fB] = r["nonparametric"]["kruskal_p"].get(fB)
                npr.append(s)
    if an:
        pd.concat(an, ignore_index=True).to_csv(out_dir / "float_anova.csv", index=False)
    if tk:
        pd.concat(tk, ignore_index=True).to_csv(out_dir / "float_posthoc_tukey.csv", index=False)
    if npr:
        pd.concat(npr, ignore_index=True).to_csv(out_dir / "float_nonparametric.csv", index=False)


def _report_floating(rep: Report, title, cells, st, figs, cfg_f, results, out_dir):
    rep.section(title)
    rep.text(f"프로토콜: CC(C/{1/cfg_f.get('c_rate',0.1):.0f}) → {cfg_f['V_hold']} V 정전압 {cfg_f['hold_hours']} h. "
             f"누설전류 = 정전압 구간 마지막 {cfg_f['leak_window_h']} h 평균; 누적 누설 용량 = ∫I dt (전체 / 초기 {cfg_f['transient_h']} h 제외); "
             f"감쇠 기울기 = log10(I) vs log10(t) 기울기 ({cfg_f['decay_fit_window_h'][0]}–{cfg_f['decay_fit_window_h'][1]} h).")
    show = [c for c in ["cell_id", "surface", "electrolyte", "format", "cv_detect", "cv_duration_h", "leak_mA_g", "leak_uA_cm2",
                        "Q_leak_mAh_g", "Q_leak_ex_transient_mAh_g", "decay_slope", "decay_r2", "notes", "error"] if c in cells]
    oc = [c for c in cells.columns if c.startswith("outlier_")]
    if oc:
        cells = cells.copy(); cells["outlier_flags"] = cells[oc].apply(lambda r: ",".join(c[8:] for c in oc if r[c] is True), axis=1)
        show.append("outlier_flags")
    rep.table(cells[show], "셀별 결과 (float_cells.csv)")
    rep.image(figs.get("curves"), "전체 셀 전류 감쇠")
    if st is None:
        rep.text("그룹 통계 생략: metadata 에 surface/electrolyte 가 없거나 분석 가능한 셀이 없음.")
        return
    fA, fB = st["factors"]
    rep.table(st["summary"], "그룹별 평균 / SD / n (float_group_summary.csv)")
    for m, r in st["metrics"].items():
        rep.section(f"지표: {METRIC_LABELS.get(m, m)}")
        if "anova" not in r:
            rep.text(f"ANOVA 불가: {r.get('error')}"); continue
        a = r["anova"]
        p_int = a.loc[f"{fA}:{fB}", "PR(>F)"]
        verdict = (f'<p>이원 분산분석 (type II): {fA} p={a.loc[fA,"PR(>F)"]:.3g}, {fB} p={a.loc[fB,"PR(>F)"]:.3g}, '
                   f'상호작용 {fA}×{fB} p=<b>{p_int:.3g}</b> → '
                   + (f'<span class="warn">SAM(표면) 효과가 전해질에 따라 다름 (상호작용 유의)</span>' if p_int < st["alpha"]
                      else '<span class="ok">상호작용 유의하지 않음 (표면 효과가 전해질 간 비슷)</span>')
                   + f'. R²={r["r2"]:.3f}. 가정검정: Shapiro p={r["shapiro_p"]:.3g}, Levene p={r["levene_p"]:.3g} '
                   + ('(정상)' if r["assumptions_ok"] else '<span class="warn">(정규성/등분산 의심 → 아래 비모수 결과 병행 참고)</span>') + '</p>')
        rep.html(verdict)
        rep.image(figs.get(f"inter_{m}"), "상호작용 그림")
        rep.image(figs.get(f"bars_{m}"), "그룹 비교")
        rep.table(a, "ANOVA 표")
        if r.get("tukey_within_B") is not None and len(r["tukey_within_B"]):
            rep.table(r["tukey_within_B"], f"사후검정 (Tukey HSD): 각 {fB} 내 {fA} 쌍 비교")
        if r.get("tukey_A") is not None:
            rep.table(r["tukey_A"], f"사후검정 (Tukey HSD): {fA} 주효과")
            rep.table(r["tukey_B"], f"사후검정 (Tukey HSD): {fB} 주효과")
        npr = r.get("nonparametric")
        if npr:
            rep.table(npr["srh"], "비모수: Scheirer–Ray–Hare (순위 기반 이원 분석)")
            rep.text("Kruskal–Wallis p: " + ", ".join(f"{k}={v:.3g}" for k, v in npr["kruskal_p"].items()))
            if "dunn_A" in npr:
                rep.table(npr["dunn_A"], f"Dunn 사후검정 (Holm 보정), {fA} 주효과 p-값 행렬")
    rep.text("전체 조합 쌍 비교(Tukey, 모든 그룹 쌍)는 float_posthoc_tukey.csv 참조.")


def validate_floating(cfg: dict, data_root: Path, out_dir: Path, report: Report | None = None) -> dict:
    rng = np.random.default_rng(int(cfg["synthetic"].get("seed", 0)) + 1)
    meta, truth = make_synthetic_floating(cfg, data_root, rng)
    meta["path"] = [str(Path(data_root) / f) for f in meta["file"]]
    out = run_floating(meta, cfg, out_dir, report, per_cell_plots=False, title="플로팅 검증 (합성 데이터: 5 표면 × 3 전해질 × n=4)")
    cells = out["cells"].merge(truth, on=["cell_id", "surface", "electrolyte"], how="left")
    cells["leak_recovery_err_pct"] = (cells["leak_mA_g"] / cells["true_I_last1h_mA_g"] - 1) * 100     # vs noise-free model, last 1 h
    cells["leak_vs_asymptote_pct"] = (cells["leak_mA_g"] / cells["true_leak_mA_g"] - 1) * 100         # vs planted steady-state leak
    cells.to_csv(out_dir / "float_validation_cells.csv", index=False)
    truth.to_csv(out_dir / "float_validation_truth.csv", index=False)
    cs = cfg["synthetic"]["floating"]
    if report is not None:
        report.section("플로팅 검증: 심어둔 효과 vs 검출")
        report.text(f"기본 누설 {cs['base_leak_mA_g']} mA/g × 표면 인자 {cs['surface_effect']} × 전해질 인자 {cs['electrolyte_effect']} × "
                    f"상호작용 {cs.get('interaction')} (곱셈), 셀 간 CV {cs['cell_cv']*100:.0f} %.")
        g = cells.groupby(["surface", "electrolyte"], observed=True).agg(planted=("group_mean_leak_mA_g", "first"),
                                                                        interaction_factor=("interaction_factor", "first"),
                                                                        true_asymptote_mean=("true_leak_mA_g", "mean"),
                                                                        true_last1h_mean=("true_I_last1h_mA_g", "mean"),
                                                                        est_mean=("leak_mA_g", "mean"), est_sd=("leak_mA_g", "std"),
                                                                        err_vs_model_pct=("leak_recovery_err_pct", "mean"),
                                                                        err_vs_asymptote_pct=("leak_vs_asymptote_pct", "mean")).reset_index()
        report.table(g, "그룹별 (mA/g): 심은 그룹 평균(planted) / 셀 참값(점근 누설) / 노이즈 없는 모델의 마지막 1 h 평균 / 추정치 / 오차")
        report.text(f"셀 단위 복원 오차 (노이즈 없는 모델의 마지막 1 h 평균 대비): 평균 {cells['leak_recovery_err_pct'].mean():.2f} %, "
                    f"절대평균 {cells['leak_recovery_err_pct'].abs().mean():.2f} %. "
                    f"점근(정상상태) 누설 대비는 평균 {cells['leak_vs_asymptote_pct'].mean():.1f} % 과대: 24 h 시점에도 t^-n 과도 전류 꼬리가 남아 있기 때문이며, "
                    f"이는 실제 플로팅 시험에서도 동일하게 생기는 체계적 편향이다(그룹 간 비교에는 영향이 작지만 절대값 해석 시 유의).")
    out["truth"] = truth; out["cells_vs_truth"] = cells; out["meta"] = meta
    return out
