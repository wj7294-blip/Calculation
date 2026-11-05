"""결과 분석: 그림(figures/)과 표(results/*.md) 생성.

  load_summary()          : results/summary.csv (+ results/database.csv 의 사용자 추가 분자) 읽기
  fig_homo_dipole_map()   : HOMO(및 IP) vs 쌍극자 지도  → figures/map_homo_vs_dipole.png
  fig_lumo_dipole_map()   : LUMO(및 EA) vs 쌍극자 지도  → figures/map_lumo_vs_dipole.png
  fig_zwitterion()        : 시리즈 C: N+–음이온 거리·쌍극자 vs m → figures/zw_distance_vs_m.png
  write_tables()          : README 용 표 (results/tables.md)
"""
import glob
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    import koreanize_matplotlib  # noqa: F401  (한글 폰트 등록)
    plt.rcParams["font.family"] = ["NanumGothic", "DejaVu Sans"]
except ImportError:
    pass
plt.rcParams["axes.unicode_minus"] = False

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

# 색: 말단기 (고정 순서, 검증된 팔레트 8색 + 기준용 회색)
TERMINAL_COLORS = {
    "H": "#52514e", "CH3": "#2a78d6", "F": "#eb6834", "CF3": "#e34948", "CN": "#6250d6",
    "OCH3": "#1baf7a", "COOH": "#eda100", "SO3H": "#008300", "NH2": "#e87ba4",
    "SO3-": "#008300", "COO-": "#eda100", "F(penta)": "#eb6834", "H(3,5-CF3)": "#e34948",
}
ANCHOR_MARKERS = {"SH": "o", "PO3H2": "s", "COOH": "^"}
SPACER_FILL = {"para-phenylene": True}


def load_summary():
    frames = []
    if os.path.exists("results/summary.csv"):
        frames.append(pd.read_csv("results/summary.csv"))
    if os.path.exists("results/database.csv"):
        frames.append(pd.read_csv("results/database.csv"))
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if len(df):
        df = df.drop_duplicates("name", keep="last")
    return df


def _scatter_map(df, ycol, ylabel, ax, title):
    d = df[(df.series != "val") & df[ycol].notna()]
    for _, r in d.iterrows():
        mk = ANCHOR_MARKERS.get(r.anchor, "D")
        if r.series == "C":
            mk = "D"
        col = TERMINAL_COLORS.get(str(r.terminal), "#8c6d4f")
        hollow = (r.spacer == "(CH2)3")
        ax.scatter(r.dipole, r[ycol], marker=mk, s=46, c="none" if hollow else col, edgecolors=col,
                   linewidths=1.4, zorder=3)
        ax.annotate(r["name"], (r.dipole, r[ycol]), fontsize=4.2, xytext=(2.5, 1.5),
                    textcoords="offset points", color="#333")
    ax.set_xlabel("쌍극자 모멘트 |μ| (D)")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.grid(True, lw=0.4, alpha=0.4)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def _legend(fig):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker=m, ls="", color="k", mfc="none", label=f"앵커 {a}") for a, m in ANCHOR_MARKERS.items()]
    h.append(Line2D([], [], marker="D", ls="", color="k", mfc="none", label="쯔비터이온(시리즈 C)"))
    h.append(Line2D([], [], marker="o", ls="", color="k", mfc="k", label="채움=페닐렌 / 빈=(CH2)3"))
    for t, c in TERMINAL_COLORS.items():
        if t in ("SO3-", "COO-", "F(penta)", "H(3,5-CF3)"):
            continue
        h.append(Line2D([], [], marker="s", ls="", color=c, label=f"말단 {t}"))
    fig.legend(handles=h, loc="lower center", ncol=5, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.0))


def fig_homo_dipole_map(df, path="figures/map_homo_vs_dipole.png"):
    fig, axes = plt.subplots(1, 2, figsize=(14, 7.5))
    _scatter_map(df, "HOMO", "HOMO (eV, GFN2-xTB)", axes[0], "HOMO vs 쌍극자")
    _scatter_map(df, "IP", "vertical IP (eV, IPEA-xTB)", axes[1], "IP vs 쌍극자")
    _legend(fig)
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    fig.savefig(path, dpi=200)
    plt.close(fig)


def fig_lumo_dipole_map(df, path="figures/map_lumo_vs_dipole.png"):
    fig, axes = plt.subplots(1, 2, figsize=(14, 7.5))
    _scatter_map(df, "LUMO", "LUMO (eV, GFN2-xTB)", axes[0], "LUMO vs 쌍극자")
    _scatter_map(df, "EA", "vertical EA (eV, IPEA-xTB)", axes[1], "EA vs 쌍극자")
    _legend(fig)
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    fig.savefig(path, dpi=200)
    plt.close(fig)


def zw_conformer_table(df):
    """시리즈 C 모든 conformer 의 N+–X 거리, 상대 에너지, 볼츠만 가중치."""
    rows = []
    for _, r in df[df.series == "C"].iterrows():
        p = f"results/molecules/{r['name']}.json"
        if not os.path.exists(p):
            continue
        d = json.load(open(p))
        for c in d["conformers"]:
            rows.append(dict(name=r["name"], n=int(r.n), m=int(r.m), terminal=r.terminal, conf=c["idx"],
                             rel_kcal=c["rel_kcal"], boltz_w=c["boltz_w"], zw_dist=c.get("zw_dist"),
                             dipole=c["dipole"], connectivity_changed=bool(c["connectivity_changed"])))
    t = pd.DataFrame(rows)
    t.to_csv("results/seriesC_conformers.csv", index=False)
    return t


def fig_zwitterion(df, path="figures/zw_distance_vs_m.png"):
    c = df[df.series == "C"].copy()
    if not len(c):
        return
    c["n"] = c.n.astype(int)
    c["m"] = c.m.astype(int)
    t = zw_conformer_table(df)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
    styles = {("SO3-", 2): ("#008300", "-"), ("SO3-", 3): ("#008300", "--"),
              ("COO-", 2): ("#eda100", "-"), ("COO-", 3): ("#eda100", "--")}
    for (term, n), (col, ls) in styles.items():
        s = c[(c.terminal == term) & (c.n == n)].sort_values("m")
        lab = f"{term}, n={n}"
        axes[0].errorbar(s.m, s.zw_dist_boltz, yerr=s.zw_dist_std, color=col, ls=ls, marker="o", ms=5,
                         capsize=2, lw=1.6, label=lab)
        axes[0].plot(s.m, s.zw_dist_best, color=col, ls="none", marker="x", ms=5, alpha=0.7)
        axes[2].plot(s.m, s.dipole, color=col, ls=ls, marker="o", ms=5, lw=1.6, label=lab)
        # 분포: conformer 전체 (점 크기 ∝ 볼츠만 가중치)
        tt = t[(t.terminal == term) & (t.n == n)]
        jitter = (0.12 if term == "SO3-" else -0.12) + (0.05 if n == 3 else -0.05)
        axes[1].scatter(tt.m + jitter, tt.zw_dist, s=4 + 60 * tt.boltz_w, color=col, alpha=0.5,
                        marker="o" if n == 2 else "^", linewidths=0, label=lab)
    axes[0].set_title("N⁺–X 거리 (볼츠만 가중 평균 ± 가중 표준편차, ×=최저 E 구조)", fontsize=9)
    axes[0].set_ylabel("거리 (Å), X = S(SO₃⁻) 또는 C(COO⁻)")
    axes[1].set_title("conformer 전체 거리 분포 (점 크기 ∝ 볼츠만 가중치)", fontsize=9)
    axes[2].set_title("최저 E 구조의 쌍극자 모멘트 (ALPB water)", fontsize=9)
    axes[2].set_ylabel("|μ| (D)")
    for ax in axes:
        ax.set_xlabel("m (N⁺와 음이온 사이 CH₂ 수)")
        ax.set_xticks([1, 2, 3, 4, 5])
        ax.grid(True, lw=0.4, alpha=0.4)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def fig_validation(df, path="figures/validation_ip.png"):
    v = df[df.series == "val"].dropna(subset=["exp_ip"])
    if not len(v):
        return
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    ax.scatter(v.exp_ip, v.IP, color="#2a78d6", label="IPEA-xTB vertical IP", zorder=3)
    ax.scatter(v.exp_ip, -v.HOMO, color="#eb6834", marker="s", label="−ε(HOMO), GFN2-xTB", zorder=3)
    for _, r in v.iterrows():
        ax.annotate(r["name"], (r.exp_ip, r.IP), fontsize=6, xytext=(3, 2), textcoords="offset points")
    lo, hi = 7, 13
    ax.plot([lo, hi], [lo, hi], color="#999", lw=0.8, ls="--")
    ax.set_xlabel("실험 IP (eV, NIST)")
    ax.set_ylabel("계산값 (eV)")
    ax.set_xlim(7.4, 9.6)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(True, lw=0.4, alpha=0.4)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def md_table(d, cols, fmt=None):
    fmt = fmt or {}
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in d.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                if np.isnan(v):
                    cells.append("")
                else:
                    cells.append(fmt.get(c, "{:.2f}").format(v))
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def write_tables(df, path="results/tables.md"):
    out = []
    main = df[df.series != "val"].copy()
    cols = ["name", "series", "solvent", "HOMO", "LUMO", "gap", "IP", "EA", "dipole", "dipole_axial",
            "q_anchor_atom", "q_terminal_group", "TPSA", "logP", "MW", "nF", "nRotB", "chain_bonds"]
    out.append("## 전체 결과 표 (최저 에너지 구조)\n")
    out.append(md_table(main.sort_values(["series", "name"]), cols, {"MW": "{:.1f}", "TPSA": "{:.1f}"}))
    # 검증
    v = df[df.series == "val"].copy()
    if len(v):
        v["neg_HOMO"] = -v.HOMO
        v = v.sort_values("exp_ip", ascending=False)
        out.append("\n\n## 검증: 실험 IP vs 계산 (기체상)\n")
        out.append(md_table(v, ["name", "exp_ip", "IP", "neg_HOMO", "HOMO", "gap"]))
        from scipy.stats import spearmanr, pearsonr
        if len(v) >= 3:
            rs = spearmanr(v.exp_ip, v.IP).correlation
            rp = pearsonr(v.exp_ip, v.IP)[0]
            rs2 = spearmanr(v.exp_ip, v.neg_HOMO).correlation
            mae = (v.IP - v.exp_ip).abs().mean()
            out.append(f"\n- IPEA-xTB IP: Spearman ρ = {rs:.2f}, Pearson r = {rp:.2f}, MAE = {mae:.2f} eV")
            out.append(f"\n- −ε(HOMO): Spearman ρ = {rs2:.2f}")
    # 결합 변화
    ch = main[main.n_conf_connectivity_changed > 0]
    out.append("\n\n## 최적화 후 결합 연결이 바뀐 conformer\n")
    if len(ch):
        out.append(md_table(ch, ["name", "n_conf_ok", "n_conf_connectivity_changed"]))
        out.append("\n(해당 conformer 는 최저 에너지 구조 선택과 볼츠만 평균에서 제외됨. 상세: results/connectivity_changes.csv)")
    else:
        out.append("없음")
    # 수렴 검사
    p50 = "results/convergence/summary_50.csv"
    if os.path.exists(p50):
        a = df.set_index("name")
        b = pd.read_csv(p50).set_index("name")
        rows = []
        for nm50, r in b.iterrows():
            nm = nm50.replace("_50conf", "")
            if nm in a.index:
                dE = (r.E_min_Eh - a.loc[nm, "E_min_Eh"]) * 627.5095
                rows.append(dict(name=nm, E20_Eh=a.loc[nm, "E_min_Eh"], E50_Eh=r.E_min_Eh, dE_kcal=dE,
                                 HOMO20=a.loc[nm, "HOMO"], HOMO50=r.HOMO, dip20=a.loc[nm, "dipole"], dip50=r.dipole,
                                 same=("예" if abs(dE) < 0.3 else "아니오")))
        cv = pd.DataFrame(rows)
        cv.to_csv("results/convergence/comparison.csv", index=False)
        out.append("\n\n## conformer 수 20 → 50 수렴 검사\n")
        out.append(md_table(cv, ["name", "E20_Eh", "E50_Eh", "dE_kcal", "HOMO20", "HOMO50", "dip20", "dip50", "same"],
                            {"E20_Eh": "{:.6f}", "E50_Eh": "{:.6f}"}))
        out.append("\n(dE = E_min(50) − E_min(20), kcal/mol. |dE| < 0.3 kcal/mol 이면 같은 최저 구조로 판정)")
    # 시리즈 C 요약
    c = df[df.series == "C"]
    if len(c):
        out.append("\n\n## 시리즈 C: N⁺–음이온 중심 거리\n")
        out.append(md_table(c.sort_values(["terminal", "n", "m"]),
                            ["name", "zw_dist_best", "zw_dist_boltz", "zw_dist_std", "zw_dist_min", "zw_dist_max", "dipole", "dipole_axial", "HOMO"]))
    open(path, "w").write("\n".join(out) + "\n")


def connectivity_report(df, path="results/connectivity_changes.csv"):
    rows = []
    for p in glob.glob("results/molecules/*.json"):
        d = json.load(open(p))
        for c in d["conformers"]:
            if c["connectivity_changed"]:
                rows.append(dict(name=d["summary"]["name"], conf=c["idx"], rel_kcal=c.get("rel_kcal"),
                                 changes=" ".join(c["connectivity_changed"])))
    pd.DataFrame(rows, columns=["name", "conf", "rel_kcal", "changes"]).to_csv(path, index=False)


def main():
    os.makedirs("figures", exist_ok=True)
    df = load_summary()
    fig_homo_dipole_map(df)
    fig_lumo_dipole_map(df)
    fig_zwitterion(df)
    fig_validation(df)
    connectivity_report(df)
    write_tables(df)
    print("figures/ 와 results/tables.md 생성 완료")


if __name__ == "__main__":
    main()
