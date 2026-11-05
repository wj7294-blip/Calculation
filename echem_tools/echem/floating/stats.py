"""Group statistics: two-way ANOVA (type II) with interaction, Tukey HSD post-hoc,
assumption checks (Shapiro-Wilk on residuals, Levene), nonparametric alternatives
(Scheirer-Ray-Hare for the two-way layout, Kruskal-Wallis per factor, Dunn post-hoc), outlier flags."""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
from scipy import stats


def flag_outliers(df: pd.DataFrame, metric: str, factors: list[str], method: str = "iqr") -> pd.Series:
    out = pd.Series(False, index=df.index)
    for _, g in df.groupby(factors, observed=True):
        x = g[metric].astype(float)
        if x.notna().sum() < 3:
            continue
        if method == "mad":
            med = x.median(); mad = (x - med).abs().median()
            if mad > 0:
                out.loc[g.index] = (0.6745 * (x - med).abs() / mad) > 3.5
        else:
            q1, q3 = x.quantile(0.25), x.quantile(0.75); iqr = q3 - q1
            out.loc[g.index] = (x < q1 - 1.5 * iqr) | (x > q3 + 1.5 * iqr)
    return out


def two_way_anova(df: pd.DataFrame, metric: str, fA: str, fB: str, alpha: float = 0.05) -> dict:
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    d = df[[metric, fA, fB]].dropna().copy()
    d = d.rename(columns={metric: "y"})
    res = {"metric": metric, "n": int(len(d))}
    if d[fA].nunique() < 2 or d[fB].nunique() < 2:
        res["error"] = "need >=2 levels per factor"
        return res
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = smf.ols(f"y ~ C({fA}) * C({fB})", data=d).fit()
        try:
            tab = sm.stats.anova_lm(model, typ=2)
        except Exception:
            tab = sm.stats.anova_lm(model, typ=1)
    tab = tab.rename(index=lambda s: s.replace(f"C({fA})", fA).replace(f"C({fB})", fB))
    ss_tot = tab["sum_sq"].sum()
    tab["partial_eta2"] = tab["sum_sq"] / (tab["sum_sq"] + tab.loc["Residual", "sum_sq"])
    tab.loc["Residual", "partial_eta2"] = np.nan
    res["anova"] = tab
    res["r2"] = float(model.rsquared)
    # assumptions
    resid = model.resid
    sw = stats.shapiro(resid) if 3 <= len(resid) <= 5000 else (np.nan, np.nan)
    groups = [g["y"].values for _, g in d.groupby([fA, fB], observed=True) if len(g) >= 2]
    lev = stats.levene(*groups, center="median") if len(groups) >= 2 else (np.nan, np.nan)
    res["shapiro_p"] = float(sw[1]); res["levene_p"] = float(lev[1])
    res["assumptions_ok"] = bool(res["shapiro_p"] > alpha and res["levene_p"] > alpha)
    # post-hoc Tukey: all cell combinations, and factor A within each level of B
    d["cell"] = d[fA].astype(str) + "|" + d[fB].astype(str)
    try:
        tk = pairwise_tukeyhsd(d["y"], d["cell"], alpha=alpha)
        res["tukey_all"] = pd.DataFrame(tk.summary().data[1:], columns=tk.summary().data[0])
    except Exception as e:
        res["tukey_all"] = pd.DataFrame({"error": [str(e)]})
    within = []
    for lvB, g in d.groupby(fB, observed=True):
        if g[fA].nunique() < 2:
            continue
        try:
            tk = pairwise_tukeyhsd(g["y"], g[fA], alpha=alpha)
            t = pd.DataFrame(tk.summary().data[1:], columns=tk.summary().data[0]); t.insert(0, fB, lvB)
            within.append(t)
        except Exception:
            pass
    res["tukey_within_B"] = pd.concat(within, ignore_index=True) if within else pd.DataFrame()
    # main-effect Tukey
    try:
        tk = pairwise_tukeyhsd(d["y"], d[fA], alpha=alpha)
        res["tukey_A"] = pd.DataFrame(tk.summary().data[1:], columns=tk.summary().data[0])
        tk = pairwise_tukeyhsd(d["y"], d[fB], alpha=alpha)
        res["tukey_B"] = pd.DataFrame(tk.summary().data[1:], columns=tk.summary().data[0])
    except Exception:
        pass
    return res


def scheirer_ray_hare(df: pd.DataFrame, metric: str, fA: str, fB: str) -> pd.DataFrame:
    """Scheirer-Ray-Hare: rank-transform, two-way ANOVA sums of squares, H = SS / MS_total ~ chi2."""
    d = df[[metric, fA, fB]].dropna().copy()
    d["r"] = stats.rankdata(d[metric])
    N = len(d)
    grand = d["r"].mean()
    ms_tot = d["r"].var(ddof=1)
    def ss(by):
        return sum(len(g) * (g["r"].mean() - grand) ** 2 for _, g in d.groupby(by, observed=True))
    ssA, ssB = ss(fA), ss(fB)
    ssAB_cells = ss([fA, fB]) - ssA - ssB
    a, b = d[fA].nunique(), d[fB].nunique()
    rows = []
    for name, s, dof in [(fA, ssA, a - 1), (fB, ssB, b - 1), (f"{fA}:{fB}", ssAB_cells, (a - 1) * (b - 1))]:
        H = s / ms_tot
        rows.append({"term": name, "H": H, "df": dof, "p": float(stats.chi2.sf(H, dof))})
    return pd.DataFrame(rows).set_index("term")


def nonparametric(df: pd.DataFrame, metric: str, fA: str, fB: str, alpha: float) -> dict:
    d = df[[metric, fA, fB]].dropna().copy()
    out = {"srh": scheirer_ray_hare(d, metric, fA, fB)}
    kw = {}
    for f in (fA, fB):
        groups = [g[metric].values for _, g in d.groupby(f, observed=True)]
        if len(groups) >= 2:
            kw[f] = float(stats.kruskal(*groups).pvalue)
    d["cell"] = d[fA].astype(str) + "|" + d[fB].astype(str)
    groups = [g[metric].values for _, g in d.groupby("cell", observed=True)]
    kw["cells"] = float(stats.kruskal(*groups).pvalue) if len(groups) >= 2 else np.nan
    out["kruskal_p"] = kw
    try:
        import scikit_posthocs as sp
        out["dunn_A"] = sp.posthoc_dunn(d, val_col=metric, group_col=fA, p_adjust="holm")
        within = {}
        for lvB, g in d.groupby(fB, observed=True):
            if g[fA].nunique() >= 2:
                within[lvB] = sp.posthoc_dunn(g, val_col=metric, group_col=fA, p_adjust="holm")
        out["dunn_A_within_B"] = within
    except Exception as e:
        out["dunn_error"] = str(e)
    return out


def group_stats(cells: pd.DataFrame, cfg_stats: dict) -> dict:
    """cells: one row per cell with factor columns and metric columns. Returns dict of results per metric."""
    fA, fB = cfg_stats.get("factors", ["surface", "electrolyte"])
    alpha = float(cfg_stats.get("alpha", 0.05))
    metrics = [m for m in cfg_stats.get("metrics", []) if m in cells.columns]
    mode = str(cfg_stats.get("nonparametric", "auto")).lower()
    out = {"factors": (fA, fB), "metrics": {}, "alpha": alpha}
    summ = (cells.groupby([fA, fB], observed=True)[metrics]
            .agg(["mean", "std", "count"]))
    summ.columns = [f"{m}_{s}" for m, s in summ.columns]
    out["summary"] = summ.reset_index()
    flags = pd.DataFrame(index=cells.index)
    for m in metrics:
        flags[f"outlier_{m}"] = flag_outliers(cells, m, [fA, fB], cfg_stats.get("outlier_method", "iqr"))
    out["outliers"] = flags
    for m in metrics:
        r = two_way_anova(cells, m, fA, fB, alpha)
        if "error" not in r:
            need_np = (mode == "always") or (mode == "auto" and not r["assumptions_ok"])
            r["nonparametric"] = nonparametric(cells, m, fA, fB, alpha) if need_np else None
        out["metrics"][m] = r
    return out
