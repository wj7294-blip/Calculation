#!/usr/bin/env python3
"""echem_tools command line.

  python analyze.py drt  FILE [FILE ...] [--meta metadata.csv] [--lambda auto|1e-3] [--out results]
  python analyze.py floating --meta metadata.csv [--data-dir .] [--transient-h 1] [--out results]
  python analyze.py floating FILE --mass-mg 8 --area-cm2 1.13        (single cell quick look)
  python analyze.py validate [--out results_synthetic]               (synthetic end-to-end check)
  python analyze.py all --meta metadata.csv                          (eis + floating from one metadata)
"""
from __future__ import annotations
import argparse
import sys
import warnings
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from echem.config import load_config  # noqa: E402
from echem.io.metadata import load_metadata  # noqa: E402
from echem.report import Report  # noqa: E402


def _cfg(args):
    cfg = load_config(args.config)
    if getattr(args, "lam", None) is not None:
        cfg["eis"]["drt"]["lambda"] = args.lam if args.lam == "auto" else float(args.lam)
    if getattr(args, "transient_h", None) is not None:
        cfg["floating"]["transient_h"] = float(args.transient_h)
    if getattr(args, "leak_window_h", None) is not None:
        cfg["floating"]["leak_window_h"] = float(args.leak_window_h)
    if getattr(args, "v_hold", None) is not None:
        cfg["floating"]["V_hold"] = float(args.v_hold)
    if getattr(args, "inductive", None):
        cfg["eis"]["inductive_tail"] = args.inductive
    if getattr(args, "fmin", None) is not None:
        cfg["eis"]["freq_min_Hz"] = args.fmin
    if getattr(args, "fmax", None) is not None:
        cfg["eis"]["freq_max_Hz"] = args.fmax
    return cfg


def _meta(args, cfg, kind=None):
    meta_path = args.meta or (ROOT / cfg["paths"]["metadata"])
    data_dir = args.data_dir or (ROOT / cfg["paths"]["data_dir"])
    m = load_metadata(meta_path, data_dir)
    if kind:
        m = m[m["type"] == kind]
    return m


def cmd_drt(args):
    from echem.pipeline_eis import run_eis
    cfg = _cfg(args)
    out = Path(args.out or cfg["paths"]["results_dir"]); out.mkdir(parents=True, exist_ok=True)
    entries = []
    if args.files:
        meta = None
        if args.meta:
            try:
                meta = _meta(args, cfg)
            except Exception as e:
                warnings.warn(f"metadata not used: {e}")
        for i, f in enumerate(args.files):
            p = Path(f)
            e = {"label": p.stem, "path": str(p), "profile": args.profile, "order": i}
            if meta is not None:
                hit = meta[meta["file"].apply(lambda x: Path(x).name == p.name)]
                if len(hit):
                    e["label"] = hit.iloc[0]["cell_id"]
                    gk = args.group_by or cfg["eis"].get("compare", {}).get("group_by")
                    if gk in hit.columns:
                        e["group"] = hit.iloc[0][gk]
            entries.append(e)
    else:
        meta = _meta(args, cfg, "eis")
        if not len(meta):
            sys.exit("no EIS files given and metadata has no rows with type=eis")
        gk = args.group_by or cfg["eis"].get("compare", {}).get("group_by")
        for i, m in meta.iterrows():
            e = {"label": m["cell_id"], "path": m["path"], "profile": m.get("profile") if isinstance(m.get("profile"), str) else None, "order": i}
            if gk in meta.columns and gk != "cell_id":
                e["group"] = m[gk]
            entries.append(e)
    rep = Report("EIS → DRT report")
    summ = run_eis(entries, cfg, out, rep, compare_key=args.group_by)
    rep.save(out / "report.html")
    pd.set_option("display.width", 200)
    print(summ.to_string(index=False))
    print(f"\nresults -> {out}/  (report.html, eis_summary.csv, eis_peaks.csv, eis/*.png)")


def cmd_floating(args):
    from echem.pipeline_float import run_floating
    cfg = _cfg(args)
    out = Path(args.out or cfg["paths"]["results_dir"]); out.mkdir(parents=True, exist_ok=True)
    if args.files:
        rows = [{"cell_id": Path(f).stem, "file": f, "path": str(Path(f)), "type": "float", "surface": None, "electrolyte": None,
                 "active_mass_mg": args.mass_mg, "area_cm2": args.area_cm2, "profile": args.profile, "note": ""} for f in args.files]
        meta = pd.DataFrame(rows)
    else:
        meta = _meta(args, cfg, "float")
        if not len(meta):
            sys.exit("metadata has no rows with type=float")
    rep = Report("Floating leakage-current report")
    res = run_floating(meta, cfg, out, rep)
    rep.save(out / "report.html")
    pd.set_option("display.width", 250)
    cols = [c for c in ["cell_id", "surface", "electrolyte", "leak_mA_g", "leak_uA_cm2", "Q_leak_mAh_g", "Q_leak_ex_transient_mAh_g", "decay_slope", "error"] if c in res["cells"]]
    print(res["cells"][cols].to_string(index=False))
    if res["stats"]:
        fA, fB = res["stats"]["factors"]
        for m, r in res["stats"]["metrics"].items():
            if "anova" in r:
                a = r["anova"]
                print(f"{m:32s} p({fA})={a.loc[fA,'PR(>F)']:.3g}  p({fB})={a.loc[fB,'PR(>F)']:.3g}  p(interaction)={a.loc[f'{fA}:{fB}','PR(>F)']:.3g}")
    print(f"\nresults -> {out}/  (report.html, float_cells.csv, float_anova.csv, float/*.png)")


def cmd_validate(args):
    from echem.pipeline_eis import validate_eis
    from echem.pipeline_float import validate_floating
    cfg = _cfg(args)
    out = Path(args.out or "results_synthetic"); out.mkdir(parents=True, exist_ok=True)
    data_root = Path(args.data_dir or cfg["synthetic"].get("out_dir", "data_synthetic"))
    rep = Report("echem_tools — synthetic validation")
    v1 = validate_eis(cfg, out, rep, write_files=data_root)
    v2 = validate_floating(cfg, data_root, out, rep)
    meta = pd.concat([pd.DataFrame(v1["meta_rows"]), v2["meta"].drop(columns=["path"])], ignore_index=True)
    meta.to_csv(data_root / "metadata.csv", index=False)
    rep.save(out / "report.html")
    pd.set_option("display.width", 250)
    print("== DRT peak recovery ==")
    print(v1["recovery"][["circuit", "noise_pct", "element", "tau_true_s", "tau_found_s", "tau_err_decades", "R_true_ohm", "R_found_ohm", "R_err_pct", "found"]].to_string(index=False))
    print("\n== floating ANOVA ==")
    for m, r in v2["stats"]["metrics"].items():
        a = r["anova"]; fA, fB = v2["stats"]["factors"]
        print(f"{m:32s} p({fA})={a.loc[fA,'PR(>F)']:.2e}  p({fB})={a.loc[fB,'PR(>F)']:.2e}  p(interaction)={a.loc[f'{fA}:{fB}','PR(>F)']:.2e}  assumptions_ok={r['assumptions_ok']}")
    print(f"\nsynthetic data -> {data_root}/ (with metadata.csv)\nresults -> {out}/report.html")


def cmd_all(args):
    from echem.pipeline_eis import run_eis
    from echem.pipeline_float import run_floating
    cfg = _cfg(args)
    out = Path(args.out or cfg["paths"]["results_dir"]); out.mkdir(parents=True, exist_ok=True)
    meta = _meta(args, cfg)
    rep = Report("echem_tools report")
    me = meta[meta["type"] == "eis"]
    if len(me):
        gk = cfg["eis"].get("compare", {}).get("group_by")
        entries = [{"label": m["cell_id"], "path": m["path"], "profile": m.get("profile") if isinstance(m.get("profile"), str) else None,
                    "order": i, **({"group": m[gk]} if gk in meta.columns and gk != "cell_id" else {})} for i, m in me.iterrows()]
        run_eis(entries, cfg, out, rep)
    mf = meta[meta["type"] == "float"]
    if len(mf):
        run_floating(mf, cfg, out, rep)
    rep.save(out / "report.html")
    print(f"done: {len(me)} EIS, {len(mf)} floating -> {out}/report.html")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=str(ROOT / "config.yaml"))
    sub = ap.add_subparsers(dest="cmd", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--meta", help="metadata.csv path")
    common.add_argument("--data-dir", help="root for relative paths in metadata (default: config paths.data_dir)")
    common.add_argument("--out", help="results directory")
    common.add_argument("--profile", help="csv profile name (generic|neware|lanhe) for CSV inputs")

    p = sub.add_parser("drt", parents=[common], help="EIS -> DRT for given files (or metadata rows with type=eis)")
    p.add_argument("files", nargs="*")
    p.add_argument("--lambda", dest="lam", help="auto (GCV) or a number")
    p.add_argument("--inductive", choices=["keep", "drop", "fit"])
    p.add_argument("--fmin", type=float); p.add_argument("--fmax", type=float)
    p.add_argument("--group-by", help="metadata column for comparison overlays (e.g. surface, cell_id)")
    p.set_defaults(func=cmd_drt)

    p = sub.add_parser("floating", parents=[common], help="floating leakage current analysis")
    p.add_argument("files", nargs="*")
    p.add_argument("--mass-mg", type=float, help="active mass (single-file mode)")
    p.add_argument("--area-cm2", type=float, help="electrode area (single-file mode)")
    p.add_argument("--transient-h", type=float); p.add_argument("--leak-window-h", type=float); p.add_argument("--v-hold", type=float)
    p.set_defaults(func=cmd_floating)

    p = sub.add_parser("validate", parents=[common], help="generate synthetic data and run both analyses")
    p.add_argument("--lambda", dest="lam"); p.add_argument("--transient-h", type=float)
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("all", parents=[common], help="eis + floating from metadata")
    p.add_argument("--lambda", dest="lam"); p.add_argument("--transient-h", type=float)
    p.set_defaults(func=cmd_all)

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
