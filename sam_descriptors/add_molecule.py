"""새 분자 1개 추가 CLI.

  python add_molecule.py --smiles "FC(F)(F)c1ccc(S)cc1" --name "4-CF3-thiophenol" [--solvent thf|water|gas]
                         [--nconf 20] [--force] [--no-refresh]

  - 결과는 results/database.csv 에 한 줄 추가 (상세 JSON: results/molecules/<name>.json, xyz/, molden/)
  - 같은 정규화 SMILES 가 database.csv 또는 summary.csv 에 이미 있으면 재계산하지 않고 저장값을 출력
    (--force 일 때만 재계산)
  - 추가 후 HOMO–쌍극자 지도(figures/)와 viewer.html 을 다시 만든다 (--no-refresh 로 생략)
"""
import argparse
import csv
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
os.environ.setdefault("OMP_NUM_THREADS", "1")

DB = "results/database.csv"
SHOW = ["name", "smiles", "solvent", "HOMO", "LUMO", "gap", "IP", "EA", "dipole", "dipole_axial",
        "q_anchor_atom", "q_terminal_group", "MW", "TPSA", "logP", "nF", "nRotB", "chain_bonds",
        "zw_dist_boltz", "n_conf_connectivity_changed", "time_s"]


def canonical(smiles):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        sys.exit(f"SMILES 파싱 실패: {smiles}")
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    return Chem.MolToSmiles(m)


def find_existing(canon):
    for path in (DB, "results/summary.csv"):
        if os.path.exists(path):
            for r in csv.DictReader(open(path)):
                if r.get("smiles") == canon:
                    return path, r
    return None, None


def append_db(summary):
    rows = []
    if os.path.exists(DB):
        rows = [r for r in csv.DictReader(open(DB)) if r["name"] != summary["name"]]
    cols = list(summary.keys())
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    os.makedirs("results", exist_ok=True)
    with open(DB, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})
        w.writerow({c: summary.get(c, "") for c in cols})


def show(summary, header):
    print(header)
    for k in SHOW:
        if k in summary and summary[k] not in ("", None):
            v = summary[k]
            try:
                v = f"{float(v):.4g}"
            except (TypeError, ValueError):
                pass
            print(f"  {k:28s} {v}")


def refresh():
    import analyze
    import viewer
    df = analyze.load_summary()
    analyze.fig_homo_dipole_map(df)
    analyze.fig_lumo_dipole_map(df)
    viewer.build()
    print("figures/map_homo_vs_dipole.png, figures/map_lumo_vs_dipole.png, viewer.html 갱신 완료")


def add(smiles, name, solvent=None, nconf=20, force=False, do_refresh=True):
    canon = canonical(smiles)
    path, old = find_existing(canon)
    if old and not force:
        show(old, f"이미 계산된 분자입니다 ({path}, 이름 {old['name']}). 재계산하려면 --force.")
        return old
    from compute import compute_molecule
    t0 = time.time()
    print(f"계산 시작: {name}  {canon}  (conformer {nconf}개)")
    s = compute_molecule(smiles, name, solvent=solvent, n_conf=nconf)
    with open("logs/progress.log", "a") as f:
        f.write(f"[add] {name} | conformer {s['n_conf_ok']}개 | 최저 E {s['E_min_Eh']:.6f} Eh | HOMO {s['HOMO']:.2f} eV | {time.time()-t0:.1f} s\n")
    append_db(s)
    show(s, f"완료 ({time.time()-t0:.1f} s) → {DB} 에 추가")
    if do_refresh:
        refresh()
    return s


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--solvent", choices=["thf", "water", "gas"], default=None,
                    help="생략하면 쯔비터이온→water, 그 외→thf")
    ap.add_argument("--nconf", type=int, default=20)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--no-refresh", action="store_true")
    a = ap.parse_args()
    os.makedirs("logs", exist_ok=True)
    add(a.smiles, a.name, a.solvent, a.nconf, a.force, not a.no_refresh)
