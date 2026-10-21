"""Cu/Ag 결과가 문헌 경향과 반대로 나온 원인 점검.

1) setup    : 슬랩 설정(격자상수, 깨끗한 슬랩 최적화, 원자 수, 층 간격, 셀) 비교
2) surface  : bulk 격자상수와 (111) 표면 에너지
3) subsurf  : Li가 표면 안으로 들어간 구조(1층 치환, 1-2층 사이 침입형)와 hollow 흡착 비교
4) cross    : crosscheck.py 결과(results/crosscheck/*.json)를 한 표로 정리

사용법:
    python diagnose.py setup
    python diagnose.py surface --models medium-0b3 medium-mpa-0 medium-omat-0
    python diagnose.py subsurf --models medium-0b3 medium-mpa-0 medium-omat-0
    python diagnose.py cross
결과는 results/diagnosis/ 에 저장된다.
"""

import argparse
import json
import os

import numpy as np
import pandas as pd
from ase.build import bulk, fcc111
from ase.constraints import FixAtoms
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.optimize import BFGS

METALS = ["Cu", "Ag"]
RUN_DIRS = {
    ("Cu", "medium-0b3"): "results",
    ("Cu", "medium-mpa-0"): "results/models/medium-mpa-0",
    ("Cu", "medium-omat-0"): "results/models/medium-omat-0",
    ("Ag", "medium-0b3"): "results/Ag111",
    ("Ag", "medium-mpa-0"): "results/Ag111/models/medium-mpa-0",
    ("Ag", "medium-omat-0"): "results/Ag111/models/medium-omat-0",
}
OUT = "results/diagnosis"
EV_PER_A2_TO_J_PER_M2 = 16.02177
# adsorption.py(fmax 0.05)에서는 깨끗한 슬랩이 0 step으로 끝나는 경우가 있어, 여기서는 모두 더 엄격하게 수렴
FMAX_TIGHT = 0.01


def mace(model):
    from mace.calculators import mace_mp

    return mace_mp(model=model, device="cpu", default_dtype="float64")


def relax_bulk(el, struct, calc):
    atoms = bulk(el, struct, cubic=True)
    atoms.calc = calc
    BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=0.01)
    return atoms.cell.lengths()[0], atoms.get_potential_energy() / len(atoms)


# ---------- 1) 설정 점검 ----------
def setup():
    rows = []
    for (metal, model), d in RUN_DIRS.items():
        info = json.load(open(os.path.join(d, "run_info.json")))
        traj = read(os.path.join(d, "clean_slab.traj"), ":")
        slab = traj[-1]
        ads = read(os.path.join(d, f"{metal}111_Li_fcc.traj"))
        tags = slab.get_tags()
        z = [slab.positions[tags == t, 2].mean() for t in range(1, info["settings"]["n_layers"] + 1)]
        a = info["a_metal (Å)"]
        fixed_slab = sorted(slab.constraints[0].get_indices()) if slab.constraints else []
        fixed_ads = sorted(ads.constraints[0].get_indices()) if ads.constraints else []
        thickness = max(slab.positions[:, 2]) - min(slab.positions[:, 2])
        rows.append({
            "metal": metal, "model": model,
            "a (Å)": a,
            "a source": "MACE bulk relax (FrechetCellFilter+BFGS)",
            "clean slab relaxed separately": f"yes ({len(traj) - 1} BFGS steps)",
            "fixed atoms slab / slab+Li": f"{len(fixed_slab)} / {len(fixed_ads)}",
            "same fixed indices": fixed_slab == fixed_ads,
            "N atoms slab": len(slab),
            "d12 (Å)": round(z[0] - z[1], 3),
            "d23 (Å)": round(z[1] - z[2], 3),
            "d34 (Å)": round(z[2] - z[3], 3),
            "ideal d (Å)": round(a / np.sqrt(3), 3),
            "d12 change (%)": round((z[0] - z[1]) / (a / np.sqrt(3)) * 100 - 100, 2),
            "|a1|=|a2| (Å)": round(np.linalg.norm(slab.cell[0]), 4),
            "c (Å)": round(slab.cell[2, 2], 3),
            "vacuum gap (Å)": round(slab.cell[2, 2] - thickness, 3),
            "settings": json.dumps({k: v for k, v in info["settings"].items() if k != "metal"}),
        })
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "setup_check.csv"), index=False)
    print(df.drop(columns=["settings"]).to_string(index=False))
    print("settings identical except metal:", df["settings"].nunique() == 1)


# ---------- 2) bulk 격자상수와 표면 에너지 ----------
def surface_energy(metal, a, e_bulk, calc, n_layers):
    """대칭 1x1 슬랩, 모든 원자 최적화. γ = (E_slab - N E_bulk) / 2A"""
    slab = fcc111(metal, size=(1, 1, n_layers), a=a, vacuum=10.0)
    area = np.linalg.norm(np.cross(slab.cell[0], slab.cell[1]))
    slab.calc = calc
    e_unrel = slab.get_potential_energy()
    BFGS(slab, logfile=None).run(fmax=0.005)
    e_rel = slab.get_potential_energy()
    g = lambda e: (e - n_layers * e_bulk) / (2 * area) * EV_PER_A2_TO_J_PER_M2
    return g(e_unrel), g(e_rel)


def surface(models):
    rows = []
    for model in models:
        calc = mace(model)
        for metal in METALS:
            a, e_bulk = relax_bulk(metal, "fcc", calc)
            for n in (5, 7, 9):
                g_unrel, g_rel = surface_energy(metal, a, e_bulk, calc, n)
                rows.append({"model": model, "metal": metal, "a (Å)": round(a, 4),
                             "E_bulk (eV/atom)": round(e_bulk, 4), "layers": n,
                             "gamma unrelaxed (J/m2)": round(g_unrel, 3),
                             "gamma relaxed (J/m2)": round(g_rel, 3)})
            print(rows[-1])
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "surface_energy.csv"), index=False)
    print(df.to_string(index=False))


# ---------- 3) 표면 안으로 들어간 Li ----------
def build_slab(metal, a):
    slab = fcc111(metal, size=(4, 4, 4), a=a, vacuum=15.0)
    slab.set_constraint(FixAtoms(indices=[at.index for at in slab if at.tag > 2]))
    return slab


def nearest_to_center(atoms, mask):
    center = (atoms.cell[0] + atoms.cell[1]) / 2
    idx = np.where(mask)[0]
    d = np.linalg.norm(atoms.positions[idx, :2] - center[:2], axis=1)
    return idx[np.argmin(d)]


def configs(metal, a):
    """(이름, 원자 구조, 1층에서 빠진 금속 원자 수)"""
    out = []
    base = build_slab(metal, a)
    tags = base.get_tags()
    z1 = base.positions[tags == 1, 2].mean()
    z2 = base.positions[tags == 2, 2].mean()
    d = z1 - z2

    # 1층 원자 1개를 Li로 치환
    s = base.copy()
    i = nearest_to_center(s, tags == 1)
    s[i].symbol = "Li"
    out.append(("substitution (layer 1)", s, 1))

    # 1-2층 사이 팔면체 자리: fcc hollow 바로 아래 (3층 원자의 xy)
    s = base.copy()
    j = nearest_to_center(s, tags == 3)
    s += base[[0]]
    s[-1].symbol = "Li"
    s[-1].tag = 0
    s.positions[-1] = [*s.positions[j, :2], z2 + d / 2]
    out.append(("subsurface octahedral (1-2)", s, 0))

    # 1-2층 사이 사면체 자리: hcp hollow 바로 아래 (2층 원자 위, 1층에서 d/4 아래)
    s = base.copy()
    j = nearest_to_center(s, tags == 2)
    s += base[[0]]
    s[-1].symbol = "Li"
    s[-1].tag = 0
    s.positions[-1] = [*s.positions[j, :2], z1 - d / 4]
    out.append(("subsurface tetrahedral (1-2)", s, 0))

    # 비교용: 표면 위 fcc hollow (adsorption.py와 같은 방식)
    from ase.build import add_adsorbate

    s = base.copy()
    add_adsorbate(s, "Li", 2.0, position="fcc")
    out.append(("adsorbed fcc hollow", s, 0))
    return out


def li_depth(atoms, base_top_z):
    li = [at.index for at in atoms if at.symbol == "Li"][0]
    return atoms.positions[li, 2] - base_top_z


def subsurf(models):
    rows = []
    os.makedirs(os.path.join(OUT, "structures"), exist_ok=True)
    for model in models:
        calc = mace(model)
        _, e_li = relax_bulk("Li", "bcc", calc)
        for metal in METALS:
            a, e_m = relax_bulk(metal, "fcc", calc)
            clean = build_slab(metal, a)
            clean.calc = calc
            BFGS(clean, logfile=None).run(fmax=FMAX_TIGHT)
            e_clean = clean.get_potential_energy()
            top_z = clean.positions[clean.get_tags() == 1, 2].mean()
            for name, atoms, n_removed in configs(metal, a):
                atoms.calc = calc
                opt = BFGS(atoms, logfile=None)
                opt.run(fmax=FMAX_TIGHT, steps=2000)
                # E = E(구조) + n_removed·μ_M(bulk) − E(깨끗한 슬랩) − μ_Li(bulk)
                e = atoms.get_potential_energy() + n_removed * e_m - e_clean - e_li
                tag = name.split(" ")[0] + ("_" + name.split(" ")[1] if "subsurface" in name else "")
                write(os.path.join(OUT, "structures", f"POSCAR_{metal}111_{model}_{tag}"),
                      atoms, format="vasp", direct=True, sort=True)
                rows.append({"model": model, "metal": metal, "config": name,
                             "E (eV)": round(e, 4),
                             "Li z - top layer (Å)": round(li_depth(atoms, top_z), 3),
                             "BFGS steps": opt.get_number_of_steps()})
                print(rows[-1])
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "subsurface.csv"), index=False)
    print(df.to_string(index=False))


# ---------- 4) 교차 확인 표 ----------
def cross():
    import glob

    rows = []
    for f in sorted(glob.glob("results/crosscheck/*.json")):
        d = json.load(open(f))
        rows.append({
            "model": d["spec"], "package": f"{d['package']} {d['version']}",
            "a Cu (Å)": d["Cu"]["a (Å)"], "a Ag (Å)": d["Ag"]["a (Å)"],
            "E_ads fcc Cu (eV)": d["Cu"]["E_ads fcc (eV)"],
            "E_ads fcc Ag (eV)": d["Ag"]["E_ads fcc (eV)"],
            "Ag - Cu, fcc (eV)": d["Ag - Cu (eV)"],
            "E_sub Cu (eV)": d["Cu"]["E_sub layer1 (eV)"],
            "E_sub Ag (eV)": d["Ag"]["E_sub layer1 (eV)"],
            "Ag - Cu, most stable (eV)": d["Ag - Cu, most stable (eV)"],
            "time (s)": d["time (s)"],
        })
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "crosscheck_summary.csv"), index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("task", choices=["setup", "surface", "subsurf", "cross"])
    p.add_argument("--models", nargs="+", default=["medium-0b3", "medium-mpa-0", "medium-omat-0"])
    args = p.parse_args()
    os.makedirs(OUT, exist_ok=True)
    {"setup": setup, "surface": lambda: surface(args.models),
     "subsurf": lambda: subsurf(args.models), "cross": cross}[args.task]()
