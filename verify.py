"""adsorption.py 결과 검증 및 모델 비교.

1) 자리별 Li-Cu 최근접 거리 3개, 최적화 전후 수평 이동, 최종 자리 판정
2) 높이 정의 비교 (맨 위 층 평균 / 배위 Cu 평균 / 바로 아래 Cu)
3) 기준 격자상수가 MACE 최적화 값인지 확인 (E(a) 스캔)
4) 여러 모델의 E_ads 비교

사용법:
    python verify.py results results/models/medium-mpa-0 results/models/medium-omat-0
"""

import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ase.build import bulk
from ase.geometry import find_mic
from ase.io import read

SITES = ["fcc", "hcp", "ontop", "bridge"]
SITE_LABELS = {"fcc": "fcc hollow", "hcp": "hcp hollow", "ontop": "top", "bridge": "bridge"}
N_COORD = {"top": 1, "bridge": 2, "fcc hollow": 3, "hcp hollow": 3}
EXP_A = {"Cu": 3.615, "Li": 3.51}  # 실온 실험 격자상수 (Å)


def xy_mic(vecs, cell):
    vecs = np.atleast_2d(np.array(vecs, dtype=float))
    vecs[:, 2] = 0.0
    return find_mic(vecs, cell, pbc=[True, True, False])


def classify_site(atoms, li):
    """최종 Li의 xy 위치에서 가장 가까운 기준점(top/bridge/fcc/hcp)을 찾는다."""
    tags = atoms.get_tags()
    pos = atoms.positions
    top = pos[tags == 1]
    # ASE fcc111: hcp 자리 = 2층 원자 위, fcc 자리 = 3층 원자 위
    refs = {"top": top, "hcp hollow": pos[tags == 2], "fcc hollow": pos[tags == 3]}
    d_nn = np.sort(np.linalg.norm(xy_mic(top[1:] - top[0], atoms.cell)[0], axis=1))[0]
    bridges = []
    for p in top:
        v, d = xy_mic(top - p, atoms.cell)
        for j in np.where(np.abs(d - d_nn) < 0.1)[0]:
            bridges.append(p + v[j] / 2)
    refs["bridge"] = np.array(bridges)
    best = {}
    for name, pts in refs.items():
        _, d = xy_mic(pts - li, atoms.cell)
        best[name] = d.min()
    site = min(best, key=best.get)
    return site, best[site]


def analyze(run_dir):
    rows = []
    for site in SITES:
        traj = read(os.path.join(run_dir, f"Cu111_Li_{site}.traj"), ":")
        init, final = traj[0], traj[-1]
        li = len(final) - 1
        cu = [a.index for a in final if a.symbol == "Cu"]
        d_cu = np.sort(final.get_distances(li, cu, mic=True))
        _, shift = xy_mic(final.positions[li] - init.positions[li], final.cell)
        label = SITE_LABELS[site]
        final_site, off = classify_site(final, final.positions[li])

        tags = final.get_tags()
        z_li = final.positions[li, 2]
        h_layer = z_li - final.positions[tags == 1, 2].mean()
        # 배위 Cu: 판정된 자리의 배위수만큼 가장 가까운 Cu
        order = np.argsort(final.get_distances(li, cu, mic=True))
        coord = [cu[k] for k in order[:N_COORD[final_site]]]
        h_coord = z_li - final.positions[coord, 2].mean()
        # 바로 아래 Cu: xy가 가장 가까운 Cu, 같은 xy에 여러 층이 겹치면 가장 위 원자
        _, dxy = xy_mic(final.positions[cu] - final.positions[li], final.cell)
        cand = [cu[k] for k in np.where(dxy < dxy.min() + 0.1)[0]]
        below = max(cand, key=lambda k: final.positions[k, 2])
        h_below = z_li - final.positions[below, 2]

        rows.append({
            "site": label,
            "d1 (Å)": round(d_cu[0], 3), "d2 (Å)": round(d_cu[1], 3), "d3 (Å)": round(d_cu[2], 3),
            "lateral shift (Å)": round(float(shift[0]), 3),
            "final site": final_site,
            "xy offset from site (Å)": round(off, 3),
            "slid?": "yes" if final_site != label else "no",
            "h_layer (Å)": round(h_layer, 3),
            "h_coord (Å)": round(h_coord, 3),
            "Cu below: layer": int(tags[below]),
            "Cu below: xy dist (Å)": round(float(dxy.min()), 3),
            "h_below (Å)": round(h_below, 3),
        })
    return pd.DataFrame(rows)


def perturb_test(run_dir, calc, info, shift=0.2, fmax=0.01):
    """대칭으로 수평 힘이 0인 문제를 피하려고 Li를 shift(Å)만큼 밀고 재최적화.

    hollow/top은 대칭축과 겹치지 않는 방향(x축에서 17°)으로 민다.
    bridge는 두 Cu를 잇는 축에 수직인 양쪽 방향(각각 fcc/hcp hollow 쪽)으로 민다.
    """
    from ase.optimize import BFGS

    e_slab = info["E_clean_slab (eV)"]
    e_bulk = info["E_adsorbate_bulk_per_atom (eV)"]
    cases = []
    for site in SITES:
        atoms = read(os.path.join(run_dir, f"Cu111_Li_{site}.traj"))
        li = len(atoms) - 1
        if site == "bridge":
            cu = [a.index for a in atoms if a.tag == 1]
            pair = np.array(cu)[np.argsort(atoms.get_distances(li, cu, mic=True))[:2]]
            axis, _ = xy_mic(atoms.positions[pair[1]] - atoms.positions[pair[0]], atoms.cell)
            perp = np.cross([0, 0, 1], axis[0])
            perp /= np.linalg.norm(perp)
            cases += [(site, "+perp", perp), (site, "-perp", -perp)]
        else:
            cases.append((site, "17deg", np.array([np.cos(np.radians(17)), np.sin(np.radians(17)), 0.0])))

    rows = []
    for site, tag, direction in cases:
        atoms = read(os.path.join(run_dir, f"Cu111_Li_{site}.traj"))
        li = len(atoms) - 1
        start = atoms.positions[li].copy()
        atoms.positions[li] += shift * direction
        atoms.calc = calc
        opt = BFGS(atoms, logfile=None,
                   trajectory=os.path.join(run_dir, f"perturbed_Cu111_Li_{site}_{tag}.traj"))
        opt.run(fmax=fmax, steps=1000)
        _, moved = xy_mic(atoms.positions[li] - start, atoms.cell)
        final_site, off = classify_site(atoms, atoms.positions[li])
        rows.append({
            "start site": SITE_LABELS[site],
            "push": tag,
            "final site (perturbed)": final_site,
            "xy offset from site (Å)": round(off, 3),
            "moved from start (Å)": round(float(moved[0]), 3),
            "E_ads (eV)": round(atoms.get_potential_energy() - e_slab - e_bulk, 4),
            "BFGS steps": opt.get_number_of_steps(),
        })
    return pd.DataFrame(rows)


def check_lattice(run_dir):
    """run_info의 격자상수가 해당 모델 에너지 최소점인지 ±1% 스캔으로 확인."""
    from mace.calculators import mace_mp

    info = json.load(open(os.path.join(run_dir, "run_info.json")))
    s = info["settings"]
    calc = mace_mp(model=info["model"], device="cpu", default_dtype="float64")
    slab = read(os.path.join(run_dir, "POSCAR", "POSCAR_clean_slab"))
    a_from_slab = np.linalg.norm(slab.cell[0]) / s["slab_size"][0] * np.sqrt(2)
    out = []
    for el, struct, a0 in [(s["metal"], "fcc", info["a_metal (Å)"]),
                           (s["adsorbate"], s["adsorbate_bulk"], info["a_adsorbate_bulk (Å)"])]:
        e = []
        for f in (0.99, 1.0, 1.01):
            b = bulk(el, struct, a=a0 * f, cubic=True)
            b.calc = calc
            e.append(b.get_potential_energy() / len(b))
        out.append({
            "element": f"{el} ({struct})", "a used (Å)": a0, "a exp (Å)": EXP_A.get(el),
            "E(0.99a)-E(a) (meV/atom)": round((e[0] - e[1]) * 1000, 2),
            "E(1.01a)-E(a) (meV/atom)": round((e[2] - e[1]) * 1000, 2),
            "minimum at a?": "yes" if e[1] < e[0] and e[1] < e[2] else "no",
        })
    out[0]["a from slab POSCAR (Å)"] = round(a_from_slab, 4)
    return pd.DataFrame(out), calc, info


def main(run_dirs):
    pd.set_option("display.width", 250)
    comp = {}
    for d in run_dirs:
        model = json.load(open(os.path.join(d, "run_info.json")))["model"]
        print(f"\n===== {model} ({d}) =====")
        geo = analyze(d)
        geo.to_csv(os.path.join(d, "geometry_check.csv"), index=False)
        print(geo.to_string(index=False))
        lat, calc, info = check_lattice(d)
        lat.to_csv(os.path.join(d, "lattice_check.csv"), index=False)
        print(lat.to_string(index=False))
        pert = perturb_test(d, calc, info)
        pert.to_csv(os.path.join(d, "perturbation_check.csv"), index=False)
        print(pert.to_string(index=False))
        e = pd.read_csv(os.path.join(d, "adsorption_energies.csv"))
        comp[model] = e.set_index("site")["E_ads (eV)"]

    df = pd.DataFrame(comp)
    df.index.name = "site"
    os.makedirs("results", exist_ok=True)
    df.to_csv("results/model_comparison.csv")
    print("\n===== E_ads (eV) by model =====")
    print(df.to_string())

    ax = df.plot.bar(figsize=(7, 3.8), rot=0, width=0.8)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("E_ads (eV)")
    ax.set_xlabel("")
    ax.set_title("Li on Cu(111) 4x4: E_ads by MACE model")
    ax.legend(title="model", fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    plt.tight_layout()
    plt.savefig("results/model_comparison.png", dpi=200)


if __name__ == "__main__":
    main(sys.argv[1:] or ["results"])
