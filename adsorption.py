"""MACE(사전학습 uMLIP)로 금속 (111) 표면 위 단일 원자 흡착 에너지 계산.

E_ads = E(slab+X) - E(slab) - E(X bulk, 원자당)

사용법:
    python adsorption.py                                   # 맨 위 MODEL, RESULTS_DIR 사용
    python adsorption.py --model medium-mpa-0 --out results/models/medium-mpa-0
결과는 RESULTS_DIR(기본 results/) 폴더에 저장된다.
"""

# ===== 사용자 설정 =====
METAL = "Cu"                 # 표면 금속 (fcc)
ADSORBATE = "Li"             # 흡착 원자
ADSORBATE_BULK = "bcc"       # 흡착 원자 기준 bulk 결정 구조
SLAB_SIZE = (4, 4)           # 표면 슈퍼셀 크기
N_LAYERS = 4                 # 슬랩 층수
VACUUM = 15.0                # 진공 두께 (Å, ASE fcc111 기준 슬랩 위/아래 각각)
N_FIXED_LAYERS = 2           # 아래에서부터 고정할 층수
SITES = ["fcc", "hcp", "ontop", "bridge"]  # 흡착 자리 (ASE fcc111 이름)
INITIAL_HEIGHT = 2.0         # 흡착 원자 초기 높이 (Å)
FMAX = 0.05                  # BFGS 수렴 기준 (eV/Å)
MODEL = "medium-0b3"         # MACE-MP 체크포인트 (mace_mp model 인자)
DEVICE = "cpu"
RESULTS_DIR = "results"
# =====================

import argparse
import json
import os
import sys
import time
from importlib.metadata import version

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from ase.build import add_adsorbate, bulk, fcc111
from ase.constraints import FixAtoms
from ase.filters import FrechetCellFilter
from ase.io import write
from ase.optimize import BFGS

SITE_LABELS = {"fcc": "fcc hollow", "hcp": "hcp hollow", "ontop": "top", "bridge": "bridge"}


def load_calculator():
    from mace.calculators import mace_mp

    try:
        return mace_mp(model=MODEL, device=DEVICE, default_dtype="float64")
    except Exception as exc:  # 네트워크 차단 등으로 체크포인트를 못 받으면 중단
        print(f"MACE 모델 '{MODEL}' 로드 실패: {exc}", file=sys.stderr)
        sys.exit(1)


def relax_bulk(atoms, calc, logfile):
    """격자상수까지 포함해 bulk를 최적화하고 원자당 에너지를 반환."""
    atoms.calc = calc
    BFGS(FrechetCellFilter(atoms), logfile=logfile).run(fmax=FMAX / 5)
    return atoms, atoms.get_potential_energy() / len(atoms)


def relax(atoms, calc, name):
    atoms.calc = calc
    t0 = time.perf_counter()
    opt = BFGS(atoms, trajectory=os.path.join(RESULTS_DIR, f"{name}.traj"),
               logfile=os.path.join(RESULTS_DIR, "logs", f"{name}.log"))
    opt.run(fmax=FMAX, steps=1000)
    elapsed = time.perf_counter() - t0
    write(os.path.join(RESULTS_DIR, "POSCAR", f"POSCAR_{name}"), atoms, format="vasp", direct=True)
    return atoms.get_potential_energy(), opt.get_number_of_steps(), elapsed


def build_slab(a):
    slab = fcc111(METAL, size=(*SLAB_SIZE, N_LAYERS), a=a, vacuum=VACUUM)
    # fcc111 태그: 맨 위 층 = 1, 맨 아래 층 = N_LAYERS
    fixed = [atom.index for atom in slab if atom.tag > N_LAYERS - N_FIXED_LAYERS]
    slab.set_constraint(FixAtoms(indices=fixed))
    return slab


def main():
    global MODEL, RESULTS_DIR
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", default=MODEL, help="mace_mp model 이름 또는 체크포인트 경로")
    parser.add_argument("--out", default=RESULTS_DIR, help="결과 폴더")
    args = parser.parse_args()
    MODEL, RESULTS_DIR = args.model, args.out

    torch.set_num_threads(os.cpu_count())
    for sub in ("", "logs", "POSCAR"):
        os.makedirs(os.path.join(RESULTS_DIR, sub), exist_ok=True)
    log = os.path.join(RESULTS_DIR, "logs", "bulk.log")

    calc = load_calculator()

    # 1) 기준 bulk: 표면 금속 격자상수, 흡착 원자 원자당 에너지
    t0 = time.perf_counter()
    metal_bulk, _ = relax_bulk(bulk(METAL, "fcc", cubic=True), calc, log)
    a_metal = metal_bulk.cell.lengths()[0]
    ads_bulk, e_ads_bulk = relax_bulk(bulk(ADSORBATE, ADSORBATE_BULK, cubic=True), calc, log)
    t_bulk = time.perf_counter() - t0
    print(f"{METAL} fcc a = {a_metal:.4f} Å, {ADSORBATE} {ADSORBATE_BULK} a = "
          f"{ads_bulk.cell.lengths()[0]:.4f} Å, E({ADSORBATE} bulk)/atom = {e_ads_bulk:.4f} eV "
          f"({t_bulk:.1f} s)")

    # 2) 깨끗한 슬랩
    e_slab, n_slab, t_slab = relax(build_slab(a_metal), calc, "clean_slab")
    print(f"clean slab: E = {e_slab:.4f} eV, {n_slab} steps, {t_slab:.1f} s")

    # 3) 자리별 흡착
    rows = []
    for site in SITES:
        atoms = build_slab(a_metal)
        add_adsorbate(atoms, ADSORBATE, INITIAL_HEIGHT, position=site)
        xy0 = atoms.positions[-1, :2].copy()
        e_tot, n_steps, elapsed = relax(atoms, calc, f"{METAL}111_{ADSORBATE}_{site}")

        top_z = np.mean([atom.z for atom in atoms if atom.tag == 1])
        ads = atoms.positions[-1]
        # 주기 경계를 고려한 수평 이동 거리 (top/bridge가 hollow로 미끄러졌는지 확인용)
        dxy = np.append(ads[:2] - xy0, 0.0)
        frac = np.linalg.solve(atoms.cell.T, dxy)
        frac[:2] -= np.round(frac[:2])
        lateral = np.linalg.norm(atoms.cell.T @ frac)

        e_ads = e_tot - e_slab - e_ads_bulk
        rows.append({
            "site": SITE_LABELS.get(site, site),
            "E_ads (eV)": round(e_ads, 4),
            "height (Å)": round(ads[2] - top_z, 3),
            "lateral shift (Å)": round(lateral, 3),
            "BFGS steps": n_steps,
            "time (s)": round(elapsed, 1),
            "E_total (eV)": round(e_tot, 4),
        })
        print(f"{site:>7}: E_ads = {e_ads:.4f} eV, h = {ads[2] - top_z:.3f} Å, "
              f"shift = {lateral:.3f} Å, {n_steps} steps, {elapsed:.1f} s")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RESULTS_DIR, "adsorption_energies.csv"), index=False)
    print()
    print(df.to_string(index=False))

    # 4) 막대그래프
    fig, ax = plt.subplots(figsize=(5, 3.5))
    bars = ax.bar(df["site"], df["E_ads (eV)"], color="#4C72B0")
    ax.bar_label(bars, fmt="%.3f", padding=3)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("E_ads (eV)")
    ax.set_title(f"{ADSORBATE} on {METAL}(111) {SLAB_SIZE[0]}x{SLAB_SIZE[1]}, MACE {MODEL}")
    ax.margins(y=0.15)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "adsorption_energies.png"), dpi=200)

    # 5) 실행 정보
    info = {
        "model": MODEL,
        "versions": {p: version(p) for p in ("mace-torch", "torch", "ase", "e3nn")},
        "device": DEVICE,
        "threads": torch.get_num_threads(),
        "settings": {
            "metal": METAL, "adsorbate": ADSORBATE, "adsorbate_bulk": ADSORBATE_BULK,
            "slab_size": SLAB_SIZE, "n_layers": N_LAYERS, "vacuum": VACUUM,
            "n_fixed_layers": N_FIXED_LAYERS, "fmax": FMAX,
        },
        "a_metal (Å)": round(a_metal, 4),
        "a_adsorbate_bulk (Å)": round(ads_bulk.cell.lengths()[0], 4),
        "E_adsorbate_bulk_per_atom (eV)": round(e_ads_bulk, 4),
        "E_clean_slab (eV)": round(e_slab, 4),
        "time_bulk (s)": round(t_bulk, 1),
        "time_clean_slab (s)": round(t_slab, 1),
    }
    with open(os.path.join(RESULTS_DIR, "run_info.json"), "w") as f:
        json.dump(info, f, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    main()
