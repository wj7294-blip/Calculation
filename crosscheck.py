"""여러 uMLIP으로 Cu(111), Ag(111) fcc hollow Li 흡착 에너지(와 1층 치환 에너지) 교차 확인.

adsorption.py와 같은 구조(4x4, 4층, 진공 15 Å, 아래 2층 고정, Li bcc bulk 기준)를 쓰되,
모든 최적화를 fmax 0.01 eV/Å로 수렴시킨다.
치환 에너지: E_sub = E(1층 원자 1개를 Li로 바꾼 슬랩) + μ_M(bulk) − E(깨끗한 슬랩) − μ_Li(bulk)
ase/numpy만 쓰므로 모델별 가상환경에서 실행한다.

사용법:
    python crosscheck.py mace:medium-0b3
    python crosscheck.py sevennet:7net-omni
    python crosscheck.py orb:orb_v3_conservative_inf_omat
    python crosscheck.py mattersim:MatterSim-v1.0.0-5M.pth
결과는 results/crosscheck/<family>_<model>.json 에 저장된다.
"""

import json
import os
import sys
import time
from importlib.metadata import version

import numpy as np
from ase.build import add_adsorbate, bulk, fcc111
from ase.constraints import FixAtoms
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS

METALS = ["Cu", "Ag"]
FMAX = 0.01
OUT = "results/crosscheck"


def get_calc(spec):
    family, model = spec.split(":", 1)
    if family == "mace":
        from mace.calculators import mace_mp

        return mace_mp(model=model, device="cpu", default_dtype="float64"), "mace-torch"
    if family == "sevennet":
        from sevenn.calculator import SevenNetCalculator

        kwargs = {"modal": "mpa"} if model in ("7net-omni", "7net-mf-ompa") else {}
        return SevenNetCalculator(model=model, device="cpu", **kwargs), "sevenn"
    if family == "orb":
        from orb_models.forcefield import pretrained
        from orb_models.forcefield.inference.calculator import ORBCalculator

        orbff, adapter = getattr(pretrained, model)(device="cpu", precision="float32-highest")
        return ORBCalculator(orbff, adapter, device="cpu"), "orb-models"
    if family == "mattersim":
        from mattersim.forcefield import MatterSimCalculator

        return MatterSimCalculator(load_path=model, device="cpu"), "mattersim"
    raise ValueError(family)


def relax_bulk(el, struct, calc):
    atoms = bulk(el, struct, cubic=True)
    atoms.calc = calc
    BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=FMAX)
    return atoms.cell.lengths()[0], atoms.get_potential_energy() / len(atoms)


def relax(atoms, calc):
    atoms.calc = calc
    opt = BFGS(atoms, logfile=None)
    opt.run(fmax=FMAX, steps=2000)
    return atoms.get_potential_energy(), opt.get_number_of_steps()


def main(spec):
    t0 = time.perf_counter()
    calc, pkg = get_calc(spec)
    a_li, e_li = relax_bulk("Li", "bcc", calc)
    out = {"spec": spec, "package": pkg, "version": version(pkg), "fmax": FMAX,
           "a_Li_bcc (Å)": round(a_li, 4), "E_Li_bulk (eV/atom)": round(e_li, 4)}
    for metal in METALS:
        a, e_m = relax_bulk(metal, "fcc", calc)
        slab = fcc111(metal, size=(4, 4, 4), a=a, vacuum=15.0)
        # ORB는 일부 축만 주기적인 셀을 지원하지 않는다. 진공(30 Å)이 cutoff보다 훨씬 커서 결과는 같다.
        slab.pbc = True
        slab.set_constraint(FixAtoms(indices=[at.index for at in slab if at.tag > 2]))
        ads = slab.copy()
        sub = slab.copy()
        e_slab, n_slab = relax(slab, calc)
        add_adsorbate(ads, "Li", 2.0, position="fcc")
        e_ads_tot, n_ads = relax(ads, calc)
        top_z = ads.positions[ads.get_tags() == 1, 2].mean()
        center = (sub.cell[0] + sub.cell[1]) / 2
        top = np.where(sub.get_tags() == 1)[0]
        i = top[np.argmin(np.linalg.norm(sub.positions[top, :2] - center[:2], axis=1))]
        sub[i].symbol = "Li"
        e_sub_tot, n_sub = relax(sub, calc)
        out[metal] = {
            "a (Å)": round(a, 4),
            "E_ads fcc (eV)": round(e_ads_tot - e_slab - e_li, 4),
            "height (Å)": round(ads.positions[-1, 2] - top_z, 3),
            "E_sub layer1 (eV)": round(e_sub_tot + e_m - e_slab - e_li, 4),
            "BFGS steps clean/ads/sub": [n_slab, n_ads, n_sub],
        }
        print(metal, out[metal], flush=True)
    out["Ag - Cu (eV)"] = round(out["Ag"]["E_ads fcc (eV)"] - out["Cu"]["E_ads fcc (eV)"], 4)
    best = {m: min(out[m]["E_ads fcc (eV)"], out[m]["E_sub layer1 (eV)"]) for m in METALS}
    out["Ag - Cu, most stable (eV)"] = round(best["Ag"] - best["Cu"], 4)
    out["time (s)"] = round(time.perf_counter() - t0, 1)
    os.makedirs(OUT, exist_ok=True)
    name = spec.replace(":", "_").replace("/", "_").removesuffix(".pth")
    with open(os.path.join(OUT, f"{name}.json"), "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False, default=float)  # numpy float32/64 → float
    print(json.dumps(out, indent=2, ensure_ascii=False, default=float))


if __name__ == "__main__":
    main(sys.argv[1])
