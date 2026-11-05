"""SMILES → conformer → GFN2-xTB 최적화 → 성질 계산 (분자 1개 단위).

주요 함수 (계산 흐름 순서)
  rdkit_descriptors(molH, info)      : RDKit 2D 기술자
  generate_conformers(molH, n)       : ETKDG conformer 생성 + MMFF 사전 최적화
  xtb_optimize(Z, pos, charge, solv) : tblite(GFN2-xTB, ALPB) + ASE LBFGS 로 구조 최적화
  xtb_properties(Z, pos, ...)        : 최적화 구조에서 HOMO/LUMO/쌍극자/부분전하
  xtb_binary_ipea(...)               : xtb 바이너리로 vertical IP/EA (IPEA-xTB), molden 출력
  check_connectivity(molH, pos)      : 최적화 후 결합 연결 변화 검사
  process_molecule(row, n_conf)      : 위 단계를 모두 묶은 분자 1개 처리 함수 (worker)
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import json
import shutil
import subprocess
import tempfile
import time

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem import rdDetermineBonds  # noqa: F401  (존재 확인용)
from ase import Atoms
from ase.data import covalent_radii
from ase.optimize import LBFGS, FIRE
from tblite.ase import TBLite

HARTREE_EV = 27.211386
AU_DEBYE = 2.541746
BOHR = 0.52917721
KT_KCAL = 0.0019872 * 298.15  # kcal/mol
HARTREE_KCAL = 627.5095

from library import parse_entry  # noqa: E402

XTB_BIN = shutil.which("xtb") or "/opt/xtb/bin/xtb"


def solvent_for(row):
    """시리즈별 ALPB 용매 선택. C 시리즈·cysteine 쯔비터이온 → water, val → 기체상, 나머지 → thf."""
    if row.get("solvent"):
        return None if row["solvent"] in ("gas", "none") else row["solvent"]
    if row["series"] == "val":
        return None
    if row["series"] == "C" or row["name"] == "L-cysteine-zwitterion":
        return "water"
    return "thf"


# ---------------------------------------------------------------- 1) RDKit 기술자
def rdkit_descriptors(molH, info):
    mol = Chem.RemoveHs(molH)
    d = dict(
        MW=Descriptors.MolWt(mol),
        TPSA=rdMolDescriptors.CalcTPSA(mol),
        logP=Descriptors.MolLogP(mol),
        nF=sum(a.GetAtomicNum() == 9 for a in mol.GetAtoms()),
        nRotB=rdMolDescriptors.CalcNumRotatableBonds(mol),
        heavy_atoms=mol.GetNumHeavyAtoms(),
    )
    a, t = info["anchor_idx"], info["terminal_idx"]
    if a is not None and t is not None:
        path = Chem.GetShortestPath(molH, a, t)
        d["chain_bonds"] = len(path) - 1
    else:
        d["chain_bonds"] = None
    return d


def group_atoms(molH, atom_idx, other_idx):
    """atom_idx 를 포함하는 작용기 원자 집합: atom_idx 와 other_idx 를 잇는 최단 경로에서
    atom_idx 의 이웃(스페이서 쪽) 결합을 끊었을 때 atom_idx 쪽 조각."""
    if atom_idx is None or other_idx is None:
        return [atom_idx] if atom_idx is not None else []
    path = Chem.GetShortestPath(molH, atom_idx, other_idx)
    if len(path) < 2:
        return [atom_idx]
    nb = path[1]
    em = Chem.RWMol(molH)
    em.RemoveBond(atom_idx, nb)
    frags = Chem.GetMolFrags(em, asMols=False, sanitizeFrags=False)
    for fr in frags:
        if atom_idx in fr:
            return list(fr)
    return [atom_idx]


# ---------------------------------------------------------------- 2) conformer 생성
def generate_conformers(molH, n_conf=20, seed=42):
    """ETKDGv3 로 n_conf 개 생성 후 MMFF94s 로 사전 최적화. (atomic numbers, [positions...]) 반환."""
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = 0.1
    params.useRandomCoords = True
    cids = list(AllChem.EmbedMultipleConfs(molH, numConfs=n_conf, params=params))
    if len(cids) < n_conf:  # 가지치기로 부족하면 보충
        params.pruneRmsThresh = -1
        params.randomSeed = seed + 1
        cids = list(AllChem.EmbedMultipleConfs(molH, numConfs=n_conf, params=params))
    try:
        AllChem.MMFFOptimizeMoleculeConfs(molH, maxIters=2000, mmffVariant="MMFF94s")
    except Exception:
        pass
    Z = [a.GetAtomicNum() for a in molH.GetAtoms()]
    pos = [molH.GetConformer(c).GetPositions() for c in cids]
    return Z, pos


# ---------------------------------------------------------------- 3) xTB 최적화
def make_calc(charge, solvent):
    kw = dict(method="GFN2-xTB", charge=charge, verbosity=0, max_iterations=400)
    if solvent:
        kw["solvation"] = ("alpb", solvent)
    return TBLite(**kw)


def xtb_optimize(Z, pos, charge=0, solvent="thf", fmax=0.05, steps=800):
    """GFN2-xTB(+ALPB) 에너지/힘으로 ASE LBFGS 최적화. (atoms, energy_eV, converged) 반환."""
    atoms = Atoms(numbers=Z, positions=pos)
    atoms.calc = make_calc(charge, solvent)
    conv = False
    try:
        opt = LBFGS(atoms, logfile=None, maxstep=0.15)
        conv = opt.run(fmax=fmax, steps=steps)
    except Exception:
        pass
    if not conv:
        try:
            atoms.calc = make_calc(charge, solvent)
            opt = FIRE(atoms, logfile=None)
            conv = opt.run(fmax=fmax, steps=steps)
        except Exception:
            pass
    e = atoms.get_potential_energy()
    return atoms, float(e), bool(conv)


# ---------------------------------------------------------------- 4) xTB 성질
def xtb_properties(atoms):
    """현재 계산기의 결과에서 궤도 에너지, 쌍극자, 부분전하를 읽는다."""
    calc = atoms.calc
    res = calc._res if hasattr(calc, "_res") else None
    if res is None:
        atoms.get_potential_energy()
        res = calc._res
    eps = np.array(res.get("orbital-energies")) * HARTREE_EV
    occ = np.array(res.get("orbital-occupations"))
    nocc = int(np.sum(occ > 1.0))
    homo, lumo = float(eps[nocc - 1]), float(eps[nocc])
    dip = np.array(res.get("dipole")) * AU_DEBYE  # Debye, 벡터
    q = np.array(res.get("charges"))
    return dict(HOMO=homo, LUMO=lumo, gap=lumo - homo, dipole_vec=dip.tolist(),
                dipole=float(np.linalg.norm(dip)), charges=q.tolist())


# ---------------------------------------------------------------- 5) 결합 연결 검사
def adjacency_from_geometry(Z, pos, scale=1.25):
    pos = np.asarray(pos)
    r = np.array([covalent_radii[z] for z in Z])
    d = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
    cut = scale * (r[:, None] + r[None, :])
    adj = (d < cut) & ~np.eye(len(Z), dtype=bool)
    return adj


def check_connectivity(molH, pos):
    """RDKit 토폴로지와 기하학적 연결을 비교. 바뀐 결합 목록 반환 (비어 있으면 정상)."""
    Z = [a.GetAtomicNum() for a in molH.GetAtoms()]
    adj = adjacency_from_geometry(Z, pos)
    ref = np.zeros_like(adj)
    for b in molH.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        ref[i, j] = ref[j, i] = True
    changed = []
    sym = [Chem.GetPeriodicTable().GetElementSymbol(z) for z in Z]
    for i in range(len(Z)):
        for j in range(i + 1, len(Z)):
            if adj[i, j] and not ref[i, j]:
                changed.append(f"+{sym[i]}{i}-{sym[j]}{j}")
            elif ref[i, j] and not adj[i, j]:
                changed.append(f"-{sym[i]}{i}-{sym[j]}{j}")
    return changed


# ---------------------------------------------------------------- 6) xtb 바이너리: IPEA, molden
def write_xyz(path, Z, pos, comment=""):
    pt = Chem.GetPeriodicTable()
    with open(path, "w") as f:
        f.write(f"{len(Z)}\n{comment}\n")
        for z, p in zip(Z, pos):
            f.write(f"{pt.GetElementSymbol(int(z)):2s} {p[0]:14.8f} {p[1]:14.8f} {p[2]:14.8f}\n")


def _parse_float(lines, key):
    for ln in lines:
        if key in ln:
            try:
                return float(ln.split(":")[-1].split()[0])
            except Exception:
                pass
    return None


def xtb_binary_ipea(Z, pos, charge, solvent, molden_out=None, workdir=None):
    """xtb 바이너리로 (1) IPEA-xTB vertical IP/EA, (2) GFN2 molden 파일을 생성."""
    out = dict(IP=None, EA=None, xtb_HOMO=None, xtb_LUMO=None)
    if not os.path.exists(XTB_BIN):
        return out
    env = dict(os.environ, OMP_NUM_THREADS="1", XTBPATH="/opt/xtb/share/xtb")
    tmp = tempfile.mkdtemp(prefix="xtb_", dir=workdir)
    try:
        xyz = os.path.join(tmp, "mol.xyz")
        write_xyz(xyz, Z, pos)
        base = [XTB_BIN, "mol.xyz", "--chrg", str(charge), "--norestart"]
        solv = ["--alpb", solvent] if solvent else []
        r = subprocess.run(base + ["--vipea"] + solv, cwd=tmp, env=env,
                           capture_output=True, text=True, timeout=600)
        lines = r.stdout.splitlines()
        out["IP"] = _parse_float(lines, "delta SCC IP (eV):")
        out["EA"] = _parse_float(lines, "delta SCC EA (eV):")
        if molden_out:
            for junk in ("xtbrestart", "charges", "wbo"):
                p = os.path.join(tmp, junk)
                if os.path.exists(p):
                    os.remove(p)
            r = subprocess.run(base + ["--gfn", "2", "--molden"] + solv, cwd=tmp, env=env,
                               capture_output=True, text=True, timeout=600)
            for ln in r.stdout.splitlines():
                if "(HOMO)" in ln:
                    out["xtb_HOMO"] = float(ln.split()[-2])
                elif "(LUMO)" in ln:
                    out["xtb_LUMO"] = float(ln.split()[-2])
            src = os.path.join(tmp, "molden.input")
            if os.path.exists(src):
                shutil.copy(src, molden_out)
    except Exception as e:  # noqa
        out["error"] = str(e)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return out


# ---------------------------------------------------------------- 7) 분자 1개 전체 처리
def process_molecule(row, n_conf=20, seed=42, outdir="results/molecules", xyzdir="xyz",
                     moldendir="molden", do_ipea=True, tag=""):
    t0 = time.time()
    name = row["name"] + tag
    molH, info = parse_entry(row)
    charge = info["charge"]
    solvent = solvent_for(row)
    desc = rdkit_descriptors(molH, info)
    Z, confs = generate_conformers(molH, n_conf, seed)

    a_idx, t_idx = info["anchor_idx"], info["terminal_idx"]
    anchor_group = group_atoms(molH, a_idx, t_idx)
    terminal_group = group_atoms(molH, t_idx, a_idx)

    conf_results = []
    for ci, pos in enumerate(confs):
        try:
            atoms, e, conv = xtb_optimize(Z, pos, charge, solvent)
            props = xtb_properties(atoms)
            fpos = atoms.get_positions()
            changed = check_connectivity(molH, fpos)
            conf_results.append(dict(idx=ci, energy_eV=e, converged=conv, positions=fpos.tolist(),
                                     connectivity_changed=changed, **props))
        except Exception as ex:  # noqa
            conf_results.append(dict(idx=ci, energy_eV=None, error=str(ex)))

    ok = [c for c in conf_results if c.get("energy_eV") is not None]
    intact = [c for c in ok if not c["connectivity_changed"]]
    pool = intact if intact else ok
    best = min(pool, key=lambda c: c["energy_eV"])
    emin = best["energy_eV"]
    for c in ok:
        c["rel_kcal"] = (c["energy_eV"] - emin) * 23.0605
        c["boltz_w"] = float(np.exp(-c["rel_kcal"] / KT_KCAL)) if not c["connectivity_changed"] else 0.0
    wsum = sum(c["boltz_w"] for c in ok) or 1.0
    for c in ok:
        c["boltz_w"] /= wsum

    fpos = np.array(best["positions"])
    q = np.array(best["charges"])
    dip = np.array(best["dipole_vec"])
    # 앵커→말단 단위 벡터와 쌍극자 성분
    if a_idx is not None and t_idx is not None:
        u = fpos[t_idx] - fpos[a_idx]
        u = u / np.linalg.norm(u)
        dip_axial = float(np.dot(dip, u))
        anchor_terminal_dist = float(np.linalg.norm(fpos[t_idx] - fpos[a_idx]))
    else:
        dip_axial, anchor_terminal_dist = None, None

    # 쯔비터이온 N+ – 음이온 중심 거리 (conformer 전체)
    zw = None
    n_idx = info["nplus_idx"]
    if n_idx is not None and t_idx is not None:
        dists = []
        for c in ok:
            p = np.array(c["positions"])
            c["zw_dist"] = float(np.linalg.norm(p[n_idx] - p[t_idx]))
            dists.append(c["zw_dist"])
        w = np.array([c["boltz_w"] for c in ok])
        dists = np.array(dists)
        zw = dict(best=float(np.linalg.norm(fpos[n_idx] - fpos[t_idx])),
                  boltz_mean=float(np.sum(w * dists)), min=float(dists.min()), max=float(dists.max()),
                  std=float(np.sqrt(np.sum(w * (dists - np.sum(w * dists)) ** 2))))

    # xyz / molden / IPEA
    os.makedirs(xyzdir, exist_ok=True)
    os.makedirs(moldendir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)
    write_xyz(os.path.join(xyzdir, f"{name}.xyz"), Z, fpos,
              f"{name} GFN2-xTB/ALPB({solvent}) E={emin/HARTREE_EV:.8f} Eh")
    ipea = {}
    if do_ipea:
        ipea = xtb_binary_ipea(Z, fpos, charge, solvent,
                               molden_out=os.path.join(moldendir, f"{name}.molden"))

    summary = dict(
        name=name, series=row["series"], anchor=row["anchor"], terminal=row["terminal"],
        spacer=row["spacer"], smiles=info["canonical_smiles"], charge=charge, solvent=solvent or "gas",
        n=row.get("n", ""), m=row.get("m", ""), exp_ip=row.get("exp_ip", ""),
        n_conf=len(confs), n_conf_ok=len(ok), n_conf_connectivity_changed=len(ok) - len(intact),
        best_conf=best["idx"], E_min_Eh=emin / HARTREE_EV, best_converged=best["converged"],
        HOMO=best["HOMO"], LUMO=best["LUMO"], gap=best["gap"], dipole=best["dipole"],
        dipole_axial=dip_axial, anchor_terminal_dist=anchor_terminal_dist,
        q_anchor_atom=float(q[a_idx]) if a_idx is not None else None,
        q_anchor_group=float(q[anchor_group].sum()) if anchor_group else None,
        q_terminal_atom=float(q[t_idx]) if t_idx is not None else None,
        q_terminal_group=float(q[terminal_group].sum()) if terminal_group else None,
        IP=ipea.get("IP"), EA=ipea.get("EA"), xtb_HOMO=ipea.get("xtb_HOMO"), xtb_LUMO=ipea.get("xtb_LUMO"),
        zw_dist_best=zw["best"] if zw else None, zw_dist_boltz=zw["boltz_mean"] if zw else None,
        zw_dist_std=zw["std"] if zw else None, zw_dist_min=zw["min"] if zw else None,
        zw_dist_max=zw["max"] if zw else None,
        time_s=time.time() - t0, **desc,
    )
    detail = dict(summary=summary, Z=Z, anchor_idx=a_idx, terminal_idx=t_idx, nplus_idx=n_idx,
                  anchor_group=anchor_group, terminal_group=terminal_group,
                  conformers=sorted(ok, key=lambda c: c["energy_eV"]),
                  failed=[c for c in conf_results if c.get("energy_eV") is None])
    with open(os.path.join(outdir, f"{name}.json"), "w") as f:
        json.dump(detail, f)
    return summary


# ---------------------------------------------------------------- 8) 분자 1개를 추가하는 공개 API
ANCHOR_SMARTS = [("SH", "[SX2;H1]"), ("PO3H2", "[PX4](=O)([OH])[OH]"), ("COOH", "[CX3](=O)[OH]"),
                 ("S-", "[SX1-]"), ("COO-", "[CX3](=O)[O-]")]


def auto_annotate(smiles):
    """atom-map 이 없는 SMILES 에 앵커(1)/말단(2)/N+(4) map 을 자동으로 붙인다.
    앵커: 티올 S > 포스폰산 P > 카복실 C 순서로 탐색. 말단: 앵커에서 결합 수가 가장 먼 무거운 원자."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit 파싱 실패: {smiles}")
    if any(a.GetAtomMapNum() for a in mol.GetAtoms()):
        return smiles, "user-mapped"
    anchor_kind, a_idx = None, None
    for kind, sma in ANCHOR_SMARTS:
        hit = mol.GetSubstructMatches(Chem.MolFromSmarts(sma))
        if hit:
            anchor_kind, a_idx = kind, hit[0][0]
            break
    if a_idx is None:
        a_idx, anchor_kind = 0, "atom0"
    dm = Chem.GetDistanceMatrix(mol)
    far = max(range(mol.GetNumAtoms()), key=lambda i: (dm[a_idx, i], mol.GetAtomWithIdx(i).GetAtomicNum()))
    # CF3, SO3H, COOH, PO3H2 처럼 중심 원자에 말단 원자가 여러 개 달린 작용기는 중심 원자를 말단으로 삼는다
    fa = mol.GetAtomWithIdx(far)
    if fa.GetDegree() == 1:
        nb = fa.GetNeighbors()[0]
        if nb.GetIdx() != a_idx and sum(x.GetDegree() == 1 for x in nb.GetNeighbors()) >= 2:
            far = nb.GetIdx()
    mol.GetAtomWithIdx(a_idx).SetAtomMapNum(1)
    if far != a_idx:
        mol.GetAtomWithIdx(far).SetAtomMapNum(2)
    for a in mol.GetAtoms():
        if a.GetFormalCharge() > 0 and a.GetAtomicNum() == 7:
            a.SetAtomMapNum(4)
            break
    return Chem.MolToSmiles(mol), anchor_kind


def make_row(smiles, name, solvent=None, series="user"):
    """compute_molecule 용 row(dict) 생성. 용매를 주지 않으면 쯔비터이온→water, 그 외→thf."""
    mapped, anchor_kind = auto_annotate(smiles)
    mol = Chem.MolFromSmiles(mapped)
    zw = any(a.GetFormalCharge() > 0 for a in mol.GetAtoms()) and any(a.GetFormalCharge() < 0 for a in mol.GetAtoms())
    if solvent is None:
        solvent = "water" if zw else "thf"
    return dict(name=name, smiles=mapped, anchor=anchor_kind, spacer="-", terminal="auto",
                series=series, n="", m="", exp_ip="", solvent=solvent)


def compute_molecule(smiles, name, solvent=None, n_conf=20, seed=42, **kw):
    """분자 하나(SMILES, 이름)를 받아 2–4 단계(RDKit 기술자 → conformer → xTB 최적화 → 성질)를
    모두 수행하고 summary dict 를 반환한다. solvent: 'thf' | 'water' | 'gas' | None(자동)."""
    row = make_row(smiles, name, solvent)
    return process_molecule(row, n_conf=n_conf, seed=seed, **kw)
