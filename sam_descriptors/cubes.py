"""molden/*.molden → cubes/*_homo.cube, *_lumo.cube (pyscf 로 궤도를 격자에 평가).

xtb 의 molden 파일은 primitive 계수에 정규화 상수가 이미 곱해져 있어 pyscf 가 다시 정규화하면
H/C/N/O/F 궤도가 왜곡된다. 여기서는 그 상수를 나눠 준 뒤 읽는다(s, p 는 정확; S, P 의 3s/3p/3d 는
xtb 내부 기저와 Gaussian 전개가 달라 근사적). 궤도 모양은 경향 비교용이다.
"""
import copy
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

REPRESENTATIVE = ["SH-Ph-H", "SH-Ph-CF3", "SH-Ph-NH2", "SH-Ph-CN", "C6F5-thiol", "PO3H2-Ph-H",
                  "COOH-Ph-H", "SH-C3-H", "ZW-SO3-n2-m3", "L-cysteine-zwitterion"]


def load_xtb_molden(path):
    from pyscf.tools import molden
    from pyscf import gto
    txt = open(path).read()
    if "[5D]" not in txt:
        txt = txt.replace("[GTO]", "[5D]\n[GTO]")
    tmp = path + ".5d"
    open(tmp, "w").write(txt)
    mol, mo_energy, mo_coeff, mo_occ, irrep, spins = molden.load(tmp)
    os.remove(tmp)
    basis = copy.deepcopy(mol._basis)
    for sym, shells in basis.items():
        for sh in shells:
            l = sh[0]
            for prim in sh[1:]:
                prim[1] /= gto.gto_norm(l, prim[0])
    mol._basis = basis
    mol.basis = basis
    mol.build(verbose=0)
    S = mol.intor("int1e_ovlp")
    for i in range(mo_coeff.shape[1]):  # 표시용 재정규화
        n = mo_coeff[:, i] @ S @ mo_coeff[:, i]
        if n > 0:
            mo_coeff[:, i] /= np.sqrt(n)
    return mol, mo_energy, mo_coeff, mo_occ


def write_cubes(name, moldendir="molden", cubedir="cubes", spacing=0.35, margin=3.5, nmax=48):
    from pyscf.tools import cubegen
    path = os.path.join(moldendir, f"{name}.molden")
    if not os.path.exists(path):
        return None
    mol, e, c, occ = load_xtb_molden(path)
    nocc = int(round(occ.sum() / 2))
    coords = mol.atom_coords() * 0.52917721
    ext = coords.max(0) - coords.min(0) + 2 * margin
    n = [int(min(nmax, max(24, round(x / spacing)))) for x in ext]
    os.makedirs(cubedir, exist_ok=True)
    out = {}
    for label, idx in (("homo", nocc - 1), ("lumo", nocc)):
        f = os.path.join(cubedir, f"{name}_{label}.cube")
        cubegen.orbital(mol, f, c[:, idx], nx=n[0], ny=n[1], nz=n[2], margin=margin / 0.52917721)
        out[label] = f
    return out


if __name__ == "__main__":
    os.chdir(HERE)
    names = sys.argv[1:] or REPRESENTATIVE
    for nm in names:
        r = write_cubes(nm)
        print(nm, "->", r)
