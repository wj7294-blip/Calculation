"""분자 라이브러리 생성 (library.csv).

SMILES 조립 규칙
  - 앵커 원자에는 atom-map 1, 말단 원자에는 atom-map 2, 말단이 붙는 스페이서 원자에는
    atom-map 3, 쯔비터이온의 N+ 에는 atom-map 4 를 붙인다.
  - RDKit 로 파싱/정규화한 뒤 map 번호로 원자 인덱스를 기록하고 map 은 제거한다.
"""
import csv
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

ANCHORS = {
    "SH": "[SH:1]",
    "PO3H2": "O[P:1](=O)(O)",   # 접두 조각: 다음 원자가 P 에 결합하도록 P 를 마지막에 둔다
    "COOH": "O[C:1](=O)",
}
# 스페이서: {A} 자리에 앵커, {T} 자리에 말단. 말단이 -H 인 경우는 별도 처리.
SPACERS = {
    "Ph": ("{A}c1cc[c:3]({T})cc1", "{A}c1cc[cH:3]cc1"),
    "C3": ("{A}CC[CH2:3]{T}", "{A}CC[CH3:3]"),
}
TERMINALS = {
    "H": "",
    "CH3": "[CH3:2]",
    "F": "[F:2]",
    "CF3": "[C:2](F)(F)F",
    "CN": "[C:2]#N",
    "OCH3": "[O:2]C",
    "COOH": "[C:2](=O)O",
    "SO3H": "[S:2](=O)(=O)O",
    "NH2": "[NH2:2]",
}
SPACER_LABEL = {"Ph": "para-phenylene", "C3": "(CH2)3"}


def series_a():
    rows = []
    for a, asmi in ANCHORS.items():
        for sp, (tmpl, tmpl_h) in SPACERS.items():
            for t, tsmi in TERMINALS.items():
                smi = tmpl_h.format(A=asmi) if t == "H" else tmpl.format(A=asmi, T=tsmi)
                rows.append(dict(name=f"{a}-{sp}-{t}", smiles=smi, anchor=a,
                                 spacer=SPACER_LABEL[sp], terminal=t, series="A"))
    return rows


def series_b():
    rows = []
    for a, asmi in ANCHORS.items():
        if a == "COOH":
            continue
        tag = "thiol" if a == "SH" else "PO3H2"
        rows += [
            dict(name=f"C6F5-{tag}", smiles=f"{asmi}c1c(F)c(F)[c:3]([F:2])c(F)c1F",
                 anchor=a, spacer="para-phenylene", terminal="F(penta)", series="B"),
            dict(name=f"4-CF3-C6H4-{tag}", smiles=f"{asmi}c1cc[c:3]([C:2](F)(F)F)cc1",
                 anchor=a, spacer="para-phenylene", terminal="CF3", series="B"),
            dict(name=f"4-F-C6H4-{tag}", smiles=f"{asmi}c1cc[c:3]([F:2])cc1",
                 anchor=a, spacer="para-phenylene", terminal="F", series="B"),
            dict(name=f"3,5-(CF3)2-C6H3-{tag}", smiles=f"{asmi}c1cc(C(F)(F)F)[cH:3]c(C(F)(F)F)c1",
                 anchor=a, spacer="para-phenylene", terminal="H(3,5-CF3)", series="B"),
        ]
    return rows


def series_c():
    rows = []
    for anion, asmi, tag in [("SO3-", "[S:2](=O)(=O)[O-]", "SO3"), ("COO-", "[C:2](=O)[O-]", "COO")]:
        for n in (2, 3):
            for m in range(1, 6):
                smi = "[SH:1]" + "C" * n + "[N+:4](C)(C)" + "C" * (m - 1) + "[CH2:3]" + asmi
                rows.append(dict(name=f"ZW-{tag}-n{n}-m{m}", smiles=smi, anchor="SH",
                                 spacer=f"(CH2){n}-N+(CH3)2-(CH2){m}", terminal=anion,
                                 series="C", n=n, m=m))
    return rows


def references():
    return [
        dict(name="L-cysteine", smiles="N[C@@H](C[SH:1])[C:2](=O)O", anchor="SH",
             spacer="CH2-CH(NH2)", terminal="COOH", series="ref"),
        dict(name="L-cysteine-zwitterion", smiles="[NH3+:4][C@@H](C[SH:1])[C:2](=O)[O-]", anchor="SH",
             spacer="CH2-CH(NH3+)", terminal="COO-", series="ref"),
        dict(name="3-mercaptopropionic-acid", smiles="[SH:1]CC[C:2](=O)O", anchor="SH",
             spacer="(CH2)2", terminal="COOH", series="ref"),
        dict(name="thiophenol", smiles="[SH:1]c1cc[cH:3]cc1", anchor="SH",
             spacer="para-phenylene", terminal="H", series="ref"),
    ]


# 실험 수직/단열 이온화 에너지 (NIST WebBook, 기체상, eV)
VALIDATION = [
    dict(name="benzene", smiles="c1ccccc1", exp_ip=9.24),
    dict(name="toluene", smiles="Cc1ccccc1", exp_ip=8.83),
    dict(name="anisole", smiles="COc1ccccc1", exp_ip=8.20),
    dict(name="phenol", smiles="Oc1ccccc1", exp_ip=8.49),
    dict(name="aniline", smiles="Nc1ccccc1", exp_ip=7.72),
    dict(name="thiophenol-gas", smiles="[SH:1]c1cc[cH:3]cc1", exp_ip=8.30),
]


def validation():
    return [dict(name=v["name"], smiles=v["smiles"], anchor="-", spacer="-", terminal="-",
                 series="val", exp_ip=v["exp_ip"]) for v in VALIDATION]


def parse_entry(row):
    """SMILES 파싱·정규화. 성공하면 (mol_with_H, info) 반환, 실패하면 예외."""
    mol = Chem.MolFromSmiles(row["smiles"])
    if mol is None:
        raise ValueError("RDKit 파싱 실패")
    Chem.SanitizeMol(mol)
    if Descriptors.NumRadicalElectrons(mol):
        raise ValueError("라디칼 전자가 있음 (SMILES 조립 오류 가능)")
    molH = Chem.AddHs(mol)
    idx = {}
    for atom in molH.GetAtoms():
        m = atom.GetAtomMapNum()
        if m:
            idx[m] = atom.GetIdx()
    anchor = idx.get(1)
    terminal = idx.get(2)
    if terminal is None and 3 in idx:  # 말단 -H: map-3 원자에 붙은 수소
        hs = [n.GetIdx() for n in molH.GetAtomWithIdx(idx[3]).GetNeighbors() if n.GetAtomicNum() == 1]
        terminal = hs[0] if hs else None
    nplus = idx.get(4)
    for atom in molH.GetAtoms():
        atom.SetAtomMapNum(0)
    canon = Chem.MolToSmiles(Chem.RemoveHs(molH))
    info = dict(anchor_idx=anchor, terminal_idx=terminal, nplus_idx=nplus,
                canonical_smiles=canon, charge=Chem.GetFormalCharge(molH))
    return molH, info


def build(path="library.csv"):
    rows = series_a() + series_b() + series_c() + references() + validation()
    failed, seen = [], {}
    out = []
    for r in rows:
        try:
            molH, info = parse_entry(r)
        except Exception as e:  # noqa
            failed.append((r["name"], r["smiles"], str(e)))
            continue
        r = dict(r)
        r.update(info)
        key = (info["canonical_smiles"], r["series"] == "val")
        r["duplicate_of"] = seen.get(key, "")
        if key not in seen:
            seen[key] = r["name"]
        out.append(r)
    cols = ["name", "smiles", "anchor", "spacer", "terminal", "series", "n", "m",
            "canonical_smiles", "charge", "anchor_idx", "terminal_idx", "nplus_idx",
            "duplicate_of", "exp_ip"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in out:
            w.writerow({c: r.get(c, "") for c in cols})
    return out, failed


if __name__ == "__main__":
    out, failed = build()
    uniq = [r for r in out if not r["duplicate_of"]]
    print(f"총 {len(out)}개 항목, 고유 분자 {len(uniq)}개, 파싱 실패 {len(failed)}개")
    for s in ["A", "B", "C", "ref", "val"]:
        print(f"  시리즈 {s}: {sum(r['series']==s for r in out)}개 (고유 {sum(r['series']==s and not r['duplicate_of'] for r in out)})")
    dups = [(r["name"], r["duplicate_of"]) for r in out if r["duplicate_of"]]
    if dups:
        print("중복(같은 정규 SMILES):")
        for d in dups:
            print(f"  {d[0]} == {d[1]}")
    if failed:
        print("파싱 실패:")
        for f in failed:
            print("  ", f)
