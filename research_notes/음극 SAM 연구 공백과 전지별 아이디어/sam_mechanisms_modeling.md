# Cross-cutting science of SAMs at battery anode interfaces: stability, mechanism, characterization, modeling (plus graphite / hard-carbon / current-collector SAMs)

Scope note for the report writer: this file covers the chemistry-agnostic fundamentals. Per-chemistry molecule lists (Li, Si, Zn, Mg/Ca/Al/Na/K) are left to other researchers. Emphasis 2018–2026.
Method/limitation note: ~45 web searches were completed; the per-turn search budget then ran out. Direct full-text fetch was blocked for every publisher host tried (ScienceDirect, ACS, Wiley, Nature, Science, Springer, RSC, PMC/PubMed, arXiv, ChemRxiv, Semantic Scholar API, NSF-PAR). All findings below are therefore abstract/snippet-level, each with a link; numbers that I could not confirm at full-text level are marked "(abstract-level, unverified in full text)". Citations are given as completely as the snippets allowed; where author lists or DOIs are missing they are flagged.

---

## Q1. Is there consensus on whether SAMs survive at ~0 V vs Li+/Li (or Na+/Na, Mg2+/Mg, −0.76 V vs SHE for Zn)? What experiments show?

### Takeaway
There is no consensus and, more importantly, almost no direct measurement: the classical electrochemistry literature shows thiolate SAMs reductively desorb (1-electron, S–metal bond cleavage) at potentials far positive of 0 V vs Li+/Li (Cu is the most robust common substrate, ~0.6 V more negative than Au), and the thiol/Li-battery literature puts thiol reduction at 1.2–1.6 V vs Li+/Li; yet most battery SAM papers implicitly assume the monolayer persists, while a minority ("electrochemically active" SAMs) explicitly design it to decompose into a LiF-rich SEI. Operando/in situ tracking of SAM fate during plating was not found anywhere.

### Cited Findings
**Reductive desorption of thiolates (classical electrochemistry, aqueous/alkaline conditions)**
- Alkanethiolate SAMs on polycrystalline Cu and Cu(111) electrodesorb at potentials ~0.6 V more negative than the same SAM on Au; DFT on methanethiolate/Cu(111) attributes the shift to (i) the energy to add an electron to the thiolate–metal system and (ii) the energy to break the S–metal bond, with adsorbate–adsorbate interactions nearly substrate-independent — [Electrodesorption Potentials of Self-Assembled Alkanethiolate Monolayers on Copper Electrodes, J. Phys. Chem. B 2003, 107, 13446 (authors per snippet not shown; believed Azzaroni/Vela/Salvarezza group — unverified)](https://pubs.acs.org/jpcbfk/article-abstract/107/48/13446/3489121/Electrodesorption-Potentials-of-Self-Assembled?redirectedFrom=fulltext); for short-chain thiols, thiolates on Cu resist reductive desorption better than on Au, Ag and Ni — [ResearchGate record of same paper](https://www.researchgate.net/publication/235934359_Electrodesorption_Potentials_of_Self-Assembled_Alkanethiolate_Monolayers_on_Copper_Electrodes_An_Experimental_and_Theoretical_Study)
- Reductive stability order Au < Pt < Cu for thiolate SAMs; desorption is a one-electron transfer; peak potential shifts negative with longer alkyl chains; a full monolayer on Au reportedly holds until ~0.33 V vs RHE (abstract-level, unverified in full text) — [Electrochemical Stability of Thiolate SAMs on Au, Pt, Cu (ChemRxiv 2022-dt7dk)](https://chemrxiv.org/doi/pdf/10.26434/chemrxiv-2022-dt7dk); [NSF-PAR copy](https://par.nsf.gov/servlets/purl/10433337); related thesis [Electrochemical Reductive Stability of SAMs (U. Michigan Deep Blue)](https://deepblue.lib.umich.edu/bitstream/handle/2027.42/174390/rebonsa_1.pdf?sequence=1)
- Reductive-desorption peak moves negative with faster sweep rate, indicating kinetically slow desorption — [Hokkaido Univ. repository, J. Phys. Chem. B 108, 6422 (2004)](https://eprints.lib.hokudai.ac.jp/repo/huscap/all/50227/JPCB108-20_6422-6428.pdf)
- On reactive metals, chalcogenide–C bond scission (not just S–metal cleavage) can occur upon thiol adsorption — [Review: Self-assembled monolayers of thiolates on metals (sulfur–metal chemistry)](https://www.academia.edu/30485024/Self_assembled_monolayers_of_thiolates_on_metals_a_review_article_on_sulfur_metal_chemistry_and_surface_structures)
- General thiol electrochemistry in Li batteries: thiols show reductive decomposition at 1.2–1.6 V vs Li/Li+, i.e., far above the plating potential — [Thiols and Their Derivatives: Electrochemical Behaviors and Applications in Lithium Batteries (IntechOpen, online-first)](https://www.intechopen.com/online-first/1263212)

**What battery SAM papers actually report about fate**
- "Electrochemically active" monolayer (1,3-benzenedisulfonyl fluoride, BDSF) on Cu is designed to decompose on charging to form a stable interfacial layer; SEI has a LiF-rich inner layer and amorphous outer layer — [Gao et al., Nature Energy 5, 534–542 (2020), DOI 10.1038/s41560-020-0640-7 (OSTI record)](https://www.osti.gov/pages/biblio/1633448); [phys.org summary](https://phys.org/news/2020-08-thin-layer-battery-cold.html); summarized in [Surface engineering toward stable lithium metal anodes, Sci. Adv. 2023, DOI 10.1126/sciadv.adf1550](https://www.science.org/doi/10.1126/sciadv.adf1550)
- Si-containing monolayer on Li: "inner Si–O–Li components" persist and promote LiF formation in the inner SEI (abstract-level) — [J. Power Sources 2025, S0378775325013254](https://www.sciencedirect.com/science/article/abs/pii/S0378775325013254)
- Silane coupling agents anchor to Li via Li–O–Si bonds formed from methoxy + LiOH, with organic groups participating in SEI — [Wang et al., Adv. Mater. 2021, DOI 10.1002/adma.202008133](https://advanced.onlinelibrary.wiley.com/doi/10.1002/adma.202008133); [Recent progress in silane coupling agents for ASSLBs, Chem. Sci. 2026, DOI 10.1039/d6sc04626f](https://pubs.rsc.org/sc/article/doi/10.1039/d6sc04626f/1291634/Recent-progress-in-silane-coupling-agents-for-high)
- Phosphonic-acid SAMs on Li bind via lithium phosphate; long chains block O2/H2O — [Nano Research 2019, DOI 10.1007/s12274-019-2565-7](https://link.springer.com/article/10.1007/s12274-019-2565-7)
- Thiol SAM on Cu: authors explicitly claim the monolayer enables electron tunneling during plating/stripping (i.e., assume persistence) — [Electrochim. Acta 2024, S0013468624002421](https://www.sciencedirect.com/science/article/abs/pii/S0013468624002421)

**Characterization tools that exist but have not been applied to SAM fate**
- Operando XPS has been run on sulfide SEs cycled to −0.05 V vs Li/Li+ (technique exists) — [Operando XPS + pressure, PMC12632433](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12632433/)
- ToF-SIMS sputter depth profiling protocol for thin interphases/coatings on Li metal (Ar+ at 500 eV recommended for thin layers) — [ToF-SIMS depth profiling on Li metal, Commun. Chem. 2025, s42004-025-01426-0](https://www.nature.com/articles/s42004-025-01426-0); in situ ToF-SIMS of Li|SE interfaces — [Otto et al., Adv. Mater. Interfaces 2022, DOI 10.1002/admi.202102387](https://advanced.onlinelibrary.wiley.com/doi/10.1002/admi.202102387)
- Static-SIMS caution for fragile monolayers: OPA on mica — static limit questioned above ~1×10^11 ions cm−2 (as reported in snippet; verify) — [ResearchGate 237846551](https://www.researchgate.net/publication/237846551_Primary_ion_fluence_dependence_in_time-of-flight_SIMS_of_a_self-assembled_monolayer_of_octadecylphosphonic_acid_molecules_on_mica_discussion_of_static_limit); OPA SAM ToF-SIMS/XPS correlation — [Langmuir 22, 9244](https://pubs.acs.org/langd5/article-abstract/22/22/9244/3593088/ToF-SIMS-Investigation-of-Octadecylphosphonic-Acid?redirectedFrom=fulltext)
- Any spatially resolved SEI analysis is a "snapshot", especially on Li metal with massive volume change — [XPS and SIMS of SEI on Li in ether electrolytes](https://www.researchgate.net/publication/321701536_XPS_and_SIMS_Analysis_of_Solid_Electrolyte_Interphases_on_Lithium_Formed_by_Ether-Based_Electrolytes)
- In situ AFM/SECM and operando Raman for Li deposition exist but no SAM application found — [PNAS 2023 in situ/operando Li review](https://www.pnas.org/doi/10.1073/pnas.2220419120); [Operando Raman at microelectrode, PMC12243104](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12243104/)

### Inferences
- Combining the classical data (thiolate desorption on Cu ≈ 0.6 V negative of Au; thiol reduction 1.2–1.6 V vs Li+/Li) with the plating potential (0 V vs Li+/Li ≈ −3.04 V vs SHE) implies an intact thiolate SAM at 0 V vs Li+/Li is thermodynamically implausible; survival, if any, must be kinetic (slow desorption, buried under SEI, or anchored by a converted species such as Li2S/Cu2S/LiF). The Zn case (−0.76 V vs SHE, aqueous) is far milder and partial survival is plausible, consistent with the many Zn SAM papers reporting thousands of hours of cycling.
- The field is bifurcated: "passivating/persistent" SAM claims (thiols, phosphonates, silanes) vs. "sacrificial/templating" SAM claims (BDSF, fluorinated silanes). The same anchoring chemistry is sometimes claimed to do both, without discriminating experiments.
- The Si–O–Li and Li–O–P anchors are more likely than Cu–S to survive as inorganic residues (Li silicate/phosphate), but no study quantifies cleavage vs. retention.

### Gaps
- No paper found that tracks a SAM's anchoring bond (S 2p, P 2p, Si 2p) by XPS/ToF-SIMS as a function of potential down to 0 V vs Li+/Li, nor any operando (EQCM, SERS, in situ AFM) observation of SAM loss during plating.
- No reductive-desorption data for thiolates on Cu in nonaqueous Li electrolytes (all classical data are aqueous/alkaline).
- No data at all for SAM stability on Na, K, Mg, or Ca metal surfaces.
- Full text of the ChemRxiv Au/Pt/Cu stability study and of the J. Phys. Chem. B 2003 Cu paper could not be accessed; the exact potentials and author lists should be confirmed.

---

## Q2. How do SAMs template SEI formation (SAM-to-SEI conversion, LiF-rich SEI from fluorinated SAMs, molecular dipole / interfacial field / desolvation effects)?

### Takeaway
Two templating routes are documented: (a) sacrificial conversion, where a reactive terminal group (sulfonyl fluoride, C–F, Si–OR) is reduced into LiF/Li-silicate-rich inner SEI (Gao Nat. Energy 2020; F-terminated silane on graphite 2024); and (b) anion/solvent steering, where a persistent terminal group (carboxylate, tetrazole, ionic dipole) changes local solvation and drives anion (TFSI−, NO3−) decomposition into an inorganic SEI (Liu Science 2022; tetrazole-Cu 2025). Dipole effects are invoked widely but quantified rarely; work-function shifts of 0.3–2 eV by SAMs are well established outside batteries but have not been measured on battery current collectors.

### Cited Findings
**Sacrificial conversion**
- BDSF monolayer on Cu decomposes on charging → LiF-rich inner SEI, enabling low-temperature and high-rate Li metal cycling — [Gao et al., Nature Energy 2020 (OSTI)](https://www.osti.gov/pages/biblio/1633448); mechanism framed as benzenesulfinate → lithiophilic sites, fluoride → in-situ LiF nuclei — [Fluorine Chemistry in Rechargeable Batteries (UTS repository PDF)](https://opus.lib.uts.edu.au/bitstream/10453/177440/3/Fluorine%20Chemistry%20in%20Rechargeable%20Batteries.pdf); associated patent describes BDSF bonded to Cu "via thiol-Cu reaction" and AIMD of PC:FEC/Li(100) with LiF nanoparticle — [US Patent 12,658,447](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12658447)
- Fluorine-terminated silane SAM grafted on spherical graphite: F groups decompose during SEI formation to give a LiF-dominated inorganic SEI layer, improving fast-charging LCO//graphite full cells — [ACS Appl. Mater. Interfaces 2024, DOI 10.1021/acsami.3c15639](https://pubs.acs.org/doi/10.1021/acsami.3c15639); [PubMed 38272467](https://pubmed.ncbi.nlm.nih.gov/38272467)
- Si-containing monolayer on Li: Si–O–Li inner components promote LiF formation within inner SEI (quasi-solid-state cells) — [J. Power Sources 2025](https://www.sciencedirect.com/science/article/abs/pii/S0378775325013254)
- Review framing: modified current collectors with active F/O species raise SEI LiF content — [Facile heat treatment of Cu for anode-free (PMC9866673), review section](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9866673/)

**Anion/solvent steering by persistent terminal groups**
- Carboxyl-terminated SAMs (on separator side) promote a LiF-rich SEI with enhanced Li+ transport in Li metal cells — [Self-assembled monolayers direct a LiF-rich interphase toward long-life lithium metal batteries, Science 375, 739–745 (2022), DOI 10.1126/science.abn1818](https://www.science.org/doi/10.1126/science.abn1818) (first author Liu per citing sources; full author list unverified)
- 5-methyl-1H-tetrazole SAM on Cu anchored via Cu–N bonds promotes preferential LiNO3 decomposition into a stable SEI; claimed 600 Wh kg−1 anode-free cells (abstract-level) — [Molecular Interfacial Engineering of Lithiophilic Protection Layer, PMC13432276](https://pmc.ncbi.nlm.nih.gov/articles/PMC13432276/)
- HMDS SAM on Cu foil mitigates electrolyte decomposition and promotes uniform Li deposition — [Electrochim. Acta 2024, S0013468624013914](https://www.sciencedirect.com/science/article/abs/pii/S0013468624013914)
- Review of SAM-inspired artificial SEI design — [Self-Assembly Monolayer Inspired Stable Artificial SEI Design, 2023 (PubMed 37079652)](https://pubmed.ncbi.nlm.nih.gov/37079652/)

**Dipole / interfacial field / desolvation**
- SAM on graphite lowers interfacial Li+ desolvation barrier "from 9.52 to 6.65 eV" (as quoted in snippet; units implausible for a single-ion barrier — likely kJ/mol or a computed cluster energy; verify) for low-temperature LIBs — [Sci. China Chem. 2026, DOI 10.1007/s11426-026-3770-8](https://link.springer.com/article/10.1007/s11426-026-3770-8)
- Dipolar-molecule-grafted separator: ionic conductivity 0.517 vs 0.308 mS cm−1 (pristine), tLi+ = 0.646, dendrite-free Li — [Dipole-mediated interfacial solvation, Chem. Sci. 2026, 17, 5145](https://pubs.rsc.org/sc/article/17/10/5145/911602/Dipole-mediated-interfacial-solvation-for)
- Interfacial polarity of dipolar SAMs used as an explicit descriptor of electrode–electrolyte interaction (cathode side, stable to 4.7 V) — [Nature Nanotechnology 2026, s41565-026-02152-x](https://www.nature.com/articles/s41565-026-02152-x)
- Interfacial layers with desolvation function for Li metal — [Nano Lett. 2024, 24(26), 8055](https://pubs.acs.org/nalefd/article-abstract/24/26/8055/148852/Interfacial-Layers-with-Desolvation-Function?redirectedFrom=fulltext)
- Field-aligned amphiphilic anions (perfluoroalkyl tail + polar head) form ordered interfacial layers promoting desolvation and F/N-rich SEI — [Supramolecular Chemistry for LMBs, CCS Chem. 2025](https://www.chinesechemsoc.org/doi/10.31635/ccschem.025.202506033)
- Work-function tuning by SAMs (outside batteries): on CuO (Φ ≈ 4.9 eV), fluorinated and non-fluorinated thiol SAMs shift Φ via molecular dipoles; Kelvin probe shifts ≈ +0.28 to +0.5 eV — [Tuning of Metal Work Functions with SAMs (ResearchGate 30475453; believed de Boer et al., Adv. Mater. 2005 — unverified)](https://www.researchgate.net/publication/30475453_Tuning_of_Metal_Work_Functions_with_Self-Assembled_Monolayers); fluorinated alkanethiolates on Ag(111) raise Φ by up to 2 eV, non-fluorinated lower it — [arXiv 0705.0651](https://arxiv.org/pdf/0705.0651); KPFM of SAM surface potentials on Au — [ResearchGate 231681464](https://www.researchgate.net/publication/231681464_Kelvin_Probe_Force_Microscopy_on_Surfaces_Investigation_of_the_Surface_Potential_of_Self-Assembled_Monolayers_on_Gold); chemical vs electrostatic XPS shifts in SAMs — [J. Phys. Chem. C 2016, acs.jpcc.5b12387](https://pubs.acs.org/doi/10.1021/acs.jpcc.5b12387)
- The 2024 thiol-on-Cu Li paper attributes nucleation changes to lowered surface energy, not to a measured work-function/dipole — [Electrochim. Acta 2024](https://www.sciencedirect.com/science/article/abs/pii/S0013468624002421)

### Inferences
- The clearest mechanistic evidence (XPS-level) exists for the sacrificial route; the anion-steering route is supported by SEI composition but rarely by direct measurement of interfacial solvation.
- Dipole/work-function arguments in battery SAM papers are almost always DFT-only (electrostatic potential maps, ESP) rather than measured (KPFM, UPS); the toolset exists in the organic-electronics SAM literature and could be transplanted directly.
- A fluorinated terminal group does double duty (dipole + LiF source), which confounds attribution of improvements to "dipole" vs "LiF".

### Gaps
- No study found that separates dipole effect from conversion-chemistry effect (e.g., isostructural SAMs differing only in terminal dipole but both non-reactive).
- No measured work function / surface potential for SAM-modified Cu/Zn current collectors in a battery context.
- Full-text details (numbers, methods) of Science 2022 abn1818, Nature Energy 2020, Nano Lett. 2024 could not be accessed.

---

## Q3. Ion transport through SAMs: tunneling vs defect-mediated, chain-length dependence, head/terminal-group effects, measured Li+/Zn2+ permeation

### Takeaway
Direct measurements of Li+ or Zn2+ permeation through a battery SAM are absent; what exists is (i) the classical SAM-junction literature establishing exponential electron-tunneling decay with chain length (β ≈ 1.25 Å−1 for alkanes) and strong dependence on packing/order, (ii) a single 2024 Li/Cu paper invoking tunneling through a short thiol spacer (MPA, 17.6 mV nucleation overpotential, 98.4% CE), and (iii) Zn studies using MD to argue the SAM reshapes interfacial water/H-bond networks (proton blocking) and Zn2+ solvation rather than acting as a sieve.

### Cited Findings
- Tunneling current through alkane(di)thiol monolayers is temperature-independent and decays exponentially with molecular length — [PNAS 2007, 10.1073/pnas.0701472104](https://www.pnas.org/doi/pdf/10.1073/pnas.0701472104); decay constant ≈ 1.25 Å−1 (ab initio, consistent with Hg–Hg junction experiments) — [Long-range electron transfer through alkanethiol mono/bilayers (ResearchGate 10771577)](https://www.researchgate.net/publication/10771577_Long-Range_Electron_Transfer_through_Monolayers_and_Bilayers_of_Alkanethiols_in_Electrochemically_Controlled_Hg-Hg_Tunneling_Junctions); tunneling efficiency falls substantially when chains become disordered — [Electron tunneling across hexadecanethiolate on Hg, J. Phys. Chem. B 103, 8544](https://pubs.acs.org/jpcbfk/article/103/40/8544/1289708/Electron-Tunneling-Across-Hexadecanethiolate); impedance/potentiometric characterization of n-alkanethiol SAMs on Hg — [J. Electroanal. Chem. 2006](https://www.sciencedirect.com/science/article/abs/pii/S0022072806002427); simultaneous nanoindentation + tunneling — [ResearchGate 7076115](https://www.researchgate.net/publication/7076115_Simultaneous_Nanoindentation_and_Electron_Tunneling_through_Alkanethiol_Self-Assembled_Monolayers)
- Li/Cu: thiol SAMs on planar Cu; "alkyl spacer length and end groups tune lithiophilicity and electron/ion transfer"; best case 3-mercaptopropionic acid: nucleation overpotential 17.6 mV, ~98.4% CE over 250 cycles; authors state SAMs permit electron tunneling through the organic spacer unlike thick insulating polymer coatings. Chain length was not isolated as a single variable (per snippet) — [Electrochim. Acta 2024, S0013468624002421](https://www.sciencedirect.com/science/article/abs/pii/S0013468624002421)
- Zn: BMIM+ cation SAM on Ti3C2Tx; MD shows reduced interfacial H2O density, broken H-bond connectivity, blocked interfacial proton transport → HER-free; symmetric cells 3000 h at 1 mA cm−2 and 1000 h at 5 mA cm−2; reported water coordination numbers ~25–27 around SAM vs 4–6 around Zn2+ (snippet-level) — [Chen et al., Angew. Chem. Int. Ed. 2024, 63, e202404825](https://doi.org/10.1002/anie.202404825); [PubMed 38647332](https://pubmed.ncbi.nlm.nih.gov/38647332/)
- Zn: fluoroalkylsilane SAM on Cu-modified Zn repels water while "facilitating Zn2+ transport" and yields compact Zn(002) deposition — [Self-assembly of super-hydrophobic and zincophilic surface monolayer, 2025 (S240582972500279X)](https://www.sciencedirect.com/science/article/abs/pii/S240582972500279X)
- Zn: PFDT SAM creates zincophilic δ−F···Zn2+ sites via C–F bonds; fluorinated SAM claimed to enhance Zn2+ kinetics through ion–dipole interactions — [PFDT SAMs for Stable Zn Anodes, ACS Omega 2026, 11(2), 3569 (PMC12824933)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12824933)
- Zn: "interfacial confinement effect" of self-adsorbed monolayer — [Huo et al., Adv. Sci. 2025, DOI 10.1002/advs.202413731](https://advanced.onlinelibrary.wiley.com/doi/10.1002/advs.202413731); thiol SAM reduces interfacial water concentration — [Monolayer Thiol Engineered Covalent Interface (ResearchGate 365966815)](https://www.researchgate.net/publication/365966815_Monolayer_Thiol_Engineered_Covalent_Interface_toward_Stable_Zinc_Metal_Anode); SAM molecular engineering for high-DOD Zn — [eScience 2023, DOI 10.1016/j.esci.2023.100205](https://www.sciengine.com/doi/10.1016/j.esci.2023.100205); nucleophilic interfacial SAM layer — [ACS AMI 2025, acsami.5c03770](https://pubs.acs.org/doi/10.1021/acsami.5c03770); multifunctional monomolecular interface — [Small 2026, smll.75025](https://onlinelibrary.wiley.com/doi/10.1002/smll.75025)
- Thickness lesson from inorganic analogues: thicker evaporated LiF/Li3N on Cu raised overpotential and aggregated Li+ flux — artificial interphases must be thin and uniform — [Re-evaluating nano-sized inorganic protective layers on Cu, 2023](https://www.sciencedirect.com/science/article/abs/pii/S1001841723010355)

### Inferences
- For Li plating on a SAM-coated collector, electron transfer must either tunnel through ≤~1.5 nm of alkyl chain or pass through defects; β ≈ 1.25 Å−1 means a C18 chain (~2.3 nm) attenuates tunneling by >10^12 relative to bare metal, so plating on long-chain SAMs is essentially certain to be defect-/pinhole-mediated or to occur after SAM breakdown. Short chains (C2–C6, e.g., MPA) are the only regime where tunneling is a credible mechanism, which matches the 2024 Li/Cu result.
- Zn work has converged on an "interfacial water/solvation editing" picture rather than ion sieving; whether Zn2+ crosses the monolayer or deposits at defects is not addressed.
- Chain-length and packing-density dependence of CE is a natural systematic experiment that no one has reported.

### Gaps
- No measured Li+ or Zn2+ permeability/flux through a battery SAM (e.g., EIS with redox probe, SECM, EQCM, or isotopic tracers).
- No controlled chain-length series (C2…C18, same anchor/terminal) with CE/nucleation overpotential outputs for Li or Zn; the Zn systematic-comparison search returned none.
- Classical pinhole/defect-density electrochemistry (redox-probe blocking) has not been transferred to battery SAMs.

---

## Q4. SAM effects on nucleation: wetting, lithiophilicity/zincophilicity, nucleation overpotential, texture/epitaxy (Zn(002), Li(110))

### Takeaway
For Zn, a coherent and DFT-supported mechanism exists: thiol/fluorothiol SAMs adsorb preferentially on high-energy Zn(100)/(101) facets, shielding them so (002) basal growth dominates (~85% (002) texture reported in 2026); nucleation overpotential is typically raised by hydrophobic SAMs. For Li, SAM-based texture control is essentially unreported (one unverified abstract claims (110)-preferred Li via electrochemically active SAM), and the very existence of strong (110) texture in Li electrodeposits is contested; lithiophilicity is quantified only via DFT binding energies on non-SAM modifiers (bare Cu ≈ −0.26 to −1.6 eV depending on facet/reference).

### Cited Findings
**Zn texture / nucleation**
- PFDT− adsorption energy on Zn(100)/Zn(101) is significantly higher than on Zn(002) → preferential Zn(002) deposition — [PFDT SAMs, ACS Omega / PMC12824933](https://pmc.ncbi.nlm.nih.gov/articles/PMC12824933)
- 2,2′-dithiobis(ethylamine) SAM: thiols selectively assemble on (100)/(101) planes, shielding their nucleation → single (002)-plane growth, ~85% (002) texture, ultra-high areal capacity — [Wang et al., Adv. Funct. Mater. 2026, DOI 10.1002/adfm.202532042](https://advanced.onlinelibrary.wiley.com/doi/abs/10.1002/adfm.202532042)
- Fluoroalkylsilane/low-surface-energy SAM gives ultrafine nucleation and compact deposition along Zn(002) — [2025, S240582972500279X](https://www.sciencedirect.com/science/article/abs/pii/S240582972500279X)
- Non-SAM mechanistic context: higher overpotential at high current selectively accelerates (100) growth, exposing (002) and enabling epitaxy — [Zn(002) texture via corrosion, Chem. Eng. J. 2025](https://www.sciencedirect.com/science/article/abs/pii/S1385894725101137); moderately hydrophobic polymer layer instead enables Zn(100) deposition — [Nano-Micro Lett. 2026, s40820-026-02208-6](https://link.springer.com/article/10.1007/s40820-026-02208-6)
- SAM-coated Zn anodes generally show a larger nucleation overpotential than bare Zn, interpreted as denser/finer nucleation (abstract-level; specific molecule in snippet labeled "CHP", source not individually resolved) — Zn SAM result set incl. [ACS AMI 2025](https://pubs.acs.org/doi/10.1021/acsami.5c03770), [Adv. Sci. 2025](https://advanced.onlinelibrary.wiley.com/doi/10.1002/advs.202413731)
- Thiol monolayer on Zn: >3000 h at 1 mA cm−2, CE 99.91% over 2900 cycles (abstract-level) — [ResearchGate 365966815](https://www.researchgate.net/publication/365966815_Monolayer_Thiol_Engineered_Covalent_Interface_toward_Stable_Zinc_Metal_Anode)
- Reviews on Zn current-collector/interphase design — [Nano-Micro Lett. 2023, s40820-023-01177-4](https://link.springer.com/article/10.1007/s40820-023-01177-4); [Chem Bio Eng 2024, cbe.4c00053](https://pubs.acs.org/doi/10.1021/cbe.4c00053)

**Li nucleation / texture**
- Thiol SAMs on Cu lower surface energy → larger "dough-like" Li nuclei, denser deposits; MPA/Cu: 17.6 mV nucleation overpotential — [Electrochim. Acta 2024](https://www.sciencedirect.com/science/article/abs/pii/S0013468624002421)
- Li plated at ultrahigh capacity reported to be strongly (110)-oriented with faster kinetics — [Nat. Commun. 2021, s41467-021-26143-9](https://www.nature.com/articles/s41467-021-26143-9); critique argues XRD does not support dominant (110) texture (nearly random with slight (110)) — [Critical evaluation (PMC9525300)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9525300/); authors' reply — [Nat. Commun. 2022, s41467-022-32950-5](https://www.nature.com/articles/s41467-022-32950-5)
- Single-crystal Cu(110) induces lateral Li growth — [ResearchGate 351546662](https://www.researchgate.net/publication/351546662_Single_crystal_Cu_110_inducing_lateral_growth_of_electrodeposition_Li_for_dendrite-free_Li_metal-based_batteries); laser-constructed 3D Cu with orientation selectivity — [Energy Environ. Mater. 2024, eem2.12768](https://onlinelibrary.wiley.com/doi/10.1002/eem2.12768)
- Unverified abstract (ResearchGate-indexed, surfaced in search): SAM of electrochemically active molecules on current collectors directs (110)-preferred Li deposition even on (200)-oriented Li — no accessible full text; treat as unconfirmed.
- Lithiophilicity DFT benchmarks (non-SAM modifiers, useful reference values): Li binding on Cu planes −0.89/−0.71/−0.26 eV vs −4.19 eV on CuOHF — [Dual Lithiophilic Structure (OSTI 2311177)](https://www.osti.gov/pages/servlets/purl/2311177); Cu foil average Li binding 1.62 eV (PBE-D3(BJ)) vs HKUST-1 — [InfoMat 2025, inf2.70067](https://onlinelibrary.wiley.com/doi/full/10.1002/inf2.70067); CuSe(102) −1.50 vs Cu(111) −0.95 eV — [EEM 2023, eem2.12534](https://onlinelibrary.wiley.com/doi/full/10.1002/eem2.12534); VN 1.98 vs Cu 0.71 eV — [Small Sci. 2022, smsc.202100110](https://onlinelibrary.wiley.com/doi/full/10.1002/smsc.202100110); heteroatom-doped carbon lithiophilicity framework — [Sci. Adv. 2019, aau7728](https://www.science.org/doi/10.1126/sciadv.aau7728)
- Growing single-crystalline seeds on lithiophobic substrates (non-SAM, relevant nucleation physics) — [Nature Energy 2023, s41560-023-01202-1](https://www.nature.com/articles/s41560-023-01202-1); salt-philic/solvent-phobic coating — [Nature Energy 2023, s41560-023-01252-5](https://www.nature.com/articles/s41560-023-01252-5)

### Inferences
- Facet-selective adsorption is the one SAM–nucleation mechanism with a clear, falsifiable computational handle (ΔE_ads across facets) and it has only been exploited for Zn; the analogous calculation for Li(110)/(100)/(200) or Na facets on Cu with a SAM has not been reported.
- Reported Li binding energies for bare Cu span 0.26–1.6 eV across papers because of facet, reference state (atom vs bulk BCC Li), and functional; any SAM lithiophilicity screening must compute its own bare-metal baseline.
- Lower surface energy (hydrophobic SAM) raises Zn nucleation overpotential but was reported to lower Li nucleation overpotential (MPA case) — the "surface energy → overpotential" relation is not monotonic and depends on terminal-group ion affinity.

### Gaps
- No quantitative, cross-paper relation between a SAM descriptor (dipole, packing, chain length, terminal group) and nucleation overpotential or CE was found for Li; for Zn only facet-adsorption-energy ↔ (002) fraction is semi-quantitative.
- No SAM-based Li or Na texture control paper was verifiable.

---

## Q5. Computational approaches: DFT, AIMD, classical MD, MLIPs (MACE/CHGNet/M3GNet/UMA-OMat24), HT/active-learning screening, LLM/agent screening — which studies model SAM + alkali/alkaline-earth metal + substrate explicitly, and any with MLIPs?

### Takeaway
Explicit SAM + alkali-metal + substrate models are rare and DFT/AIMD-only (BDSF patent AIMD; ionic-liquid monolayer on Li AIMD 2024; facet-resolved DFT of thiolates on Zn); classical MD is widely used for Zn SAM/water interfaces. MLIPs (fine-tuned MACE-MP-0, Deep Potential, other MLFFs) have reached Li-metal/liquid-electrolyte SEI formation (2025–2026) but no MLIP study of an organic monolayer on a metal with Li/Na/Zn was found; MACE-OFF covers only organic elements and uMLIP benchmarks warn that low test error does not guarantee reliability at reconstructing battery interfaces. HT/ML screening exists for electrolyte additives (GNN over 75k molecules, Zn) and for F-containing interphase molecules (PubChem down-selection), and LLM agents have been applied to electrolyte molecule ranking/idea generation (2025–2026) — but none specifically to SAM molecules.

### Cited Findings
**Explicit SAM/monolayer + metal interface models (DFT/AIMD)**
- AIMD of PC:FEC electrolyte (283 atoms) on 10 layers Li(100), and of 380-atom electrolyte + 5 layers Li(100) + LiF nanoparticle to represent BDSF-derived interphase on Cu — [US Patent 12,658,447 (Stable metal anodes)](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12658447)
- First-principles/AIMD of an ionic-liquid monolayer on Li metal comparing two self-assembled configurations — [J. Phys. Chem. C 2024, DOI 10.1021/acs.jpcc.4c04828](https://pubs.acs.org/doi/10.1021/acs.jpcc.4c04828)
- Facet-resolved DFT adsorption of PFDT− on Zn(002)/(100)/(101) — [PMC12824933](https://pmc.ncbi.nlm.nih.gov/articles/PMC12824933); selective (100)/(101) shielding — [AFM 2026](https://advanced.onlinelibrary.wiley.com/doi/abs/10.1002/adfm.202532042)
- DFT on methanethiolate/Cu(111) desorption energetics (electron-addition + S–metal bond-breaking decomposition) — [J. Phys. Chem. B 2003](https://pubs.acs.org/jpcbfk/article-abstract/107/48/13446/3489121/Electrodesorption-Potentials-of-Self-Assembled?redirectedFrom=fulltext)
- Thiolate coverage DFT: max coverage 1/3 on (111), 1/2 on (100), 2/3 on nanoparticles — [Thiolate SAM review (academia.edu)](https://www.academia.edu/30485024/Self_assembled_monolayers_of_thiolates_on_metals_a_review_article_on_sulfur_metal_chemistry_and_surface_structures); DFT-fitted force field for SAM formation thermodynamics on Au — [ResearchGate 381396417](https://researchgate.net/publication/381396417_The_Thermodynamics_of_Self_Assembled_Monolayers_Formation_A_computational_and_Experimental_Study_of_Thiols_on_Flat_Gold_Surfaces)
- Non-SAM AIMD baselines: Li film on Cu, charge rearrangement localized at Li–Cu interface — [arXiv 2001.01236](https://arxiv.org/pdf/2001.01236); early SEI AIMD on graphite — [arXiv 1009.4154](https://arxiv.org/pdf/1009.4154); voltage-dependent interfacial processes — [arXiv 1312.2945](https://arxiv.org/pdf/1312.2945)

**Classical MD of SAM–electrolyte interfaces (Zn-dominated)**
- BMIM+ SAM/MXene: MD of interfacial water density, H-bond network, coordination numbers — [Angew. 2024, e202404825](https://doi.org/10.1002/anie.202404825)
- MD + DFT showing additive edits Zn2+ solvation sheath and excludes H2O from IHP; vermiculite channel MD — [Microstructures 2024](https://www.oaepublish.com/articles/microstructures.2024.114); [Chem. Eng. J. 2022, S1385894722034076](https://www.sciencedirect.com/science/article/abs/pii/S1385894722034076); zincophilic–hydrophobic interface — [PMC13065845](https://pmc.ncbi.nlm.nih.gov/articles/PMC13065845/)

**MLIPs at metal/electrolyte interfaces (no SAM yet)**
- Fine-tuned MACE-MP-0 for early SEI formation at Li metal/carbonate electrolyte, Li self-diffusion in metal/SEI/electrolyte, coatings of SEI components suppress reactions — [Towards atomistic simulations of the SEI in Li metal batteries, 2026, DOI 10.1186/s41313-026-00092-6](https://link.springer.com/article/10.1186/s41313-026-00092-6)
- Fine-tuned MACE-MP0: EC ring-opening near Li atoms/on Li metal; inaccurate functionals can overestimate rates by up to 9 orders of magnitude — [arXiv 2509.14067](https://arxiv.org/html/2509.14067)
- MLFF-MD with automated reaction analysis on ~4000-atom Li/LiFSI-EC interface reproduces Li2O near Li and LiF in outer SEI — [Univ. Tokyo press release 2026-10-28](https://www.t.u-tokyo.ac.jp/en/press/pr2026-10-28-001); MLFF-MD on 1296-atom Li(100) slab, anion effects, liquid Madelung potential — [Faraday Discuss. 2026, DOI 10.1039/d6fd00023a](https://pubs.rsc.org/fd/article/doi/10.1039/d6fd00023a/1282488/SEI-chemistry-enlightened-by-liquid-Madelung)
- Deep Potential MLMD of early SEI nucleation in concentrated LiTFSI/LiPF6 — [ResearchGate 400505672](https://www.researchgate.net/publication/400505672_Predictive_Machine_Learning_Molecular_Dynamics_of_SEI_Formation_in_Concentrated_LiTFSI_and_LiPF6_Electrolytes_for_Lithium_Metal_Batteries); Na SEI via interface reactor sampling — [arXiv 2604.07077](https://arxiv.org/pdf/2604.07077)
- "Effective" MLIPs for SEI materials; notes a major challenge applying MACE to SEI systems — [Mater. Horiz. 2025, 12, 10770](https://pubs.rsc.org/mh/article/12/24/10770/314371/Enabling-accurate-modelling-of-materials-for-a)
- FIRE: fine-tunes MACE-MPA-0 for solid–solid battery interfaces; argues low benchmark error ≠ trust at a reconstructing interface — [arXiv 2601.17847](https://arxiv.org/html/2601.17847); Li|β-Li3PS4 with ML-potential MD — [ORNL](https://impact.ornl.gov/en/publications/visualizing-the-sei-formation-between-lithium-metal-and-solid-sta/); Li/sulfide mechanochemistry — [arXiv 2609.01088](https://arxiv.org/pdf/2609.01088)
- Fine-tuned MLFF + Markov state models for Li transport: <2 meV atom−1, <30 meV Å−1 — [JCTC 2026, 22, 5373](https://pubs.acs.org/jctcce/article/22/11/5373/5185437/Bridging-Atomistic-and-Mesoscale-Lithium-Transport)
- uMLIP landscape: UMA trained on OMat24 + OC20 + ODAC25 + OMC25 + OMol25 — [Optimizing cross-domain transfer, Nat. Commun. 2026 / arXiv 2510.11241](https://www.nature.com/articles/s41467-026-70195-8); MACE trained with OMAT beats MACE-MP-0a on OC20 metal-surface benchmark — [arXiv 2510.25380](https://arxiv.org/html/2510.25380v1); MACE-OFF covers H,C,O,N,P,S,F,Cl,Br,I only (no metals) — [JACS 2025, 147, 17598](https://pubs.acs.org/jacsat/article/147/21/17598/3736678/MACE-OFF-Short-Range-Transferable-Machine-Learning); OMol25-trained potential for electrolyte solvation structure (no electrode) — [arXiv 2603.20183](https://arxiv.org/pdf/2603.20183); uMLIPs for metal defects (CHGNet weakest, OMat24-trained eqV2 strongest) — [Mach. Learn.: Sci. Technol. 2025, 10.1088/2632-2153/adea2d](https://iopscience.iop.org/article/10.1088/2632-2153/adea2d); MOFSimBench (metal–organic bonding benchmark) — [npj Comput. Mater. 2025](https://www.nature.com/articles/s41524-025-01872-3); UniFFBench — [arXiv 2508.05762](https://arxiv.org/pdf/2508.05762); oxide/water interface MLP precedent — [arXiv 2109.14068](https://arxiv.org/pdf/2109.14068); practical MLIP guide — [arXiv 2503.09814](https://arxiv.org/pdf/2503.09814)

**High-throughput / ML screening of interfacial molecules**
- F-containing molecules for inorganic–organic hybrid interphase on Li: PubChem high-throughput down-selection, DFT + ML predicting Li+ diffusion barrier from electronic descriptors — [Data-driven discovery of artificial hybrid interphase, Matter 2023 (S2590238523002989)](https://www.sciencedirect.com/science/article/pii/S2590238523002989); [review: ML for high-safety Li metal anode, PMC10820802](https://pmc.ncbi.nlm.nih.gov/articles/PMC10820802)
- GNN screened 75,024 organic molecules → 48 candidates; cyanoacetamide and hydantoin validated for Zn; EMolCurator combines HT screening with generative models — [Xu et al., Angew. 2025, DOI 10.1002/anie.202511389](https://onlinelibrary.wiley.com/doi/abs/10.1002/anie.202511389); molecular-descriptor ML → potassium L-aspartate, 99.86% avg CE over 5000 cycles — [Cai et al., Angew. 2026, DOI 10.1002/anie.202521382](https://onlinelibrary.wiley.com/doi/10.1002/anie.202521382); robotic HT platform: 180 small molecules, ML predicts CE — [Sci. Adv. 2025, adu4391](https://www.science.org/doi/10.1126/sciadv.adu4391)
- GNN screening of >10^9 solvents — [PubMed 42003686](https://pubmed.ncbi.nlm.nih.gov/42003686/); billion-scale motif hierarchy — [arXiv 2608.16364](https://arxiv.org/pdf/2608.16364); ~10,000-candidate HT-DFT SPE solvent screening with HOMO/LUMO, Li binding, DN descriptors — [Nano Lett. 2025, 25, 7801](https://pubs.acs.org/nalefd/article/25/19/7801/3752098/Machine-Learning-Guided-Screening-of-Advantageous)
- AI ecosystems for SSB interfaces (multistep HT screening of interface coatings) — [Sci. Adv. 2026, aea0638](https://www.science.org/doi/10.1126/sciadv.aea0638); ML paradigms for battery science — [npj Comput. Mater. 2025](https://www.nature.com/articles/s41524-025-01575-9); ML in electrolyte design review — [Chem 2026](https://www.cell.com/chem/abstract/S2451-9294(26)00090-2)

**LLM / agent-based molecule discovery**
- OmniScience: dual-agent setup grades ~700 candidate molecules (one LLM grades, GPT-o1 reviews) for LMB solvents/additives — [arXiv 2503.17604](https://arxiv.org/pdf/2503.17604)
- LLMB: AI agent for Li-metal battery research combining LLM text mining, ML prediction, molecular simulation — [ACS Cent. Sci. 2026, 12(4), 484](https://pubs.acs.org/acscii/article/12/4/484/5140991/LLMB-AI-Agent-for-Lithium-Metal-Battery-Research)
- Multi-agent idea generator for Zn deep-eutectic electrolytes — [Adv. Mater. 2025, DOI 10.1002/adma.202502649 (PMC12355564)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12355564/)
- LLM + simulation-tool agent, inverse design of electrolyte additives by MW and frontier-orbital targets — [J. Mater. Inf. 2025](https://www.oaepublish.com/articles/jmi.2025.69); multi-target generative additive design — [PMC12392859](https://pmc.ncbi.nlm.nih.gov/articles/PMC12392859); VAE-based additive generation — [Appl. Sci. 2025, 10.3390/app15073640](https://doi.org/10.3390/app15073640)

### Inferences
- The ingredients for a uMLIP SAM screening pipeline exist separately — OMol25 (organics incl. S/P/F), OMat24/OC20 (metal surfaces, adsorbates), UMA (joint) — but no one has validated a uMLIP on thiolate/phosphonate/silane–Cu/Zn/Li bonding, on reductive desorption barriers, or on Li adatom insertion into a monolayer. MACE-OFF is unusable for metals; MACE-MP-0 needs fine-tuning even for liquid electrolytes on Li.
- Existing Li-SEI MLIP work (2025–2026) demonstrates the needed system sizes (1300–4000 atoms) and reaction-analysis tooling; adding a SAM is a tractable extension, not a conceptual leap.
- Screening literature treats molecules as solvated additives (HOMO/LUMO, binding energy, DN); SAM-specific descriptors (anchor bond strength, packing/tilt, dipole normal to surface, desorption potential) have not been used as targets.

### Gaps
- Zero MLIP studies of an organic SAM on a battery metal; zero uMLIP benchmarks for thiolate/phosphonate/silane–metal bonds under reduction.
- No active-learning or Bayesian-optimization campaign over SAM chemical space for any battery anode (search for perovskite-SAM HT precedents was not completed due to budget).
- No LLM/agent workflow targeted at SAM or anchored-monolayer molecules; all agent work is on electrolyte solvents/additives.
- Mg/Ca/Na SAM computational studies: none found.

---

## Q6. SAMs on graphite, hard carbon, and Cu/Al foil current collectors for conventional Li-ion/Na-ion

### Takeaway
Graphite-side SAM work is small but real: a fluorinated silane SAM on graphite (2024) converts to a LiF-dominated SEI for fast charging, and a 2026 SAM on graphite lowers Li+ desolvation energy for low-temperature operation. No SAM study on hard carbon (Na/K) was found. For Al/Cu current collectors, SAM corrosion-inhibition and adhesion-promotion data exist only in non-battery corrosion/adhesion contexts (phosphonic acids, OTS, silane coupling agents); battery-cell validation of SAM-treated foils was not found.

### Cited Findings
- F-terminated silane SAM grafted on spherical graphite → LiF-dominated SEI inorganic layer → improved fast-charging/cycling in LCO//graphite full cells — [ACS AMI 2024, acsami.3c15639](https://pubs.acs.org/doi/10.1021/acsami.3c15639)
- SAM on graphite regulates interfacial solvation; desolvation barrier lowered (reported "9.52 → 6.65 eV"; units to verify) for low-temperature LIBs — [Sci. China Chem. 2026](https://link.springer.com/article/10.1007/s11426-026-3770-8)
- Fluoroalkylsilane SAM (<2 nm) on LNMO cathode (cathode analogue, shows silane SAMs are chemically compatible with carbonate electrolytes) — [PMC4995502](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4995502/)
- Hard carbon SEI reviews (no SAM): inorganic-rich/NaF-rich homogeneous SEI critical — [Frontiers Chem. 2022](https://www.frontiersin.org/journals/chemistry/articles/10.3389/fchem.2022.986541/pdf); [Small 2025, smll.202500645](https://onlinelibrary.wiley.com/doi/full/10.1002/smll.202500645); [Chem. Sci. 2024, 15, 6244](https://pubs.rsc.org/sc/article/15/17/6244/868826/Hard-carbon-for-sodium-ion-batteries-progress); [PMC9572906](https://pmc.ncbi.nlm.nih.gov/articles/PMC9572906/)
- Al current collector corrosion/passivation reviews (SAMs not covered; MXene self-assembled layer on Al improved cycling/rate) — [npj Mater. Degrad. 2024, s41529-024-00453-x](https://www.nature.com/articles/s41529-024-00453-x); [The Innovation Materials 2023](https://www.the-innovation.org/data/article/materials/preview/pdf/XINNMATERIALS-2023-0008.pdf); [J. Energy Storage 2021](https://www.sciencedirect.com/science/article/abs/pii/S2352152X21009233); primed Al collectors — [Electrochem. Commun. 2021](https://www.sciencedirect.com/science/article/pii/S1388248121000928); carbon coating on Al removes contact-resistance-dominated impedance — [SN Appl. Sci. 2022](https://link.springer.com/article/10.1007/s42452-022-05103-y)
- Non-battery SAM corrosion/adhesion data: phosphonic-acid SAMs on Cu cut icorr ~50× in NaCl, WCA >160° — [Copper Protection by Phosphonic Acid SAMs](https://www.researchgate.net/publication/216826536_Copper_Protection_by_Phosphonic_Acid_Self-assembled_Monolayers); OTS SAM on Cu: Ecorr +1.02 V, icorr 0.59 µA cm−2 in HCl — [PMC5745662](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5745662/); phosphonic acid/silane bonding to oxidized Al by XPS, fluorinated gives lowest adhesion/friction — [Springer chapter](https://link.springer.com/chapter/10.1007/978-3-540-74083-4_10); silane coupling monolayers for Al–polymer adhesion — [PMC11434525](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11434525/)
- SEI formation and Li plating on bare Cu in anode-free cells (baseline for collector SAMs) — [PMC8392351](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8392351/); anode-free progress review — [Energy Mater. Adv. 2024, 0168](https://spj.science.org/doi/10.34133/energymatadv.0168)

### Inferences
- For intercalation anodes the SAM is applied to the active particle, not the collector, and functions as a sacrificial SEI precursor (F→LiF) or a solvation modifier; this is mechanistically identical to the Cu-collector case, so stability questions from Q1 carry over (graphite lithiation at ~0.1 V vs Li+/Li).
- Low-potential cathodic protection of Cu and the native-oxide barrier on Al mean the collector-SAM value proposition in Li-ion is adhesion/corrosion during slurry processing and high-voltage Al corrosion, not SEI control; this is an under-explored transfer from the corrosion literature.

### Gaps
- No hard-carbon or Na-ion/K-ion SAM paper found (confirmed negative across two search phrasings).
- No cell-level data for phosphonate/silane-treated Al or Cu foils in Li-ion (impedance, adhesion, corrosion at >4.3 V, slurry wetting).
- Full texts of the two graphite SAM papers were not accessible; desolvation-energy units unverified.

---

## Q7. Which SAM descriptors (dipole, packing, chain length, terminal group) have been linked quantitatively to CE or nucleation overpotential? Is there a standardized SAM quality control in battery papers?

### Takeaway
No quantitative structure–property dataset exists: individual papers vary one molecule at a time and report CE/overpotential for one or two SAMs; packing density, tilt, and thickness are almost never reported in battery SAM papers, whereas the classical Langmuir literature established a multi-technique standard (ellipsometry + XPS + contact angle + IRRAS/AFM) decades ago.

### Cited Findings
- Only paired numbers found for Li: MPA/Cu, 17.6 mV nucleation overpotential, ~98.4% CE (250 cycles) — [Electrochim. Acta 2024](https://www.sciencedirect.com/science/article/abs/pii/S0013468624002421)
- Zn: thiol SAM 99.91% CE over 2900 cycles — [ResearchGate 365966815](https://www.researchgate.net/publication/365966815_Monolayer_Thiol_Engineered_Covalent_Interface_toward_Stable_Zinc_Metal_Anode); PFDT: facet ΔE_ads qualitative — [PMC12824933](https://pmc.ncbi.nlm.nih.gov/articles/PMC12824933); ~85% (002) with dithiobis(ethylamine) — [AFM 2026](https://advanced.onlinelibrary.wiley.com/doi/abs/10.1002/adfm.202532042)
- Non-SAM additive comparison hinting at chain-length effect: benzyl-group additive vs dodecyl surfactant (higher polarization, CE <99%) — [Microstructures 2024](https://www.oaepublish.com/articles/microstructures.2024.114)
- Classical SAM QC standard: ellipsometry, AFM, WCA and XPS routinely combined; XPS gives packing density and anchor state, ellipsometry gives thickness/formation kinetics — [SAMs as platforms for nanobiotechnology, ACS Appl. Nano Mater. 2025, 8(17), 8570](https://pubs.acs.org/aanmf6/article/8/17/8570/3646906/Self-Assembled-Monolayers-as-Platforms-for); thickness ~20 Å ellipsometry + XPS packing + PM-IRRAS + contact angle to rank film order — [Langmuir 2000, 16, 2220 (Lee group)](http://lee.chem.uh.edu/2000/Langv16p2220.pdf); contact angles reproducible to ~1° — [Langmuir 1998, 14, 5821](http://lee.chem.uh.edu/1998/Langv14p5821.pdf); XPS-based direct packing-density calculation for chelating thiols — [Langmuir 2005, la0475573](https://pubs.acs.org/doi/10.1021/la0475573); perfluoroalkylsilane SAM quality optimization via WCA, ellipsometry, XPS, AFM — [Langmuir 2012, 28, 11790](https://pubs.acs.org/langd5/article-abstract/28/32/11790/1670414/Optimizing-the-Quality-of-Monoreactive?redirectedFrom=fulltext); SAMs on Ru, dynamic contact angle sensitive to heterogeneity/disorder — [Langmuir 2017, acs.langmuir.7b01068](https://pubs.acs.org/doi/10.1021/acs.langmuir.7b01068)
- Work-function/dipole descriptors measured by Kelvin probe in organic electronics (0.28–0.5 eV on CuO; up to 2 eV on Ag) — [ResearchGate 30475453](https://www.researchgate.net/publication/30475453_Tuning_of_Metal_Work_Functions_with_Self-Assembled_Monolayers); [arXiv 0705.0651](https://arxiv.org/pdf/0705.0651)
- Reductive desorption charge (CV) is itself a classical coverage/packing metric and shifts with chain length — [ChemRxiv 2022](https://chemrxiv.org/doi/pdf/10.26434/chemrxiv-2022-dt7dk)

### Inferences
- Battery SAM papers overwhelmingly use contact angle (± XPS survey) as the only SAM QC; ellipsometric thickness, XPS-derived coverage, IRRAS order (CH2 stretch position), and reductive-desorption charge are rarely reported, so inter-paper comparisons of "a SAM" conflate coverage differences with chemistry differences.
- A minimal standard could be: ellipsometric thickness vs theoretical length (tilt), XPS anchor/substrate ratio (coverage), νas(CH2) position (order), advancing/receding WCA (hysteresis = defects), and reductive-desorption charge on the actual substrate.

### Gaps
- No review or dataset found that tabulates SAM descriptors vs electrochemical outcomes for any anode; a dedicated search for a SAM-QC best-practice review in batteries returned none.
- Could not verify which QC techniques each battery SAM paper used (full text blocked).

---

## Q8. What do reviews identify as the main open questions for molecular monolayers at anodes? (plus consolidated gap list for this topic)

### Takeaway
Dedicated SAM-in-battery reviews exist (Small 2024; Sci. Adv. 2023 surface-engineering review; 2023 SAM-inspired artificial SEI review) but their open-question sections could not be read (access blocked); the gap list below is therefore synthesized from the primary findings above and is explicitly labeled as inference.

### Cited Findings
- Review: Self-Assembled Molecular Layers as Interfacial Engineering Nanomaterials in Rechargeable Battery Applications — [He et al., Small 2024, DOI 10.1002/smll.202403537](https://onlinelibrary.wiley.com/doi/10.1002/smll.202403537) (attributes BDSF work to Gao et al. and classes it as "electrochemically active")
- Review: Surface engineering toward stable lithium metal anodes — [Sci. Adv. 2023, DOI 10.1126/sciadv.adf1550; PubMed 37018409](https://pubmed.ncbi.nlm.nih.gov/37018409/)
- Review: Self-Assembly Monolayer Inspired Stable Artificial SEI Design for Next-Generation LMBs — [PubMed 37079652](https://pubmed.ncbi.nlm.nih.gov/37079652/)
- Review: Supramolecular chemistry for LMBs (field-driven interfacial ordering) — [CCS Chem. 2025](https://www.chinesechemsoc.org/doi/10.31635/ccschem.025.202506033)
- Review: Thiols and derivatives in Li batteries (1.2–1.6 V reduction) — [IntechOpen](https://www.intechopen.com/online-first/1263212)
- Review: lithiophilic modification of carbon hosts — [PMC12113940](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12113940/)
- Review: unified interphase engineering (SEI in batteries and supercapacitors) — [arXiv 2511.23005](https://arxiv.org/pdf/2511.23005)

### Inferences (consolidated gap list — synthesized, not quoted from reviews)
1. Mechanistic: whether a given SAM at 0 V vs Li+/Li is (a) intact, (b) desorbed, (c) converted to an inorganic residue (Li2S/Cu2S, Li3PO4/Li-phosphonate, Li-silicate/LiF) is unresolved for every anchor class; the field uses "passivating" and "sacrificial" language interchangeably without discriminating tests.
2. Characterization: no operando SAM-fate study (EQCM mass loss, in situ SERS/IRRAS of anchor bands, in situ AFM) and no potential-resolved ex situ XPS/ToF-SIMS of anchor-element core levels; static-SIMS dose limits and Li volume change make ex situ snapshots unreliable.
3. Transport: no measured Li+/Zn2+ permeation through a SAM; tunneling-vs-defect question unaddressed; β ≈ 1.25 Å−1 implies long-chain SAMs must plate via defects or after breakdown.
4. Descriptors: no systematic chain-length / terminal-group / packing-density series with CE or nucleation overpotential; only one Li data point (MPA, 17.6 mV, 98.4%) and facet-energy-only Zn arguments.
5. Dipole vs chemistry confound: fluorinated groups supply both dipole and LiF; no isostructural, non-reactive dipole series; no KPFM/UPS work-function measurements on battery collectors.
6. QC: contact angle is the de facto sole metric; thickness/coverage/order/desorption-charge rarely reported — hampers reproducibility and meta-analysis.
7. Computation: no MLIP (fine-tuned MACE, UMA/OMat24, CHGNet, M3GNet) study of SAM–metal–alkali systems; no uMLIP validation of S/P/Si–metal anchor bonds under reduction; SEI-MLIP methodology (1300–4000 atoms, reaction analysis) is ready to extend.
8. Screening: HT/GNN/LLM-agent screening exists for electrolyte additives but not for anchored monolayers; SAM-specific targets (anchor bond energy, desorption potential, facet selectivity, normal dipole, Li adatom insertion barrier) are unused.
9. Chemistry coverage: no SAM work found on hard carbon, Na/K metal collectors, Mg/Ca (computational), or Al/Cu foils under real Li-ion cell conditions.
10. Texture: facet-selective SAM adsorption is validated only for Zn(002); Li texture itself is contested (110 vs random), and no SAM-driven Li/Na texture study is verifiable.

### Gaps (process)
- Open-question sections of the Small 2024, Sci. Adv. 2023 and 2023 artificial-SEI reviews could not be read; the writer should treat the list above as researcher inference.
- Searches not completed due to budget: Science 2022 abn1818 details; Nano Lett. 2024 desolvation layer; classical pinhole/redox-probe SAM permeation; phosphonate cathodic stability; Zn SAM dipole DFT; graphite SAM 2025 follow-ups; MLIP studies of thiolate/Au; perovskite SAM HT-screening precedents; Na-metal collector SAMs; thiol SAM S 2p after Li cycling. A follow-up turn could complete these.
