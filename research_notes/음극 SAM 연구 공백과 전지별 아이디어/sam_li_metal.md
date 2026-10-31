# Self-Assembled Monolayers (SAMs) at Li-Metal and Anode-Free Current-Collector Interfaces: Prior Work and Gaps

> **Research-session caveats (important for the report writer).**
> 1. In this cloud session almost every primary-source host (sciencedirect.com, onlinelibrary.wiley.com, pubs.acs.org, nature.com, pmc/pubmed.ncbi.nlm.nih.gov, science.org, arxiv.org, chemrxiv.org, osti.gov, doi.org, api.semanticscholar.org, researchgate.net) was **outside the network allowlist**, so full texts could not be fetched. All findings below come from search-result snippets/abstract excerpts surfaced by WebSearch (extended mode), which quotes abstract text. Numbers are therefore abstract-level and should be re-checked against the papers.
> 2. The per-turn WebSearch budget was exhausted after ~20 queries, so several planned follow-ups (BTA/Cu primary paper, MBI/Cu metrics, HMDS/Cu metrics, Gao 2020 DFT details, "Self-Assembled Monolayers for Batteries" 2021 review) could not be completed. These are listed under Gaps with explicit "could not verify" flags.
> 3. Items marked **[memory, unverified]** are from the researcher's background knowledge and were *not* confirmed in this session; treat as leads, not citations.

---

## Key Question 1: Which SAMs have been reported on Cu (and other current collectors) for anode-free / Li||Cu cells, and what did they achieve?

### Takeaway
The verified body of ex situ SAM-on-Cu work for Li plating is small (fewer than ten primary papers found) and is dominated by (i) thiol-anchored SAMs (alkanethiol-carboxylic acids, aromatic thiols, 2-mercaptobenzimidazole), (ii) N-heterocycle Cu–N anchored SAMs (benzotriazole, 5-methyl-1H-tetrazole), (iii) one electrochemically active aryl-sulfonyl-fluoride SAM (Gao/Wang, Nature Energy 2020) and (iv) one silane-type (HMDS) SAM. Only two of these (Gao 2020; Liang 2026) report pouch-level or ≥2 mAh cm⁻² data; Ag-coated Cu, Ni and stainless steel SAM studies for Li plating were **not found at all**.

### Cited Findings (study-by-study record)

**(1) Electrochemically active SAM: 1,3-benzenedisulfonyl fluoride on Cu — Gao, Rojas, Wang, Liu, Wang, Chen, Wang, Ngo, D. Wang (Penn State / Argonne), *Nature Energy* 2020, 5, 534–542, DOI 10.1038/s41560-020-0640-7**
- Molecule/substrate/method: SAM of electrochemically active 1,3-benzenedisulfonyl fluoride formed on the Cu current collector; the molecule forms a lithiophilic benzenesulfinate anion on Cu that guides Li nucleation when interfacial Li⁺ concentration is low (low T / high rate) — [Penn State Pure listing](https://pure.psu.edu/en/publications/low-temperature-and-high-rate-charging-lithium-metal-batteries-en/); [TechXplore](https://techxplore.com/news/2020-07-lithium-metal-batteries-temperatures.amp)
- Claimed mechanism: SAM reduction yields an SEI with a LiF-rich inner layer and amorphous outer layer that "seals" the Li surface; conventional SEI is non-passivating at low T — [Penn State Pure listing](https://pure.psu.edu/en/publications/low-temperature-and-high-rate-charging-lithium-metal-batteries-en/)
- Key metric (as relayed by secondary sources): Li plating/stripping on SAM-Cu retained 98.6% CE over 250 cycles at 6.0 mA cm⁻² / 6.0 mAh cm⁻² at −15 °C, bare Cu shorted; a 2.0 mAh cm⁻² full cell lasted 200 cycles at −15 °C with 45-min charge — [Sci. Adv. 2023 review (adf1550) via search excerpt](https://www.science.org/doi/10.1126/sciadv.adf1550); [X-mol listing](https://www.x-mol.com/paper/1275293115372564480); [Nature Energy PDF](https://www.nature.com/articles/s41560-020-0640-7.pdf)
- Related patent: US 12,658,447 "Stable metal anodes and batteries utilizing the same" describes electrochemically active (incl. sulfonate-terminated) SAMs on Cu that in situ form LiF nuclei — [USPTO](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12658447)
- Verification status: title/authors/DOI/journal verified; numeric metrics come from secondary excerpts (flag).

**(2) Thiol-anchored alkyl-carboxylic SAMs: 3-mercaptopropionic acid (MPA) and 11-mercaptoundecanoic acid (MUA) on planar Cu — *Electrochimica Acta* 2024 (article PII S0013468624002421), "Self-assembled monolayer regulates lithium nucleation and growth for stable lithium metal anodes"**
- Method: SAMs anchored by spontaneous thiol-to-Cu reaction, lowering surface energy; chain length and terminal group varied to tune electron/ion transfer and lithiophilicity — [ScienceDirect abstract](https://www.sciencedirect.com/science/article/abs/pii/S0013468624002421)
- Metric: Cu-MPA gives nucleation overpotential 17.6 mV and average CE 98.4% over >250 cycles at 0.5 mA cm⁻² / 1 mAh cm⁻² (Li||Cu half cell); short alkyl spacer + –COOH terminal was optimal, producing "large lithium nuclei" at lower overpotential — [ScienceDirect abstract](https://www.sciencedirect.com/science/article/abs/pii/S0013468624002421)
- Mechanism claimed: SAM "induces large lithium nuclei" / stable Cu–Li interface "from the very beginning of nucleation."
- Not reported in abstract: electrolyte, full/anode-free cell, post-plating SAM characterization, DFT. Authors/DOI not captured (flag).

**(3) Silane-type SAM: hexamethyldisilane (HMDS) on Cu foil — *Electrochimica Acta* 2024 (PII S0013468624013914), "Facile self-assembled monolayer deposition on copper foil for high-performance lithium-metal batteries"**
- Molecular single-layer modification of Cu with HMDS; homogeneous Li plating persists after extended cycling; authors attribute benefit to reduced gas-phase by-products of electrolyte decomposition — [ScienceDirect abstract](https://www.sciencedirect.com/science/article/abs/pii/S0013468624013914)
- Metrics/electrolyte not captured in snippet (flag; search budget exhausted before follow-up).

**(4) N-heterocycle SAM: 5-methyl-1H-tetrazole (MTZ) on Cu — Liang et al., *Advanced Materials* 2026, 38(43), e73830, DOI 10.1002/adma.73830, "Molecular Interfacial Engineering of Lithiophilic Protection Layer Realizing 600 Wh kg⁻¹ Anode-Free Lithium Metal Batteries"**
- MTZ bonds to Cu via strong Cu–N bonds forming a lithiophilic SAM; provides "high-density nucleation sites with molecular precision" and promotes preferential decomposition of LiNO₃ from the electrolyte to form a stable SEI; at 0.1 mAh cm⁻² plated Li is uniform on MTZ-Cu vs. irregular/aggregated on bare Cu — [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13432276/); [Wiley](https://advanced.onlinelibrary.wiley.com/doi/10.1002/adma.73830)
- Metrics: anode-free pouch cells of 8.44 Ah and 9.96 Ah (NCM811) reach 505.31 and 605.96 Wh kg⁻¹(total) for 100 cycles; 2.66 Ah LFP anode-free pouch at 292.01 Wh kg⁻¹ — [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13432276/)
- This is the only SAM-on-Cu study found that couples SAM chemistry to a *specific electrolyte additive* (LiNO₃) decomposition pathway.

**(5) 2-Mercaptobenzimidazole (MBI) on Cu foil — *ACS Omega* 2026, 11(28), 42854, "2-Mercaptobenzimidazole-Modified Copper Foil for Regulating Lithium Deposition in Anode-Free Lithium Batteries"**
- Thiol/azole bifunctional molecular layer improves interfacial stability and regulates Li deposition in anode-free cells — [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13393032/); [ACS Omega](https://pubs.acs.org/acsodf/article/11/28/42854/5202533/2-Mercaptobenzimidazole-Modified-Copper-Foil-for)
- Metrics not captured (flag).

**(6) Benzotriazole (BTA) on Cu — "Kang et al., 2020" as cited in the Frontiers in Chemistry 2022 review**
- BTA with lithiophilic N atoms used on Cu to guide homogeneous Li plating/stripping; anode-free capacity retention improved from 13.3% to ~73.3% after 50 cycles — [Frontiers Chem. 2022 review, PMC9158430](https://pmc.ncbi.nlm.nih.gov/articles/PMC9158430/)
- Primary paper (journal/DOI) **not located** before budget ran out (flag).

**(7) Reviews that catalogue the above**
- He et al., *Small* 2024, DOI 10.1002/smll.202403537, "Self-Assembled Molecular Layers as Interfacial Engineering Nanomaterials in Rechargeable Battery Applications": covers SAMs at electrode–electrolyte, solid/liquid-electrolyte and electrode–substrate interfaces; notes densely packed layers shield electrodes from H₂O/O₂ and metal-binding tail groups promote dendrite-free deposition — [Wiley](https://onlinelibrary.wiley.com/doi/10.1002/smll.202403537)
- *ACS Nano* 2026, 20(34), 23608 Perspective, "Self-Assembled Molecular Interfaces in Batteries: From Surface Passivation to Programmable Interphase Engineering": classifies SAM functions as (a) interfacial energy regulation, (b) passivation, (c) local micro-environment control, (d) interphase (SEI/CEI) construction; emphasizes molecular-scale thickness, negligible mass penalty — [ACS Nano](https://pubs.acs.org/ancac3/article-abstract/20/34/23608/5278442/Self-Assembled-Molecular-Interfaces-in-Batteries?redirectedFrom=fulltext)
- "Surface engineering toward stable lithium metal anodes", *Sci. Adv.* 2023 (adf1550): cites Gao 2020 SAM-on-Cu as low-temperature design; notes most artificial SEIs are thick (10–160 μm) and single-component — [Science](https://www.science.org/doi/10.1126/sciadv.adf1550); [PubMed](https://pubmed.ncbi.nlm.nih.gov/37018409/)
- "Rational Engineering of Anode Current Collector for Dendrite-Free Lithium Deposition", *Front. Chem.* 2022: states bare Cu is intrinsically lithiophobic with high nucleation barrier → spatially inhomogeneous nucleation; lists BTA among molecular modifiers — [PMC9158430](https://pmc.ncbi.nlm.nih.gov/articles/PMC9158430/)
- "Supramolecular Chemistry for High-Performance Lithium Metal Batteries", *CCS Chem.* 2025: describes "coordination-assisted self-assembly" of heteroatom molecules adsorbing on Li without major chemical change; cites an engineered Cu substrate reaching 99.82% CE in anode-free cells — [CCS Chemistry](https://www.chinesechemsoc.org/doi/10.31635/ccschem.025.202506033)
- A 2021 item titled "Self-Assembled Monolayers for Batteries" exists (ResearchGate listing) — journal/authors not verified [could not fetch] — [ResearchGate](https://www.researchgate.net/publication/353835000_Self-Assembled_Monolayers_for_Batteries)

**(8) Non-SAM molecular-layer benchmarks on Cu (useful comparators, not SAMs)**
- 2D polyamide / lithiated-Nafion interphase with Li⁺-adsorbing groups on Cu for anode-free cells — [Nature Materials 2025](https://www.nature.com/articles/s41563-025-02339-y)
- Monolayer amorphous carbon grown on Cu changes Li nucleation — [PMC12806237](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12806237/)
- Ag/Cu ion-exchange seed layer lowers nucleation barrier — [ACS AMI 2022](https://pubs.acs.org/doi/abs/10.1021/acsami.2c01980); Ag–Cu alloy collector, 100 cycles in LFP full cell — [Electrochim. Acta 2025](https://www.sciencedirect.com/science/article/abs/pii/S0013468625009843); Ag-modified Cu reduces nucleation overpotential 240→50 mV — [Front. Chem. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9158430/)
- Subsurface-vacancy-engineered, oxidation-resistant Cu for anode-free — [arXiv 2508.00236](https://arxiv.org/pdf/2508.00236)
- Current-collector comparison (slurry-based solid-electrolyte, in situ Li deposition): stainless steel gives best CE ≈98% vs Al, Ni, C-coated Al — [Batteries 2023, 9, 412](https://doi.org/10.3390/batteries9080412)

**(9) SAMs on other current collectors (Ni, stainless steel, Ag-coated Cu)**
- Targeted search returned **no** SAM study on Ni or stainless steel for Li plating; only oxide-modified Ni foam (97.7% CE over 330 cycles) — [Mater. Adv.](https://pubs.rsc.org/ma/article/6/20/7609/896287/Nickel-oxide-modified-nickel-foam-current); encyclopedic note that Cu, Ni, SS are stable Li hosts — [encyclopedia.pub](https://encyclopedia.pub/entry/1206)
- Alkanethiol SAM literature on Ni exists only for tribology — [Appl. Surf. Sci. 2021](https://www.sciencedirect.com/science/article/abs/pii/S0169433221010394)

### Inferences
- Across the six primary SAM-on-Cu papers, the anchoring chemistries used are only **thiol (–SH→Cu–S)**, **azole N (Cu–N)**, **sulfonyl fluoride (electro-reducible)** and **one silazane/silane**. No phosphonic-acid, carboxylic-acid-anchored (on CuOₓ), or catechol SAM on Cu for Li plating was found, even though these are standard on oxidized Cu in corrosion/electronics literature.
- Only Gao 2020 and Liang 2026 report ≥2 mAh cm⁻² or pouch data; the Electrochim. Acta 2024 thiol work is at 1 mAh cm⁻² / 0.5 mA cm⁻². "High-areal-capacity (≥3 mAh cm⁻²) + lean electrolyte" SAM data are essentially limited to Gao 2020 (6 mAh cm⁻², low T).
- Two different design philosophies coexist and are not reconciled: *lithiophilic head-out* SAMs that lower nucleation overpotential (MPA, MTZ, BTA, MBI) vs. *sacrificial/electroactive* SAMs that are consumed into a LiF-rich SEI (Gao 2020). No paper directly compares the two on the same Cu/electrolyte.

### Gaps
- Could not retrieve electrolyte formulations, full-cell data, or post-mortem characterization for the Electrochim. Acta 2024 (MPA/MUA; HMDS), ACS Omega 2026 (MBI) and BTA (Kang 2020) papers — publisher domains blocked, search budget exhausted.
- No SAM study on **Ag-coated Cu, Ni, stainless steel, or Al** for Li plating was located; no study on SAM + 3D Cu host.
- Could not confirm whether any SAM-on-Cu work used **localized high-concentration electrolytes (LHCE)** or **fluorinated ether electrolytes**; abstracts seen mention only LiNO₃-containing (Liang 2026) or unspecified electrolytes.

---

## Key Question 2: Is there direct evidence of whether SAMs survive Li plating (operando / XPS / ToF-SIMS)? What is unknown?

### Takeaway
No paper was found that directly tracks a SAM on Cu through Li plating with XPS/ToF-SIMS/cryo-TEM and establishes whether it persists, is displaced to the Li/SEI top surface, or is reductively desorbed. The only mechanistic claim with chemistry is Gao 2020, where the SAM is *designed* to be reduced (sulfonyl fluoride → sulfinate + LiF) and thus converted into SEI; survival of "inert" thiol/azole SAMs at 0 V vs Li⁺/Li is assumed, not demonstrated, in the abstracts seen.

### Cited Findings
- Gao 2020: SAM is "electrochemically active"; its reduction produces a benzenesulfinate anion on Cu and a LiF-rich inner SEI — i.e., SAM→SEI conversion is the intended mechanism, not survival — [Penn State listing](https://pure.psu.edu/en/publications/low-temperature-and-high-rate-charging-lithium-metal-batteries-en/); [TechXplore](https://techxplore.com/news/2020-07-lithium-metal-batteries-temperatures.amp)
- Liang 2026 (MTZ): claims SAM "promotes preferential decomposition of LiNO₃" to form SEI, implying the SAM is present during early SEI formation, but no statement about SAM fate after plating was captured — [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13432276/)
- Reductive stability of thiolate SAMs (aqueous): "Electrochemical Stability of Thiolate SAMs on Au, Pt, Cu" (ChemRxiv 2022): reductive stability order **Au < Pt < Cu**; on Cu, desorption competes with hydrogen evolution; on Au desorption proceeds via Au–S cleavage, so Cu–S may differ — [ChemRxiv](https://chemrxiv.org/doi/pdf/10.26434/chemrxiv-2022-dt7dk); [NSF PAR copy](https://par.nsf.gov/servlets/purl/10433337)
- Electrochemically assembled thiol SAMs on Cu: "potential-dependent processes such as SAM reductive desorption and Cu substrate oxidation should be considered" (HKUST thesis) — [HKUST](https://researchportal.hkust.edu.hk/en/studentTheses/electrochemical-assembly-of-thiol-based-monolayer-on-copper-for-e/)
- XPS reference for thiolate on Cu: thiolate S 2p ≈162.9 eV on mixed Cu/Cu₂O surfaces; thiols also reduce Cu oxide — [J. Electroanal. Chem. 2009](https://www.sciencedirect.com/science/article/abs/pii/S0022072809003672); [Thiol adsorption on/reduction of CuOₓ](https://www.researchgate.net/publication/299569044_Thiol_Adsorption_on_and_Reduction_of_Copper_Oxide_Particles_and_Surfaces)
- Baseline SEI-on-Cu chemistry without SAM: Menkin et al. (Cambridge), "Toward an Understanding of SEI Formation and Lithium Plating on Copper in Anode-Free Batteries" (J. Phys. Chem. C 2021 **[journal/author from memory, unverified]**): SEI on Cu contains Cu compounds (Cu oxides, LiF, Li₂CO₃) and plating is preferentially located where SEI composition favors it; XPS mapping used — [PMC8392351](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8392351/); [ChemRxiv 2020](https://chemrxiv.org/doi/pdf/10.26434/chemrxiv.12839792.v1)
- Operando/SIMS tools exist but have not been applied to SAM-Cu: in situ liquid SIMS + cryo-TEM of Cu anode in 1 M LiPF₆/EC-DMC during SEI formation — [OSTI 2305717](https://www.osti.gov/pages/servlets/purl/2305717); cryogenic ToF-SIMS resolving sub-nm SEI gradients and LiH on Li plated on Cu — [Microsc. Microanal. 2026](https://academic.oup.com/mam/article/32/4/ozag065/8736929); ToF-SIMS depth profiling of coatings on Li — [PMC11790834](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11790834/); in situ Li|SE ToF-SIMS notes very thin interphases may be below depth-profiling resolution — [Adv. Mater. Interfaces 2022](https://advanced.onlinelibrary.wiley.com/doi/10.1002/admi.202102387)
- Thermal/chemical stability data from other fields: alkanethiol SAM on Cu degradation by atomic H decreases with chain length (C6<C12<C18) (area-selective ALD context) — [Langmuir 2026 / PMC13421965](https://pmc.ncbi.nlm.nih.gov/articles/PMC13421965/)

### Inferences
- Thiolate SAMs on Cu are reductively desorbed in aqueous media at roughly −1 V vs Ag/AgCl **[memory, unverified]**; 0 V vs Li⁺/Li is ≈ −3 V vs SHE, far beyond that window, so persistence of an intact Cu–S–R monolayer *under* plated Li is thermodynamically doubtful. The more plausible outcomes are (a) desorption/reduction to Li thiolate (Li–S–R) incorporated into the SEI, (b) Li plating *beneath* the SAM with the SAM floating to the Li/SEI interface, or (c) partial survival in SAM defects. None of these have been discriminated experimentally for Li systems.
- The Gao 2020 strategy sidesteps the survival problem by making SAM consumption the mechanism; the MPA/MTZ/BTA/MBI papers implicitly need the SAM to act only during nucleation (first few nm), so their claimed mechanism is consistent even if the SAM does not survive — but this has not been tested.

### Gaps
- **No direct operando/XPS/ToF-SIMS evidence of SAM fate on Cu after Li plating** found (search: "thiol monolayer copper lithium plating XPS survive reductive desorption SEI" and "ToF-SIMS OR cryo-TEM self-assembled monolayer copper lithium deposition").
- No measured reductive-desorption potential of thiol/azole/phosphonate SAMs on Cu in **aprotic Li electrolytes** was found.
- No study of SAM behaviour across the **Cu-oxide reduction window (~1.5–0.5 V vs Li⁺/Li, [memory, unverified])** where the native Cu₂O/CuO anchoring layer for carboxylate/phosphonate SAMs would itself be reduced.
- Whether SAMs re-form or re-order after stripping (cycle-to-cycle persistence) is unaddressed.

---

## Key Question 3: Are there studies combining SAMs with Ag/Zn/Au lithiophilic seed layers, 3D hosts, or solid/polymer electrolytes?

### Takeaway
None were found. Seed-layer (Ag, Au, Pt, Zn) and 3D-host literatures for anode-free Li are mature, and SAM literature on Cu exists, but no paper combining them was located; the only SAM + solid-electrolyte item is an older (pre-2010) polymer-electrolyte/Li SAM study, and no SAM-on-garnet/sulfide work exists for Li wetting.

### Cited Findings
- Seed layers without SAMs: Ag/Cu ion exchange — [ACS AMI 2022](https://pubs.acs.org/doi/abs/10.1021/acsami.2c01980); Ag–Cu alloy — [Electrochim. Acta 2025](https://www.sciencedirect.com/science/article/abs/pii/S0013468625009843); Au, Pt, C seed layers for anode-free *solid-state* nucleation dynamics — [PMC10788862](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10788862/); ZnF₂-derived LiZn/LiF layer on Cu, >900 h at 3 mA cm⁻² — [PMC12097123](https://pmc.ncbi.nlm.nih.gov/articles/PMC12097123)
- 3D hosts without SAMs: graphene-modified Cu foam — [PMC6890847](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6890847/); super-3D lithiophilic Cu — [Mater. Today Energy 2023](https://www.sciencedirect.com/science/article/abs/pii/S2468606923000977); Cu foam with Cu₂O/LiF — [J. Energy Storage 2025](https://www.sciencedirect.com/science/article/abs/pii/S2352152X25030804)
- Solid electrolytes: garnet wetting is addressed with ALD Al₂O₃/ZnO or Li-alloy layers, not SAMs — [OSTI](https://www.osti.gov/servlets/purl/1433807); [Sci. Adv. 2017](https://www.science.org/doi/10.1126/sciadv.1601659); [SusMat 2024 review](https://onlinelibrary.wiley.com/doi/full/10.1002/sus2.187). Searching "self-assembled monolayer garnet OR solid electrolyte lithium wetting" returned **no** SAM study.
- Polymer electrolyte: "Characterization of self-assembled molecular layers at the polymer electrolyte/lithium electrode interface" (OSTI record; reports SAMs greatly slow Li passivation) — [OSTI 20134055](https://www.osti.gov/etdeweb/biblio/20134055) (date/journal not captured; likely pre-2010)
- Silane layer + PEO for all-solid-state Li: propyl-acrylate silane on Li combined with PEO (Delaporte et al., ChemistryOpen 2019) — [Wiley](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/open.201900021)
- Only fluoroalkylsilane SAM in Li batteries found is on **cathode** spinel (FAS17 / sub-2 nm SAM; 94–97% retention) — [Sci. Rep. 2016](https://www.nature.com/articles/srep31999); [Sci. Rep. 2018](https://www.nature.com/articles/s41598-018-30135-z)

### Inferences
- A SAM on a **Ag-coated Cu** collector is chemically natural (thiol–Ag SAMs are the best-characterized SAM system) yet absent in the Li literature; the same holds for Au-seeded collectors. This is a clear open space.
- Mixed SAM / patterned SAM to create lithiophilic–lithiophobic nucleation patterns (common in Cu electrodeposition on SAM-patterned Au) has not been transferred to Li.

### Gaps
- No SAM + Ag/Zn/Au seed, SAM + 3D host, or SAM + LLZO/sulfide/polymer (post-2010) study found.
- No SAM study on **Al or carbon-coated Al** collectors for anode-free solid-state cells.

---

## Key Question 4: What computational studies (DFT, MD, MLIP) exist on SAM–Cu–Li interfaces?

### Takeaway
Computation specific to SAM–Cu–Li is nearly nonexistent: the only verified example is the DFT in Gao et al. 2020 (co-authors Rojas and Ngo at Argonne; details not retrievable). DFT of thiol SAMs on Cu exists only in corrosion/solvent contexts, and ML-potential work on Li nucleation does not include organic monolayers or Cu.

### Cited Findings
- Gao et al. 2020 includes computational co-authors (T. Rojas, A. T. Ngo) — the paper is known to contain DFT of the electroactive SAM on Cu **[content of DFT from memory, unverified]** — [Penn State listing](https://pure.psu.edu/en/publications/low-temperature-and-high-rate-charging-lithium-metal-batteries-en/)
- DFT of n-octadecanethiol SAM on Cu(111) (corrosion context): adsorption energy −3.37 eV, stronger than solvent molecules; solvent polarity controls SAM quality — [PMC6017570](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6017570/)
- DFT of D-alaninol self-assembly on Cu(100) (surface science) — [Surf. Sci. 2008](https://www.sciencedirect.com/science/article/abs/pii/S0749603608001973)
- ML-enhanced constant-potential simulation of dendrite nucleation at Li–electrolyte interface (no Cu, no SAM) — [arXiv 2406.14025](https://arxiv.org/pdf/2406.14025)
- ML potentials for Li–sulfide-electrolyte mechanochemistry (no SAM) — [arXiv 2609.01088](https://arxiv.org/pdf/2609.01088); rough Cu–water interfaces with MLIP (no Li) — [arXiv 2509.17833](https://arxiv.org/pdf/2509.17833)
- Nucleation on SAM-covered metals (Cu electrodeposition on thiol/Au) is defect-mediated at SAM–substrate interface — [JACS 2009, ja907003w](https://pubs.acs.org/doi/10.1021/ja907003w) (analogy, not Li)
- Electrochemical SEI-formation modeling on Cu without SAM — [arXiv 1903.09593, "Design principles for self-forming interfaces"](https://arxiv.org/pdf/1903.09593)
- Review statement: mechanisms by which molecular layers tune interfacial electronic structure and the design rules "are still poorly understood" — [Small 2024](https://onlinelibrary.wiley.com/doi/10.1002/smll.202403537)

### Inferences
- There is **no published MLIP (e.g., DeePMD/MACE/NequIP) model of a SAM/Cu/Li ternary interface**, no AIMD or ML-MD of Li⁺ desolvation/transport through a SAM on Cu, and no computed Li nucleation energetics (adatom binding, cluster formation energy, work of adhesion) on SAM-terminated Cu vs bare Cu. Searching "DFT self-assembled monolayer Cu Li adsorption nucleation" and "machine learning interatomic potential lithium copper interface organic monolayer" produced none.
- Existing DFT infrastructure for thiol/Cu(111) (adsorption geometries, c(4×2)-type packing) could be reused directly to build such models.

### Gaps
- Could not retrieve what exactly was computed in Gao 2020 (reduction potentials? Li binding to sulfinate?).
- No computational treatment found of: SAM reductive desorption potentials on Cu in aprotic media; SAM→Li-thiolate/Li-alkoxide conversion thermodynamics; Li–SAM–Cu work of adhesion / wetting; SAM effect on electron tunneling (SEI growth) vs Li⁺ permeation.

---

## Key Question 5: SAMs on Li metal itself and electrolyte-formed monolayers; what do reviews list as open questions?

### Takeaway
SAMs directly on Li metal are better represented than on Cu, but almost exclusively via Li-reactive headgroups (phosphonic acids, chlorosilanes, alkoxysilanes); fluorinated SAMs on Li and in situ additive-formed monolayers on Cu were not found. Reviews consistently name mechanism (electronic-structure tuning), design rules, characterization of sub-nm layers, and thick/single-component artificial SEIs as the open problems.

### Cited Findings — SAMs on Li metal
- Octadecylphosphonic acid (ODPA) SAM on Li–CNT composite via lithium phosphate bonding; air/humid-air stable, slurry-coatable; LFP full cell 250 cycles, CE 98.8% at 33% DOD, 1 C; 82.5% retention after 250 cycles at 2:1 anode/cathode ratio; self-assembly studied by SFG — "Self-Assembled Monolayer Enables Slurry-Coating of Li Anode", *ACS Cent. Sci.* 2019, 5(3), 468, DOI 10.1021/acscentsci.8b00845 (Cui group **[author list from memory, unverified]**) — [ACS](https://pubs.acs.org/doi/10.1021/acscentsci.8b00845); [PubMed 30937374](https://pubmed.ncbi.nlm.nih.gov/30937374/)
- SAM-inspired silane dip-coating of Li with "MPDMS" to form inorganic-rich SEI; stable Li plating/stripping >500 cycles in carbonate electrolyte vs pristine Li failing after ~300 — "Self-Assembly Monolayer Inspired Stable Artificial SEI Design for Next-Generation LMBs", 2023 — [PubMed 37079652](https://pubmed.ncbi.nlm.nih.gov/37079652/) (journal and MPDMS expansion not captured; flag)
- Chlorosilane coatings on Li: self-terminating R₃Si–Cl reaction with surface Li–OH; shorter alkyl chains give lower R_CT, stable after 7 days immersion — "Stabilization of lithium metal anodes using silane-based coatings" (Langmuir 2007, Marchioni/Dunn/Wudl **[memory, unverified]**) — [ResearchGate](https://www.researchgate.net/publication/241279798_Stabilization_of_lithium_metal_anodes_using_silane-based_coatings); Frontiers in Materials 2019 silane-modified Li (LiCl + H₂ by-products via XPS) — [Frontiers](https://www.frontiersin.org/journals/materials/articles/10.3389/fmats.2019.00267/full)
- Silane + PEO for all-solid-state Li — [ChemistryOpen 2019](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/open.201900021)
- Older patents: chlorotrimethylsilane two-step Li treatment — [US 8,383,190](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/8383190); hydrophobic layer under metal-fluoride layer — [US 6,964,829](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/6964829)
- Search for fluoroalkylsilane/fluorinated SAM directly on **Li metal** returned only cathode examples — none on Li (gap) — [Sci. Rep. 2016](https://www.nature.com/articles/srep31999)

### Cited Findings — electrolyte-additive-formed monolayers vs ex situ SAMs
- Only indirect: Liang 2026 shows an ex situ SAM steering additive (LiNO₃) decomposition — [PMC13432276](https://pmc.ncbi.nlm.nih.gov/articles/PMC13432276/); the CCS Chem 2025 review frames additive adsorption as "coordination-assisted self-assembly" on Li — [CCS Chem.](https://www.chinesechemsoc.org/doi/10.31635/ccschem.025.202506033)
- Searching "electrolyte additive adsorbs copper forms molecular layer ... self-assembled anode-free" found **no** paper demonstrating an in situ additive-derived SAM on Cu with comparison to an ex situ SAM.

### Cited Findings — analog systems (out of scope but mechanistically informative)
- n-Alkanethiol SAMs on Fe electrodes suppress HER >90%, chain-length dependent, via reduced accessible area — [JACS 2012](https://pubs.acs.org/doi/10.1021/ja3095119)
- PFDT SAM on Zn: self-discharge 0.15% vs 7.62% at 48 h — [ACS Omega 2026](https://pubs.acs.org/acsodf/article/11/2/3569/5077631/PFDT-Self-Assembled-Monolayers-for-Stable-Zinc)

### Cited Findings — open questions stated by reviews
- Small 2024: detailed mechanisms by which SAMs tune interfacial electronic structure, and design rules, remain poorly understood — [Wiley](https://onlinelibrary.wiley.com/doi/10.1002/smll.202403537)
- Sci. Adv. 2023: most artificial SEIs are too thick (10–160 μm) and single-ingredient; multi-component, thin layers needed — [Science](https://www.science.org/doi/10.1126/sciadv.adf1550)
- ACS Nano 2026 Perspective: calls for "programmable" interphase engineering—moving from passive passivation to SAMs that dictate SEI chemistry — [ACS Nano](https://pubs.acs.org/ancac3/article-abstract/20/34/23608/5278442/Self-Assembled-Molecular-Interfaces-in-Batteries?redirectedFrom=fulltext)
- Design-principles review (Matter/Cell Rep. Phys. Sci.): artificial SEI needs mechanical stability, uniform Li⁺ transport, chemical passivation — [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S2666386420301181); mechanics review — [ChemElectroChem 2024](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/celc.202400278)

### Inferences
- The Li-metal SAM literature uses **reactive** anchoring (P–O–Li, Si–O–Li) and targets air-stability/processing (Cui 2019) rather than electrochemical nucleation control; the Cu SAM literature targets nucleation. Nobody has used the same molecule on both Li and Cu to separate "nucleation" from "SEI/passivation" effects.
- "Electrolyte-additive-formed monolayer" (in situ adsorbed, dynamic) vs "ex situ covalent SAM" (static, may desorb) has never been compared head-to-head on Cu under identical electrolytes.

### Gaps
- Fluorinated SAMs (–CF₃-terminated thiols/silanes/phosphonates) on Li or Cu for Li plating: **none found**.
- Carboxylic-acid-anchored (on CuOₓ) and phosphonic-acid-anchored SAMs on Cu for Li: none found (phosphonic acid SAMs on Cu are known in etch-inhibition — [Langmuir 2024/2025](https://pubs.acs.org/doi/10.1021/acs.langmuir.4c05097) — but not for Li).
- Thiol SAMs directly on Li (Li–thiolate monolayers): none found.
- Journal/authors for the MPDMS 2023 paper and Frontiers 2019 silane paper not captured.

---

## Consolidated Gap Map (for the report writer)

| Axis | Studied (verified above) | Not found / open |
|---|---|---|
| Molecule class on Cu | thiol-COOH (MPA/MUA), aromatic thiol-azole (MBI), azoles (BTA, MTZ), aryl-SO₂F (Gao), HMDS | alkanethiol chain-length series for Li, fluorinated thiols/silanes, phosphonic acids, carboxylic acids on CuOₓ, catechols, mixed/patterned SAMs, zwitterionic/ionic SAMs, Li-salt-anion-terminated SAMs (other than sulfonate) |
| Substrate | planar Cu foil only | Ag-coated Cu, Au/Pt-seeded Cu, Ni, stainless steel, Al, 3D Cu foams, carbon-coated collectors |
| Electrolyte | unspecified carbonate/ether; LiNO₃-containing (MTZ) | LHCE, fluorinated-ether, high-concentration, ionic-liquid, polymer, LLZO, sulfide; electrolyte–SAM compatibility (SAM dissolution in ethers) |
| Conditions | ≥3 mAh cm⁻² + low T only in Gao 2020; pouch only in Gao 2020 & Liang 2026 | lean electrolyte (≤3 g Ah⁻¹) with SAM; low T for non-electroactive SAMs; high current (≥5 mA cm⁻²) for thiol/azole SAMs; calendar aging / self-discharge of SAM-Cu anode-free cells |
| Mechanism | nucleation-overpotential lowering (MPA, MTZ, BTA); SAM→LiF-rich SEI (Gao) | SAM fate at 0 V vs Li⁺/Li (desorption vs Li-thiolate vs burial); reductive-desorption potentials in aprotic media; Cu-oxide reduction effect on anchoring; SAM effect on Li wetting angle / adhesion energy; cycle-to-cycle persistence; head-out vs sacrificial SAM comparison |
| Characterization | ex situ SEM, CE curves; XPS of SEI in Gao 2020 | operando XPS/ToF-SIMS/cryo-TEM of SAM through first plating; isotopic labelling (²H/¹³C SAM) to trace SAM fragments in SEI; EQCM of SAM desorption |
| Modeling | DFT in Gao 2020 (details unretrieved); thiol/Cu DFT from corrosion field | DFT/MLIP of SAM–Cu–Li ternary interface; Li adatom/cluster energetics on SAM-Cu; Li⁺ permeation through SAM; constant-potential ML-MD of SAM reduction |

## Sources consulted that could not be opened (network allowlist) — flagged for verification
- sciencedirect.com (Electrochim. Acta 2024 ×2, Surf. Sci. 2008, J. Electroanal. Chem. 2009), onlinelibrary.wiley.com (Small 2024, Adv. Mater. 2026, ChemistryOpen 2019), pubs.acs.org (ACS Nano 2026, ACS Cent. Sci. 2019, ACS Omega 2026, Langmuir), nature.com (Nat. Energy 2020, Nat. Mater. 2025), pmc/pubmed (PMC13432276, PMC13393032, PMC9158430, PMC8392351, PMID 37079652), science.org (Sci. Adv. 2023), arxiv.org, chemrxiv.org, osti.gov, par.nsf.gov, frontiersin.org, researchgate.net, x-mol.com, doi.org, api.semanticscholar.org, encyclopedia.pub.
