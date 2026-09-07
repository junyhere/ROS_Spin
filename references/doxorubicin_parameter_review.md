# Doxorubicin semiquinone–triplet-oxygen targeted evidence-gap review

**Audience:** ROS_Spin developers and scientific reviewers  
**Review date:** 2026-09-06  
**Scope:** doxorubicin (adriamycin) semiquinone, ground-state triplet O2, their putative solution encounter, and downstream superoxide/H2O2 validation. Established background is omitted. “Unavailable” means no defensible condition-matched number was found; it does not mean zero. “Not reported in the accessible source” and “not investigated” are used separately.

## Direct answer

The evidence fills bulk aqueous electron-transfer kinetics, an isotropic semiquinone EPR signal and linewidth, environment-specific context for O2, species-resolved aqueous superoxide dismutation, several SOD rates, and quantitative ex-vivo/cellular ROS validation datasets. It does **not** fill the parameters that would make the spin encounter predictive: accessible numerical doxorubicin hyperfine tensors, encounter-specific O2 g/D/E, doxorubicin-semiquinone T1/T2/Tm, geometry, J, dipolar coupling, residence/escape rates, a prepared coherent state, or spin-resolved kD and kQ. No defensible physical range for kQ/kD exists. Gas, matrix, protein, membrane, closed-shell, and related-anthracycline evidence remains inactive context or sensitivity-only evidence.

The implementation therefore has two policy-enforced modes. `evidence_backed`
uses only measured, calculated, fitted, or explicitly approved surrogate values
and refuses encounter-level yield predictions. `sensitivity` requires explicit
opt-in plus a named scenario, and labels all output non-predictive. Disabled
zero-valued terms are algebraic placeholders rather than evidence of zero.

## Parameter provenance

The complete row-level ledger is `configs/parameter_provenance.csv`. Charge descriptions follow the source; missing protonation or formal charge is not inferred.

| Symbol / definition | Value and unit | Species; environment; T; pH | Method; uncertainty; status | Source | Suitability, limitation, intended location |
|---|---:|---|---|---|---|
| gSQ, isotropic g | 2.0035 | Doxorubicin semiquinone, protonation not assigned; N2-purged 100 mM phosphate; T not reported; pH 7.5 | CW EPR; not reported; measured | [Kalyanaraman et al. 1991](https://doi.org/10.1016/0003-9861(91)90023-C) | Initialization/sensitivity only; no tensor. `active_model.encounter.g_radical`. |
| ΔHpp,SQ, one-line width | 3.25 G | Same; 100 μM drug, 400 μM xanthine, 0.2 U XO in 2 mL; T not reported; pH 7.5 | CW EPR, 1 G modulation, 40 G scan, 1 mW; uncertainty not reported; measured | same | Spectral validation only. Unresolved hyperfine plus other broadening prevents unique T2. Evidence metadata, not relaxation. |
| ΔHpp,SQ, independent width | 3.2 G | Doxorubicin semiquinone; hypoxic 50 mM Tris, 1 mM drug, 400 μM hypoxanthine, 0.1 U/mL XO; 310 K; pH 7.4 | CW EPR; uncertainty not reported; measured | [Hasinoff et al. 2015](https://doi.org/10.1124/mol.115.098798) | Spectral validation only; same T2 caveat. Evidence metadata. |
| ASQ, hyperfine | Numeric values unavailable in accessible source | Adriamycin and daunomycin semiquinones; solvent, T, pH and charge not reported in accessible abstract | EPR, ENDOR, TRIPLE resonance, partial deuteration; uncertainty not reported; measured study identified | [Jülich et al. 1988](https://doi.org/10.1002/mrc.1260260812) | Assignments were investigated, but numerical extraction is unavailable. Keep `effective_hyperfine_rad_s` disabled. |
| kO2,total | (3.5 ± 0.4)×10^8 and (1.7 ± 0.2)×10^8 M^-1 s^-1 | Adriamycin semiquinone + neutral O2; aqueous phosphate/formate or borate/OH-; room T; pH 6.0 and 11.5 | Electron pulse radiolysis/transient absorption; measured total-channel rates | [Land et al. 1985](https://doi.org/10.1038/bjc.1985.74) | Quantitative only for matched bulk solution; no spin resolution or encounter conversion. `measured_parameters`, validation only. |
| kO2,total, earlier study | 4.4×10^7 M^-1 s^-1 | Adriamycin semiquinone + O2; aqueous buffer; exact buffer, T and pH not reported in accessible source | Pulse radiolysis; uncertainty not reported; measured | [Svingen & Powis 1981](https://doi.org/10.1016/0003-9861(81)90263-0) | Differs 4–8× from conditioned 1985 values; conditions prevent reconciliation. Never kD/kQ. |
| t1/2,chemical; travel | 50 μs; <0.6 μm anaerobic; calculated 8 μs and <0.1 μm in air-saturated buffer | Same; aqueous; T/pH not reported in accessible source | Pulse radiolysis; uncertainty not reported; lifetime measured, air values/distance author-calculated | same | Chemical loss/travel, not spin T1/T2, pair separation, encounter radius, or cage lifetime. Context only. |
| λ0; γ0; derived D, gas O2 | 1.9847530; -8.42930×10^-3; 3.969506 cm^-1 | Neutral 16O2 X3Σg-, v=0 gas; pH N/A; spectroscopy T not extracted | Microwave/optical fit; uncertainty unavailable in accessible table; fitted; D=2λ convention conversion | [NIST/HITRAN](https://hitran.org/media/refs/HITRAN-1973.pdf) | Quantitative only for gas rovibronic spectroscopy; prohibited as active solution values. |
| Dmatrix | 3.572(3) cm^-1, axial | Neutral O2 in solid air condensed at 5 K; pH N/A | 94–550 GHz EPR; parenthetical fit uncertainty; measured/fitted | [Pardi et al. 2000](https://doi.org/10.1006/jmre.2000.2175) | Solid-only; do not transfer. |
| Dmatrix,N2 | 3.497 cm^-1 (104.9 GHz); g=2 fit | O2 in solid N2; 15 K; pH N/A | HF-EPR; uncertainty not reported; measured/fitted | [van der Horst & van Bentum 2001](https://doi.org/10.1016/S0921-4526(00)00615-3) | Temperature dependent; signal vanishes above 25 K. Matrix-only. |
| T1,O2,solution | approximately 7.5 ps | Neutral O2 in water, water/glycerol and seven small-molecule solvents; room T; pH not reported/not applicable | Proton magnetic-relaxation dispersion; uncertainty not reported in abstract; indirectly fitted | [Teng et al. 2001](https://doi.org/10.1006/jmre.2000.2219) | Bulk-fluid scale only; not encounter relaxation or static ZFS. |
| KO2,protein; τrot | 48 ± 7 M^-1 (Kd 21 ± 3 mM); 0.164 ± 0.006 and 1.41 ± 0.02 ps | Neutral O2 in T4 lysozyme L99A cavities; 50 mM phosphate/25 mM NaCl; 298 K; pH 5.5 | Gas-pressure NMR + MD; K measured/fitted, τ MD-fitted | [Kitahara et al. 2016](https://doi.org/10.1038/srep20534) | Protein context: rapid motion supports averaging, not a doxorubicin/O2 tensor. |
| pKa,HO2 | 4.88 ± 0.10 | HO2•/O2•- in oxygen-saturated formic acid/formate-buffered water; T not reported; pH 0–13 examined | Pulse radiolysis; measured | [Behar et al. 1970](https://doi.org/10.1021/j100711a009) | Quantitative aqueous speciation under matched conditions. Downstream evidence metadata. |
| kHH / kHA / kAA | 0.76×10^6 / 8.5×10^7 / <100 M^-1 s^-1 | HO2•+HO2• / HO2•+O2•- / O2•-+O2•- in same system | Pulse-radiolysis decay; rate uncertainties not reported; measured/bounded | same | Use species-resolved pH-dependent chemistry; neutral-pH effective k is not an anion–anion elementary constant. |
| kSOD half-reactions | (1.2 ± 0.2)×10^9 and (2.2 ± 0.4)×10^9 M^-1 s^-1 | Bovine CuZnSOD oxidation/reduction by O2•-; aqueous O2/EDTA/formate; T/pH not reported in accessible record | Pulse radiolysis, 650/300 nm; measured | [Klug-Roth et al. 1973](https://doi.org/10.1021/ja00790a007) | Quantitative for this enzyme only; surrogate otherwise. Keep inactive. |
| kSOD,overall | (2.37 ± 0.18)×10^9 M^-1 s^-1 | Native bovine CuZnSOD + O2•-; aqueous; 298 K; pH/buffer not reported in accessible abstract | Pulse radiolysis plus optical/EPR; measured | [Fielden et al. 1974](https://doi.org/10.1042/bj1390049) | Exact preparation benchmark, not universal SOD1/SOD2. Keep inactive. |

## Eight-topic evidence-gap matrix

| # | Topic | Best evidence and disconfirming check | Status | Remaining gap |
|---:|---|---|---|---|
| 1 | Semiquinone hyperfine | Exact EPR/ENDOR/deuteration paper identified; aqueous CW EPR is an unresolved 3.2–3.25 G line. | **Partly filled** | Numerical table, signs, tensors, isotope assignments and full conditions remain unavailable from accessible full text. |
| 2 | Triplet O2 by environment | Gas λ/γ and frozen-matrix D are nontransferable; solution O2 has ~ps T1; protein O2 rotates sub-ps/ps; membrane work measures probe collisions, not O2 g/D/E. | **Context filled; encounter open** | No liquid, protein/membrane, or transient doxorubicin-SQ/O2 encounter g/D/E tensor. |
| 3 | Semiquinone T1/T2/dephasing | Only unresolved CW linewidths found; chemical lifetime is not spin relaxation. | **Unavailable** | Direct pulse-EPR T1 and T2/Tm under matched chemistry. |
| 4 | Geometry, separation, J, dipolar, encounter/diffusion | Bulk lifetime/travel bounds exist. Land found no SQ–SQ/parent association at pH 5, but did not test SQ–O2. | **Unavailable for exact encounter** | r/orientation, J(r), dipolar tensor, Kassoc, kon/koff, ksep and cage time. |
| 5 | Coherent spin-correlated D/Q encounter | No exact-system transient polarization, quantum beats, CIDEP/CIDNP, field-dependent branching, or D/Q populations. | **Not investigated in retrieved exact-system studies; no direct paper found** | Spin-polarized transient EPR or field-dependent state-resolved products. |
| 6 | Spin-resolved kD/kQ | Exact-system rates are total bulk rates; degeneracy does not determine kinetics. | **Unavailable** | Independent state populations and channel yields/rates; no physical ratio interval. |
| 7 | Spontaneous/SOD dismutation | Primary HO2/O2•- and bovine CuZnSOD rates available. Recent reevaluation rejects commonly misassigned large anion–anion rates and singlet-O2 claims. | **Aqueous filled; biological transfer open** | Matched human SOD1/SOD2, compartment, pH, ionic strength, crowding and competing sinks. |
| 8 | Superoxide/H2O2 validation | Absolute organelle/cell fluxes and PC3 H2O2 exist; 5-iminodaunorubicin is a negative control; raw MitoSOX/DCF is nonspecific. | **Partly filled** | Clinically matched in-vivo time courses with calibrated species-specific analytics. |

Searches stopped after exact-title/DOI chasing and disconfirming queries repeatedly returned these papers, reviews repeating them, inaccessible tables, or wrong-species analogues. More generic quinone/radical-pair values would be weaker evidence.

## Minimum defensible Hamiltonian and D/Q mixing

For semiquinone spin **s** = 1/2 and O2 spin **S** = 1, in angular-frequency units:

`H/ℏ = βe B·(gSQ·s + gO2·S) + J s·S + s·Ddip·S + S·DO2·S + Σk s·Ak·Ik`.

The electronic space is 2×3 = 6 states, with `PD = (0.5 I - s·S)/1.5` and `PQ = (s·S + I)/1.5`. The defensible active baseline is the six-state space, exact projectors, a declared field, and measured/explicitly chosen Zeeman information. Spin-selective loss `K = kD PD + kQ PQ` and escape `kesc I` are model structure, but their numerical values are unavailable. J, Ddip, DO2, explicit Ak, and relaxation must remain disabled unless parameterized.

Isotropic `J s·S` and a common isotropic Zeeman term commute with total spin and do not mix D/Q. Unequal local Zeeman terms, nuclear hyperfine terms, O2 anisotropy/ZFS, and anisotropic electron–electron coupling can mix them. Their encounter magnitudes are unavailable, so no mixing frequency or field-effect magnitude is defensible. Also, `ρ = I6/6` is invariant under unitary evolution; oscillations from it require state-selective loss, relaxation, preparation, or another nonunitary process.

## Initial-state scenarios

1. **Unpolarized, uncorrelated benchmark:** `ρ0 = I6/6`; pD = 1/3 and pQ = 2/3 from degeneracy only, not experiment.
2. **Pure-doublet benchmark:** `ρ0 = PD/2`.
3. **Pure-quartet benchmark:** `ρ0 = PQ/4`.
4. **Bounded mixture:** 0 ≤ pD ≤ 1, sensitivity only.
5. **Chemically prepared coherent state:** unavailable.

## kQ/kD conclusion

`kD`, `kQ`, and `kQ/kD` are unavailable. The bulk `kO2,total` cannot be decomposed without state populations and channel yields. No defensible physical range or ordering exists. Preserve only dimensionless, explicitly illustrative benchmarks such as 0, 0.01, 0.1 and 1; these are not a confidence interval, prior, or chemical range and require opt-in. Equality is the spin-independent null model.

## Reaction network

1. `Q + e- -> SQ` (formation is preparation specific).
2. `SQ + O2 <=> [SQ···O2]D,Q` (association/separation unavailable).
3. `[SQ···O2]D -> Q + O2•-` with unavailable `kD`.
4. `[SQ···O2]Q -> Q + O2•-` with unavailable `kQ`.
5. `[SQ···O2]D,Q -> SQ + O2` with unavailable `kesc`.
6. `2 HO2• -> H2O2 + O2` and `HO2• + O2•- + H+ -> H2O2 + O2`, with explicit speciation.
7. `SODox + O2•- -> SODred + O2`; `SODred + O2•- + 2H+ -> SODox + H2O2`.
8. `H2O2 -> loss/transport`, compartment-specific and unavailable.

H2O2 is downstream, requires two superoxide equivalents plus removal/transport, and is not a spin population.

## Validation datasets

| Observable | Value; uncertainty | Dose, preparation, normalization | Time, T, pH/buffer, oxygen; method | Status, suitability, limitation | Source |
|---|---|---|---|---|---|
| Superoxide, beef-heart SMP | 1.6 ± 0.2 to 69.6 ± 2.7 nmol min^-1 mg^-1 | Control vs 90 μM doxorubicin; bovine-heart SMP | Initial rate; T, pH/buffer, oxygen not reported in accessible abstract; SOD-inhibitable cytochrome-c reduction | Measured mean ± SE; matched-preparation quantitative only | [Doroshow & Davies 1986](https://doi.org/10.1016/S0021-9258(17)35747-2) |
| H2O2, beef-heart SMP | undetectable to 2.2 ± 0.3 nmol min^-1 mg^-1 | Control vs 200 μM doxorubicin; per mg SMP | Timing/T/pH/buffer/oxygen not reported in accessible abstract; catalase-released O2 | Measured mean ± SE; detection limit not reported | same |
| Superoxide, Ehrlich microsomes | 0.51 ± 0.26 (n=3) to 14.71 ± 1.43 (n=6) nmol min^-1 mg^-1 | Control vs 135 μM; 200 μg protein in 1 mL; NADPH | Initial rate, 310 K; 150 mM potassium phosphate pH 7.4, 100 μM EDTA; oxygen not separately reported; acetylated cytochrome-c assay | Measured mean ± SE; exact ex-vivo quantitative validation | [Doroshow 2019](https://doi.org/10.1155/2019/9474823) |
| Superoxide, Ehrlich mitochondria | 1.12 ± 0.15 (n=7) to 7.29 ± 0.61 (n=10) nmol min^-1 mg^-1 | Control vs 135 μM; 100 μg/mL; NADH; 4 μM rotenone, 5-min preincubation | Initial rate, 310 K; 250 mM sucrose/20 mM HEPES pH 8.2/100 μM EDTA; oxygen not separately reported; same assay | Measured mean ± SE; alkaline/rotenone-blocked, not unqualified physiology | same |
| Superoxide, Ehrlich nuclei | 0.36 ± 0.05 (n=7) to 3.29 ± 0.33 (n=12) nmol min^-1 mg^-1 | Control vs 135 μM; 200 μg nuclear protein; NADPH | Initial rate, 310 K; sucrose/20 mM HEPES pH 7.4/EDTA; oxygen not separately reported; same assay | Measured mean ± SE; isolated nuclei only | same |
| H2O2, Ehrlich microsomes / mitochondria / cells | 3.42 ± 0.48 (n=4) / 2.6 ± 0.6 (n=3) nmol min^-1 mg^-1 / 0.64 ± 0.04 (n=3) nmol min^-1 per 10^7 cells; controls undetectable | 135 / 135 / 400 μM; fractions / 1.5×10^7 cells in 3 mL | Linear 10–30 min, 310 K; fraction buffers above; cells in Chelex-PBS, pH not reported, air-bubbled 30 min; catalase-released O2 | Measured mean ± SE; cell signal is extracellularly accessible H2O2; detection limits not reported | same |
| Intracellular apparent H2O2, PC3 | 13 ± 4 pM basal to 51 ± 13 pM | 1 μM doxorubicin; PC3 human prostate cancer cells | 30 min, 310 K; culture pH/O2 not reported; catalase–aminotriazole kinetic assay | Measured/derived; medium Amplex Red did not rise and catalase overexpression did not improve resistance | [Wagner et al. 2005](https://doi.org/10.1016/j.abb.2005.06.015) |
| MitoSOX, H9c2 | 1.40 ± 0.04, 1.80 ± 0.08, 2.80 ± 0.08 fold control | 10, 20, 50 μM doxorubicin; rat H9c2; n=3 | 1 h, 310 K; HBSS + Ca/Mg + 1% BSA, pH/O2 not reported; 5 μM MitoSOX | Mean ± SD; sensitivity only, not absolute; drug autofluorescence observed | [Mukhopadhyay et al. 2007](https://doi.org/10.1016/j.bbrc.2007.04.106) |

At 135 μM, quinone-disabled 5-iminodaunorubicin gave 0.45 ± 0.15 nmol min^-1 mg^-1 in Ehrlich microsomes versus doxorubicin 14.71 ± 1.43, a negative control supporting quinone dependence. Conversely, bulk MitoSOX red fluorescence includes nonspecific products and requires chromatographic resolution of the superoxide-specific 2-hydroxy product ([Zielonka et al. 2009](https://doi.org/10.1016/j.freeradbiomed.2008.10.031)); it must not calibrate an absolute source.

## Supported claims

- Doxorubicin semiquinone forms in the reported reducing systems and transfers an electron to O2 in bulk solution.
- Exact hyperfine-assignment work exists, but its numeric table was inaccessible; biological aqueous preparations show an unresolved line.
- O2 magnetic behavior and relaxation are environment dependent.
- The 1/2⊗1 space contains doublet and quartet manifolds of dimensions 2 and 4.
- Condition-specific aqueous dismutation/SOD kinetics and preparation-specific ROS data can validate downstream chemistry.

## Prohibited claims

- A measured coherent, entangled, pure-D, pure-Q, or other D/Q birth state.
- Any numerical exact-system T1, T2/Tm, J, dipolar tensor, geometry, encounter time, kesc, kD, kQ, or physical kQ/kD range.
- Transfer of gas/matrix O2, protein-cavity, closed-shell doxorubicin, or daunomycin values into the active solution encounter.
- Conversion of CW linewidth to T2 without line-shape decomposition, or bulk M^-1 s^-1 rates/chemical lifetimes to intracomplex s^-1 rates.
- Treating a neutral-pH effective dismutation coefficient as an O2•-+O2•- elementary rate, or claiming singlet-O2 production here.
- Absolute superoxide from unseparated MitoSOX/DCF fluorescence, or quantitative field effects/ROS/biological outcomes from current inputs.

## Complete bibliography

- Arroyo, C. M.; Carmichael, A. J. *Free Radic. Biol. Med.* 9 (1990) 191–197. [doi:10.1016/0891-5849(90)90027-G](https://doi.org/10.1016/0891-5849(90)90027-G).
- Behar, D. et al. *J. Phys. Chem.* 74 (1970) 3209–3213. [doi:10.1021/j100711a009](https://doi.org/10.1021/j100711a009).
- Bielski, B. H. J. et al. *J. Phys. Chem. Ref. Data* 14 (1985) 1041–1100. [doi:10.1063/1.555739](https://doi.org/10.1063/1.555739).
- Doroshow, J. H. *Oxid. Med. Cell. Longev.* 2019, 9474823. [doi:10.1155/2019/9474823](https://doi.org/10.1155/2019/9474823).
- Doroshow, J. H.; Davies, K. J. A. *J. Biol. Chem.* 261 (1986) 3068–3074. [doi:10.1016/S0021-9258(17)35747-2](https://doi.org/10.1016/S0021-9258(17)35747-2).
- Fielden, E. M. et al. *Biochem. J.* 139 (1974) 49–60. [doi:10.1042/bj1390049](https://doi.org/10.1042/bj1390049).
- Forman, H. J.; Fridovich, I. *Arch. Biochem. Biophys.* 158 (1973) 396–400. [doi:10.1016/0003-9861(73)90636-X](https://doi.org/10.1016/0003-9861(73)90636-X).
- Hasinoff, B. B. et al. *Mol. Pharmacol.* (2015). [doi:10.1124/mol.115.098798](https://doi.org/10.1124/mol.115.098798).
- Jülich, T. et al. *Magn. Reson. Chem.* 26 (1988) 701–706. [doi:10.1002/mrc.1260260812](https://doi.org/10.1002/mrc.1260260812).
- Kalyanaraman, B. et al. *Arch. Biochem. Biophys.* (1991). [doi:10.1016/0003-9861(91)90023-C](https://doi.org/10.1016/0003-9861(91)90023-C).
- Kitahara, R. et al. *Sci. Rep.* 6 (2016) 20534. [doi:10.1038/srep20534](https://doi.org/10.1038/srep20534).
- Klug-Roth, D. et al. *J. Am. Chem. Soc.* 95 (1973) 2786–2790. [doi:10.1021/ja00790a007](https://doi.org/10.1021/ja00790a007).
- Land, E. J. et al. *Br. J. Cancer* 51 (1985) 515–523. [doi:10.1038/bjc.1985.74](https://doi.org/10.1038/bjc.1985.74); [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC1977136/).
- Mukhopadhyay, P. et al. *Biochem. Biophys. Res. Commun.* 358 (2007) 203–208. [doi:10.1016/j.bbrc.2007.04.106](https://doi.org/10.1016/j.bbrc.2007.04.106).
- Pardi, L. A. et al. *J. Magn. Reson.* 146 (2000) 375–378. [doi:10.1006/jmre.2000.2175](https://doi.org/10.1006/jmre.2000.2175).
- Stanbury, D. M. *J. Phys. Chem. B* 130 (2026) 6654–6659. [doi:10.1021/acs.jpcb.6c01603](https://doi.org/10.1021/acs.jpcb.6c01603).
- Svingen, B. A.; Powis, G. *Arch. Biochem. Biophys.* 209 (1981) 119–126. [doi:10.1016/0003-9861(81)90263-0](https://doi.org/10.1016/0003-9861(81)90263-0).
- Teng, C. L. et al. *J. Magn. Reson.* 148 (2001) 31–34. [doi:10.1006/jmre.2000.2219](https://doi.org/10.1006/jmre.2000.2219).
- van der Horst, E.; van Bentum, P. J. M. *Physica B* 294–295 (2001) 87–90. [doi:10.1016/S0921-4526(00)00615-3](https://doi.org/10.1016/S0921-4526(00)00615-3).
- Wagner, B. A. et al. *Arch. Biochem. Biophys.* 440 (2005) 181–190. [doi:10.1016/j.abb.2005.06.015](https://doi.org/10.1016/j.abb.2005.06.015).
- Zielonka, J. et al. *Free Radic. Biol. Med.* 46 (2009) 329–338. [doi:10.1016/j.freeradbiomed.2008.10.031](https://doi.org/10.1016/j.freeradbiomed.2008.10.031).

## Repository mapping

- Review: `references/doxorubicin_parameter_review.md`
- Provenance: `configs/parameter_provenance.csv`
- Evidence and sensitivity inputs: `configs/doxorubicin_parameters.json`
- Evidence gate and full six-state reference implementation: `spin_chemistry.py`
- Dual-mode command-line interface: `ROS.py`
- Scientific and policy validation: `tests/test_spin_chemistry.py`
