# Doxorubicin semiquinone–triplet-oxygen parameter review

**Audience:** ROS_Spin developers and scientific reviewers  
**Review date:** 2026-09-06  
**Scope:** an encounter between doxorubicin (adriamycin) semiquinone, spin 1/2, and ground-state O2, spin 1, followed by superoxide/H2O2 chemistry.

## Executive answer

The old two-qubit Bell-state/CZ model is not a representation of this chemical system. The electronic encounter space is 2 x 3 = 6 dimensional and decomposes into a doublet (rank 2) and quartet (rank 4). The corrected implementation therefore uses exact total-spin projectors, coherent density-matrix evolution, phenomenological relaxation, spin-selective loss, encounter escape, and a separate downstream kinetic network.

The literature supports doxorubicin semiquinone formation, transfer of an electron to oxygen, primary superoxide formation, an aqueous isotropic semiquinone signal at g = 2.0033, and a bulk semiquinone + O2 rate of 4.4e7 M^-1 s^-1. It does **not** presently provide an exact-system spin-resolved kD or kQ, encounter lifetime, T1/T2, resolved hyperfine tensor, exchange/dipolar coupling, or solution encounter geometry. Consequently this repository now supports mechanistic sensitivity analysis, not quantitative prediction of field effects or ROS yields.

## Evidence and parameter provenance

The machine-readable version is in `configs/parameter_provenance.csv`. “Unavailable” means that no exact, condition-matched primary value was located; it is not zero.

| Symbol | Result | Species and conditions | Method/status | Model use and limitation |
|---|---:|---|---|---|
| gSQ | 2.0033 | Adriamycin semiquinone in water; temperature/pH not reported in the abstract | ESR, measured ([Akman et al., 1990](https://doi.org/10.1016/0891-5849(90)90027-G)) | Isotropic initialization only; no tensor/hyperfine resolution |
| kbulk | 4.4e7 M^-1 s^-1 | Adriamycin semiquinone + O2, aqueous pulse radiolysis | Measured ([Kalyanaraman et al., 1981](https://doi.org/10.1016/0003-9861(81)90263-0)) | Bulk validation target; **not** an encounter-complex first-order kD |
| t1/2,chemical | 50 us | Adriamycin semiquinone, anaerobic aqueous medium | Measured in the same pulse-radiolysis study | Chemical lifetime, not encounter lifetime and not electron-spin T1/T2 |
| pKa,SQ | 2.9 | Adriamycin semiquinone in aqueous solutions over varied pH | Pulse radiolysis, measured ([Land et al., 1985](https://doi.org/10.1038/bjc.1985.74)) | Supports protonation context; the radical persisted 10–20 ms around pH 6–11 in that study |
| kspont | approx. 2e5 M^-1 s^-1 at pH 7.4 | Effective HO2/O2- dismutation in water near 298 K | Critically evaluated kinetics ([Bielski et al., 1985](https://doi.org/10.1063/1.555739)) | Downstream only; recompute/specify for actual pH |
| kSOD | 2e9 M^-1 s^-1 | Bovine erythrocyte CuZnSOD, pH 7.8–8.5 | Competitive kinetics, measured ([Klug et al., 1972/1973](https://doi.org/10.1016/0003-9861(73)90636-X)) | Surrogate unless the biological compartment contains that isoform/condition |
| hyperfine tensors | unavailable | Exact doxorubicin semiquinone under relevant encounter conditions | Literature gap | Explicit nuclei cannot be parameterized; effective fields are sensitivity-only |
| O2 g/ZFS in encounter | unavailable | Ground-state O2 in a chemically relevant fluid encounter | Literature gap | Gas/matrix values are not silently transferred to solution |
| kD, kQ/kD | unavailable | Spin-resolved doxorubicin-semiquinone/O2 encounter | Literature gap | Required for spin-yield prediction; use declared benchmarks only |
| kesc | unavailable | Same encounter | Literature gap | Required to convert dynamics into cage yields |
| T1/T2 | unavailable | Same encounter | Literature gap | Radical chemical decay cannot substitute for spin relaxation |
| J, Ddip, geometry | unavailable | Same encounter | Literature gap | Required for quantitative mixing/coupling calculations |
| kH2O2,loss | unavailable | Compartment-specific H2O2 | Out of scope without a compartment | Needed for H2O2 concentration predictions |

Related anthracycline or quinone measurements are not inserted as doxorubicin constants. For context only, daunorubicin pulse radiolysis showed a pH-7 semiquinone/hydroquinone/parent pseudo-equilibrium lasting hundreds of milliseconds ([Houée-Levin et al., 1985](https://doi.org/10.1016/0014-5793(85)80188-5)); the chemical species differs and this is not a spin-relaxation surrogate.

## Minimum defensible Hamiltonian

With radical spin **s** (s=1/2), oxygen spin **S** (S=1), and angular-frequency units:

H/hbar = beta_e B·(gSQ·s + gO2·S) + J s·S + s·Ddip·S + S·DO2·S + s·Ahf·I.

Required for a defined baseline:

- Zeeman terms and an explicitly chosen external field.
- A six-state Hilbert space and the exact total-spin projectors
  PD = (0.5 I - s·S)/1.5 and PQ = (s·S + I)/1.5.
- Spin-selective loss K = kD PD + kQ PQ and encounter escape kesc I.

Optional only when parameterized:

- Isotropic exchange J, anisotropic exchange/dipolar coupling Ddip, O2 zero-field tensor DO2, explicit semiquinone hyperfine tensors Ahf, and local relaxation superoperators.
- The code includes a three-component effective hyperfine field for sensitivity work. It is not a substitute for explicit nuclei and is labelled accordingly.

Terms capable of doublet–quartet mixing are those that fail to commute with total spin: unequal/anisotropic Zeeman interactions, semiquinone hyperfine coupling to nuclei, O2 ZFS, anisotropic exchange/dipolar coupling, and local relaxation. Isotropic exchange J s·S commutes with total spin and only separates doublet/quartet energies. A common isotropic Zeeman interaction likewise does not mix manifolds.

No exact-system values were found for the mixing terms. Therefore magnetic-field-effect magnitude, coherent oscillation frequency, and quantitative branching are prohibited conclusions at present.

## Initial-state scenarios

1. **Unpolarized oxygen-triplet encounter (default):** rho0 = I6/6, giving initial doublet/quartet weights 1/3 and 2/3 from degeneracy. This is the least-committal thermal encounter, not proof that the reacting cage is statistically formed.
2. **Pure doublet benchmark:** rho0 = PD/2. Tests the spin-allowed limiting channel.
3. **Pure quartet benchmark:** rho0 = PQ/4. Tests blocking, mixing, relaxation, or nonzero kQ.
4. **Chemically prepared state:** none supported for this exact encounter by the located evidence. Do not claim one without spin-polarized transient-EPR evidence.

## kQ/kD recommendation

No exact-system spin-resolved rate ratio was located. Use the dimensionless benchmark set **0, 0.01, 0.1, 1**:

- 0: ideal quartet blocking;
- 0.01 and 0.1: partial leakage benchmarks;
- 1: spin-independent null model.

This is a purely illustrative sensitivity set, not a confidence interval or chemically measured range. The qualitative spin-selection argument makes kQ <= kD a testable hypothesis, but the available doxorubicin evidence does not quantify it. Claims about selectivity require the results to remain distinguishable from the kQ/kD=1 null model over plausible relaxation and escape rates.

## Reaction-network specification

Keep chemical formation and downstream fate outside the encounter Hamiltonian:

1. Q + e- --kform--> SQ (enzyme/cofactor-dependent semiquinone formation).
2. SQ + O2 <=> [SQ...O2] (encounter formation; bulk diffusion/association model not yet parameterized).
3. [SQ...O2]D --kD--> Q + O2- (primary superoxide, doublet channel).
4. [SQ...O2]Q --kQ--> Q + O2- (quartet leakage/sensitivity channel).
5. [SQ...O2]D,Q --kesc--> SQ + O2 (escape/separation).
6. 2 O2- + 2 H+ --kspont--> H2O2 + O2.
7. 2 O2- + 2 H+ --SOD--> H2O2 + O2, represented as kSOD[SOD][O2-] in the dilute-substrate limit.
8. H2O2 --kloss--> products/transport (catalase, peroxidases, and compartment transport must be specified).

The encounter solver reports primary superoxide yield. It never maps a spin population directly to H2O2. `downstream_ros` integrates steps 6–8 separately and uses the stoichiometric factor: two superoxide give one H2O2.

## Supported claims

- Doxorubicin/adriamycin can form a semiquinone and transfer an electron to oxygen, producing superoxide; direct ESR/spin-trapping and pulse-radiolysis studies support the sequence.
- The radical–triplet encounter has doublet and quartet total-spin manifolds with dimensions two and four.
- A bulk aqueous second-order rate and isotropic aqueous g signal are available under the reported experimental conditions.
- Doxorubicin stimulates superoxide and H2O2 production in isolated cardiac systems; for example, submitochondrial measurements found matched initial O2 consumption and superoxide formation under comparable conditions ([Davies & Doroshow, 1986](https://pubmed.ncbi.nlm.nih.gov/3005279/)).
- H2O2 is downstream of superoxide dismutation and must not be equated with a spin-manifold population.

## Unsupported claims

- That the biological encounter is born coherent, entangled, pure-doublet, or pure-quartet.
- A numerical kQ/kD, encounter lifetime, escape rate, T1, T2, J, dipolar coupling, or resolved hyperfine tensor for this exact pair.
- That gas-phase/matrix O2 constants apply unchanged in aqueous or membrane encounters.
- That a particular interaction dominates doublet–quartet mixing.
- Quantitative magnetic-field effects, absolute superoxide yield, H2O2 concentration, cardiotoxicity, or therapeutic response from the current parameter set.
- That quantum hardware noise represents chemical relaxation, or that a quantum circuit provides computational advantage for this six-state model.

## Classical-versus-circuit validation

The implementation uses the classical six-state propagator as the reference. For coherent validation it embeds the qutrit as `00`, `01`, `10` in two qubits, reserves `11` as leakage, and embeds the 6x6 unitary into an 8x8 three-qubit unitary. The test compares amplitudes to machine precision. Reaction, relaxation, and escape remain in the density-matrix solver; hardware amplitude damping is not used as a chemical parameter.

## Parameter gaps and stopping rule

Searches covered exact-name combinations for doxorubicin/adriamycin semiquinone EPR, hyperfine, pulse radiolysis, oxygen electron transfer, radical–triplet doublet/quartet selection, encounter kinetics, spin relaxation, couplings, SOD, and anthracycline ROS validation. Exact-system searches repeatedly converged on the same ESR and bulk-kinetic papers and did not produce the missing spin-resolved quantities. Further generic-quinone values would violate the requested substitution rules, so the review stops with explicit gaps.

The most valuable next measurements are time-resolved EPR/ENDOR of the encounter under a specified solvent, pH, temperature and field; spin-resolved transient kinetics; and structural/diffusion information sufficient to infer encounter residence time and distance. Until then, only sensitivity analysis is defensible.

## Claim-to-source ledger

- Akman et al., “ESR study of electron transfer reactions between gamma-irradiated pyrimidines, adriamycin and oxygen,” *Free Radical Biology & Medicine* (1990), DOI [10.1016/0891-5849(90)90027-G](https://doi.org/10.1016/0891-5849(90)90027-G). Access: PubMed abstract.
- Kalyanaraman et al., “Pulse radiolysis studies of antitumor quinones: Radical lifetimes, reactivity with oxygen, and one-electron reduction potentials,” *Archives of Biochemistry and Biophysics* (1981), DOI [10.1016/0003-9861(81)90263-0](https://doi.org/10.1016/0003-9861(81)90263-0). Access: publisher abstract.
- Land et al., “Possible intermediates in the action of adriamycin—a pulse radiolysis study,” *British Journal of Cancer* (1985), DOI [10.1038/bjc.1985.74](https://doi.org/10.1038/bjc.1985.74). Access: full text via PMC.
- Bielski et al., “Reactivity of HO2/O2- radicals in aqueous solution,” *Journal of Physical and Chemical Reference Data* 14 (1985), DOI [10.1063/1.555739](https://doi.org/10.1063/1.555739). Access: NIST/AIP reference-data review.
- Klug, Rabani & Fridovich, “A direct demonstration of the catalytic action of superoxide dismutase through the use of pulse radiolysis,” *Journal/Archives of Biochemistry and Biophysics* (1972/1973 indexing), DOI [10.1016/0003-9861(73)90636-X](https://doi.org/10.1016/0003-9861(73)90636-X). Access: publisher metadata/abstract.
- Davies & Doroshow, “Redox cycling of anthracyclines by cardiac mitochondria. II,” *Journal of Biological Chemistry* (1986), [PubMed PMID 3005279](https://pubmed.ncbi.nlm.nih.gov/3005279/). Access: PubMed abstract.
- Houée-Levin et al., “Pulse-radiolysis study of daunorubicin redox cycles,” *FEBS Letters* (1985), DOI [10.1016/0014-5793(85)80188-5](https://doi.org/10.1016/0014-5793(85)80188-5). Access: publisher abstract; context-only surrogate.

## Repository mapping

- This review: `references/doxorubicin_parameter_review.md`
- Machine-readable model and sensitivity inputs: `configs/doxorubicin_parameters.json`
- Row-level provenance: `configs/parameter_provenance.csv`
- Six-state model and kinetic solver: `spin_chemistry.py`
- CLI: `ROS.py`
- Invariants and classical/circuit validation: `tests/test_spin_chemistry.py`
