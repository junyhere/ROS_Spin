# ROS_Spin: evidence-gated doxorubicin semiquinone–oxygen chemistry

The scientific refactor separates the following stages:

1. classical semiquinone formation;
2. classical encounter association;
3. coherent doublet/quartet spin evolution in the complete 2×3 electronic space;
4. spin-selective electron transfer and encounter escape;
5. 1:1 primary superoxide formation from reacted encounters;
6. species-resolved spontaneous and SOD-mediated dismutation; and
7. compartment-dependent H2O2 loss.

Semiquinone formation and association are never forced into the density matrix.
The exact six-state classical density-matrix propagation uses a constant-system
matrix exponential and remains the reference implementation. Parameter authority is limited to
[`references/doxorubicin_parameter_review.md`](references/doxorubicin_parameter_review.md),
[`configs/doxorubicin_parameters.json`](configs/doxorubicin_parameters.json), and
[`configs/parameter_provenance.csv`](configs/parameter_provenance.csv).

## Paper reproduction

The active implementation is deliberately small: `spin_chemistry.py` contains
the six-state equations, propagation, downstream kinetics, validation, and
evidence policy; `ROS.py` runs individual cases; `sensitivity_analysis.py`
defines reusable one- and two-dimensional sensitivity sweeps;
`paper_analysis.py` generates the complete paper bundle; and `plotting.py`
contains shared noninteractive Matplotlib rendering. All encounter-level paper
outputs are labelled non-predictive.

Set up the pinned environment and run the complete tests:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m unittest discover -s tests -v
```

Generate a rapid complete bundle or the publication-grid bundle. An existing,
nonempty output directory is refused unless `--overwrite` is explicit.

```bash
.venv/bin/python paper_analysis.py --mode quick --execute-circuit \
  --output-dir results/paper-quick
.venv/bin/python paper_analysis.py --mode full --execute-circuit \
  --output-dir results/paper
```

The default formats are 300-dpi PNG plus vector PDF and SVG. Useful bounded
overrides include `--reference-rate-s`, `--samples`, `--grid-size`, `--config`,
`--provenance`, `--formats`, `--dpi`, and `--seed`. Grid bounds are illustrative
computational bounds—not measured ranges, plausible ranges, confidence
intervals, or priors.

The final paper run contains Figures 1-11 as PNG, SVG, and PDF, one source CSV
and caption file per figure, Tables 1-9 as CSV and Markdown, and reproducibility
metadata.

Run one individual encounter and one mixture benchmark:

```bash
.venv/bin/python ROS.py --mode sensitivity --allow-sensitivity \
  --sensitivity-scenario local_field_mixing_probe --initial-state unpolarized
.venv/bin/python ROS.py --mode sensitivity --allow-sensitivity \
  --sensitivity-scenario local_field_mixing_probe --initial-state mixture \
  --p-doublet 0.25
```

Run literature arithmetic and an exact condition-matched bulk calculation:

```bash
.venv/bin/python ROS.py --bulk-sq-m 1e-6 --bulk-o2-m 2e-4 \
  --bulk-rate-key doxorubicin_semiquinone_plus_oxygen_pH6 \
  --bulk-use-scope literature_arithmetic
.venv/bin/python ROS.py --bulk-sq-m 1e-6 --bulk-o2-m 2e-4 \
  --bulk-rate-key doxorubicin_semiquinone_plus_oxygen_pH6 \
  --bulk-use-scope condition_matched_prediction \
  --bulk-condition-profile land_1985_pH6
```

Circuit execution is optional and covers only coherent simulator/embedding
consistency. If Qiskit cannot load, the pipeline writes the exact dependency
failure and does not create an empty validation figure. On macOS, open all PNG
figures from a completed run with:

```bash
open results/paper/figures/*.png
```

Each run contains `data/`, `figures/`, `tables/`, `metadata/`, and a generated
`README.md`. Every figure is rendered only from its saved same-named CSV source
table; manifests record Git state, hashes, versions, invocation, grids, sample
counts, numerical tolerances, skipped features, and limitations.

## Execution modes

`evidence_backed` is the default. It refuses an encounter-level yield because
`k_D`, `k_Q`, escape/association, duration, and a prepared D/Q state are
unavailable. Bulk `M^-1 s^-1` constants are never converted into encounter
`s^-1` constants. Each bulk parameter has an explicit reaction ID and allowed
use scopes; a same-unit SOD constant is rejected. `literature_arithmetic`
evaluates the published law at supplied concentrations but is not called a
condition-matched prediction. `condition_matched_prediction` requires exact,
explicit species/protonation, pH, temperature, buffer/solvent, and oxygen
conditions. Unknown or unequal conditions are refused without invented
tolerances.

```bash
python3 ROS.py
python3 ROS.py --bulk-sq-m 1e-6 --bulk-o2-m 2e-4 \
  --bulk-rate-key doxorubicin_semiquinone_plus_oxygen_pH6 \
  --bulk-use-scope literature_arithmetic

# Separate condition-specific bulk-law reference calculation
python3 ROS.py --bulk-sq-m 1e-6 --bulk-o2-m 2e-4 \
  --bulk-rate-key doxorubicin_semiquinone_plus_oxygen_pH6 \
  --bulk-use-scope condition_matched_prediction \
  --bulk-condition-profile land_1985_pH6
```

`sensitivity` requires both an explicit opt-in and a named scenario. Every
result is labelled `NON-PREDICTIVE SENSITIVITY OUTPUT`.

```bash
python3 ROS.py --mode sensitivity \
  --sensitivity-scenario baseline_no_optional_interactions \
  --allow-sensitivity --initial-state unpolarized
python3 ROS.py --mode sensitivity \
  --sensitivity-scenario local_field_mixing_probe \
  --allow-sensitivity --initial-state doublet --execute-circuit
python3 ROS.py --mode sensitivity \
  --sensitivity-scenario local_field_mixing_probe \
  --allow-sensitivity --initial-state mixture --p-doublet 0.25
```

Unavailable interactions remain implemented and disabled in the evidence-backed
profile. Their zero-valued placeholders are algebraic omissions, never physical
evidence of zero. Named sensitivity scenarios may override them with explicitly
illustrative values.

## Equations and validation scope

The encounter Hamiltonian is

`H/ℏ = βe B·(gSQ s + gO2 S) + J s·S + s·Ddip·S + S·DO2·S + Ωlocal·s`,

with `PD = (0.5 I - s·S)/1.5`, `PQ = (s·S + I)/1.5`, and

`dρ/dt = -i[H/ℏ,ρ] - {kD PD + kQ PQ + kesc I,ρ}/2 + Lrelax(ρ)`.

Reaction and escape yields are integrated from `kD Tr(PDρ)`, `kQ Tr(PQρ)`,
and `kesc Tr(ρ)`. `PD/2` and `PQ/4` are manifold-restricted mixtures, not pure
wavefunctions; variable mixtures use `pD PD/2 + (1-pD) PQ/4`. A chemically
prepared state remains unsupported.

Unitary evolution alone leaves `I6/6` invariant. Unequal spin-selective loss can
first move the surviving ensemble away from `I6/6`, after which Hamiltonian
mixing can change subsequent yields; an initially unpolarized *reactive* model
is therefore not automatically insensitive to mixing. Isotropic exchange and
common Zeeman evolution do not mix D/Q. `local_field_proxy_rad_s` is an
illustrative local electronic-field term, not explicit nuclear hyperfine
dynamics. The local isotropic Lindblad depolarization is also illustrative.

Downstream chemistry reports HO2 and O2-minus separately and consumes two
radical equivalents per H2O2. It assumes fixed pH, rapid acid-base equilibrium,
constant SOD, constant rates, a single radical pulse, and optional first-order
H2O2 loss. The analytic radical-pool solution plus integrating-factor
quadrature preserves the distinction between reaction-event and species-loss
rates and reports the accumulated H2O2 loss in its balance.

The matrix exponential is the reference. A separately constructed adaptive
Dormand–Prince 5(4) solver validates full density matrices, D/Q populations,
reaction/escape yields, survival, and probability balance for coherent mixing,
relaxation, and selective loss. Sampling-density checks are only output
sampling/roundoff checks, not independent convergence evidence.

Declared verification tolerances are: independent-solver `atol=1e-12` and
`rtol=1e-10`; maximum absolute full-density error `2e-9`; D/Q, survival, and
yield error `3e-9`; analytic reaction/escape limits `1e-8`; probability balance
`1e-8`; circuit statevector and D/Q observable error `1e-12`; and circuit
leakage probability `1e-24`. These numerical tolerances are not physical
uncertainties.

Create the isolated environment used for circuit-inclusive verification, then
run the full suite:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install numpy==2.3.2 qiskit==2.1.1
.venv/bin/python -m unittest discover -s tests -v
```

The three-qubit validation uses Qiskit's statevector executor and one inserted
dense 8×8 `UnitaryGate` containing the 6×6 physical block; it is not a
decomposition into elementary gates. It checks Qiskit basis ordering,
physical-subspace leakage, the statevector, and D/Q observables. This is
simulator/embedding consistency, not an independent Hamiltonian derivation,
and it never validates relaxation, reaction, escape, or chemistry.

## Reproducible sensitivity outputs

Named encounter benchmarks use `ROS.py` as shown above. The small sweep is
one-factor-at-a-time rather than a large Cartesian product. It takes an
explicit positive reference rate, uses `kD=kref`, normalizes local mixing,
relaxation, and escape by `kref`, and varies `kQ/kD`. It includes zero-mixing,
fast-relaxation, and spin-independent (`kQ=kD`) controls and compares
unpolarized, doublet, quartet, and `pD=0.25/0.75` mixtures.

```bash
.venv/bin/python sensitivity_analysis.py --reference-rate-s 1e6 \
  --samples 81 --output-dir results/demo
```

This writes a CSV, a dependency-free SVG summary figure, and a JSON manifest
with commit/dirty status, source fingerprint, config/provenance hashes, resolved
grid, units, numerical settings, and limitations. Bounds are illustrative—not
measured ranges or probability priors. Outputs are populations and
per-encounter yields, not concentrations or fluxes; no arbitrary amplitude fit
or conversion to biological molarity is performed.
