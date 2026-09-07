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

## Execution modes

`evidence_backed` is the default. It can evaluate a measured bulk law
`r = k_bulk[SQ][O2]` when concentrations are supplied, but it refuses an
encounter-level yield because `k_D`, `k_Q`, escape/association, duration, and a
prepared D/Q state are unavailable. Bulk `M^-1 s^-1` constants are never
converted into encounter `s^-1` constants.

```bash
python3 ROS.py
python3 ROS.py --bulk-sq-m 1e-6 --bulk-o2-m 2e-4 \
  --bulk-rate-key doxorubicin_semiquinone_plus_oxygen_pH6
```

`sensitivity` requires both an explicit opt-in and a named scenario. Every
result is labelled `NON-PREDICTIVE SENSITIVITY OUTPUT`.

```bash
python3 ROS.py --mode sensitivity \
  --sensitivity-scenario baseline_no_optional_interactions \
  --allow-sensitivity --initial-state unpolarized
python3 ROS.py --mode sensitivity \
  --sensitivity-scenario hyperfine_mixing_probe \
  --allow-sensitivity --initial-state doublet --check-unitary-embedding
```

Unavailable interactions remain implemented and disabled in the evidence-backed
profile. Their zero-valued placeholders are algebraic omissions, never physical
evidence of zero. Named sensitivity scenarios may override them with explicitly
illustrative values.

## Equations and validation scope

The encounter Hamiltonian is

`H/ℏ = βe B·(gSQ s + gO2 S) + J s·S + s·Ddip·S + S·DO2·S + Ahf·s`,

with `PD = (0.5 I - s·S)/1.5`, `PQ = (s·S + I)/1.5`, and

`dρ/dt = -i[H/ℏ,ρ] - {kD PD + kQ PQ + kesc I,ρ}/2 + Lrelax(ρ)`.

Reaction and escape yields are integrated from `kD Tr(PDρ)`, `kQ Tr(PQρ)`,
and `kesc Tr(ρ)`. Downstream chemistry consumes two HO2/O2- radical equivalents
per H2O2 and treats spontaneous, SOD-mediated, and H2O2-loss terms separately.

Run the full tests with:

```bash
python3 -m unittest discover -s tests -v
```

The local Qiskit installation cannot load on this machine because its native
extension has the wrong architecture. Therefore the repository accurately
retains only a 6×6-to-8×8 unitary-embedding consistency check. It is not an
independently implemented three-qubit circuit and covers only closed coherent
Hamiltonian evolution—not Lindblad relaxation, reaction, escape, upstream
preparation/association, or downstream chemistry.

## Legacy data

The processed semiconductor spin-shuttling data below are retained for Git
history and reproducibility but are not used by the chemical spin model.

# Original project overview

This repository uses processed datasets for simulations regarding quantum entanglement generation in open quantum systems. The raw data originate from the Mokeev dataset "Spin-based remote entanglement generation in open quantum systems" (DOI: 10.4121/d0d1007f-c27d-491d-b7e1-cc60e38047b4). Only relevant subsets/folders that can be viewed in their DATASET.zip folder are included here.

The original dataset README on 4TU.ResearchData states that the files are released under the Creative Commons Attribution 4.0 license (CC BY 4.0). The CSV and JSON files in this repository were obtained from the Mokeev datast directly but some files were not used.

## Folder summaries

### `dataset/fig03`
* `Sim07_20240405_stochastic_field_axis_3.json` – arrays `t` and `x` defining a 201×101 grid of time and position values.
* `Sim07_20240405_stochastic_field_vals_3.csv` – stochastic magnetic field values on that grid. Each row corresponds to a fixed `x` value and contains 201 comma-separated entries for the different time points.

### `dataset/fig06`
Four CSV files (`Sim08_20240405_S1_OW_vs_FB_OU_1.csv` … `_4.csv`) giving shuttling fidelities for different ranges of the dimensionless parameter `beta`.
Columns:
* `beta` – scaled shuttling speed.
* `f_fb` – final fidelity for the feed-back protocol.
* `f_ow` – final fidelity for the one-way protocol.

### `dataset/fig08`
* `Sim06_20231019_dephasing_surface_axis.csv` – 401 rows of `(v, T0)` values describing the scan parameters.
* `Sim06_20231019_dephasing_surface_vals.csv` – a 400×400 grid of fidelities associated with those parameter pairs.

### `dataset/fig09`
CSV/JSON pairs named `realistic_tau=<value>` describing two-spin shuttling simulations for various dephasing times `tau`.
* CSV columns: `v` (velocity) and `chi_ni` (numerical infidelity).
* JSON files record the simulation settings such as magnetic field parameters and sample sizes.

### `dataset/fig10-12`
Files describing numerical error scaling for Monte Carlo (`*_MC.csv`) and numerical integration (`*_NI.csv`) approaches. Each CSV contains columns
* `N` – number of time steps used in the solver,
* `epsilon` – absolute error of the method,
* `delta` – Monte Carlo sampling error (for MC files),
* `cpu_time` – execution time in seconds,
* `ram_bytes` – memory usage.
Companion JSON files list the common parameters (`T`, `L`, `M`, `corr_t`, `corr_x`).
