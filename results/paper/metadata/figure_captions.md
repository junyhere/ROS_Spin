# Generated figure captions

All encounter coordinates are normalized by the declared kref and are illustrative computational bounds, not measured ranges, confidence intervals, or priors.

## figure01_model_overview

Figure 1. Corrected model overview. Formal electron-transfer bookkeeping is AQ + e- -> AQ radical anion, followed by AQ radical anion + ground-state triplet O2 -> AQ + O2 radical anion; the modeled downstream pool uses HO2/O2 radical-anion equilibrium and two radical equivalents per H2O2. Closed-shell AQ has zero unpaired electrons, semiquinone and superoxide each have one (doublet), and ground-state O2 has two (triplet); singlet oxygen is excluded. The semiquinone protonation/formal charge is condition-dependent and is not assigned as a measured active-model input. A semiquinone-triplet-O2 encounter enters the six-state doublet/quartet stage, followed by competing conditional reaction or escape. Conceptual workflow created by the authors; not a simulation output. The diagram is non-quantitative and does not assert coherent preparation or a measured encounter.

## figure02_benchmark_trajectories

Figure 2. Benchmark encounter trajectories. Surviving doublet and quartet populations, total survival, cumulative doublet and quartet reaction yields, primary-superoxide yield, and escape yield are shown against normalized time for five benchmark mixtures. At every sampled time, YD + YQ + Yescape + Psurvival = 1 within numerical tolerance. The settings kD/kref=1, kQ/kD=0.1, kescape/kref=1, omega_local/kref=1, and gammaR/kref=gammaO/kref=0.1 are illustrative sensitivity coordinates, not measured parameters.

## figure03_initial_state_comparison

Figure 3. Initial-state comparison. Final reaction, escape, and unresolved survival outcomes are compared for the unpolarized encounter, doublet-manifold benchmark, quartet-manifold benchmark, and pD=0.25 and pD=0.75 mixtures. The doublet and quartet cases are computational limiting benchmarks rather than demonstrated chemically prepared states. Primary superoxide is the sum of integrated D and Q reaction yields, not a final spin population.

## figure04_mixing_escape_heatmaps

Figure 4. Mixing-versus-escape heatmaps. Primary-superoxide yield per encounter is calculated over normalized local D/Q mixing and escape rates for kQ/kD=0, 0.1, 1, and 10. The kQ/kD=1 panel is the spin-independent null condition. Zero and logarithmic nonzero coordinates are illustrative; no displayed range is experimentally established.

## figure05_mixing_relaxation_heatmaps

Figure 5. Mixing-versus-relaxation heatmaps. Primary-superoxide yield is calculated for an unpolarized benchmark over normalized local mixing and equal radical/O2 relaxation rates. The fixed settings are kD/kref=1, kescape/kref=1, and duration=8/kref; panels show kQ/kD=0, 0.1, 1, and 10, with equality marked as the spin-independent null. All coordinates are illustrative and non-predictive.

## figure06_reaction_selectivity

Figure 6. Reaction-selectivity sensitivity. Primary-superoxide yield is plotted against kQ/kD for unpolarized, doublet, quartet, pD=0.25, and pD=0.75 benchmarks across three mixing and escape settings. The vertical line marks kQ/kD=1 as the spin-independent null. The kQ/kD axis is a sensitivity coordinate, not a measured or chemically defensible probability distribution or physical range.

## figure07_controls

Figure 7. Controls and limiting cases. Baseline sensitivity, zero mixing, kQ=kD spin-independent reaction, fast relaxation, rapid escape, doublet benchmark, and quartet benchmark cases show primary reaction, escape, and unresolved survival. Stacked outcomes make probability accounting visible. All rates and initial-state restrictions are computational controls, not established encounter parameters or preparations.

## figure08_solver_validation

Figure 8. Independent numerical validation. A separately constructed adaptive Dormand-Prince 5(4) solver is compared with the constant-generator matrix-exponential reference for full density matrices, doublet and quartet populations, survival, reaction yields, escape yield, and probability balance. Dashed lines show declared absolute tolerances. The symmetric-log display includes exact zero without replacing it by an artificial logarithmic floor. Agreement validates numerical implementation only, not physical encounter inputs or the coherent circuit.

## figure09_downstream_kinetics

Figure 9. Downstream ROS kinetics. Total superoxide-family radical pool, HO2, O2-minus, accumulated H2O2, and accumulated H2O2 loss are shown for spontaneous dismutation, SOD-mediated dismutation, and SOD plus H2O2 loss. The calculation starts from an illustrative 10 micromolar single radical pulse at fixed pH 7.4, assumes rapid acid-base equilibrium and constant rates/SOD, and consumes two radical equivalents per H2O2. No spin population is mapped directly to H2O2, and the trajectories are not cellular predictions.

## figure10_bulk_rate_comparison

Figure 10. Literature bulk-rate comparison. Actual numerical values stored in the configuration JSON and matched provenance CSV are plotted with reported uncertainty, pH, temperature, environment, species/protonation, method, and source identifiers. Missing conditions remain visibly marked and records are not combined into a universal rate. These M^-1 s^-1 bulk total constants are not encounter-level kD or kQ in s^-1.

## figure11_circuit_validation

Figure 11. Three-qubit coherent-embedding validation. Physical six-state basis indices [0,1,2,4,5,6] are embedded in an eight-state register; indices [3,7] are unused. A normalized pseudorandom six-component statevector (seed 1729) is supplied numerically, not prepared by gates. One dense 8x8 UnitaryGate contains U=exp(-iH t) for t=1/kref; no Trotter, delay, noise, reaction, measurement, or sampling gates are executed. Statevector, numerical D/Q projector expectations, basis ordering, and unused-state leakage are compared with the direct six-state result. This validates coherent simulator and encoding consistency only; it does not independently validate H, the open-system model, chemistry, hardware performance, or quantum advantage.
