| model_element | original_implementation | corrected_implementation | scientific_consequence |
| --- | --- | --- | --- |
| Physical species | Generic two-spin or semiconductor shuttling model | Doxorubicin/adriamycin semiquinone with S=1/2 and ground-state O2 with S=1 | The chemistry is tied to the stated redox pair. |
| Species multiplicity and scope | Species and multiplicities were not consistently separated | Closed-shell AQ/H2O2: 0 unpaired; semiquinone/superoxide: 1, doublet; ground-state O2: 2, triplet; singlet oxygen excluded | The retained model cannot be read as a singlet-oxygen sensitization mechanism. |
| Spin dimensions | Two qubits; dimension 4 | Spin-1/2 x spin-1 electronic space; dimension 6 | All physical electronic states are retained. |
| Manifold classification | Singlet/triplet or basis-parity labels | Doublet (dimension 2) and quartet (dimension 4) projectors | Projectors follow angular-momentum addition. |
| Initial states | Bell states or computational-basis preparations | I6/6, PD/2, PQ/4, and bounded pD mixtures | All cases are computational benchmarks; no chemical preparation is asserted. |
| Hamiltonian | Repeated CZ or identity gates used as chemical dynamics | Zeeman, exchange, dipolar, O2 ZFS, and local-field terms in the six-state space | Unavailable interactions remain disabled or sensitivity-only. |
| Relaxation | Semiconductor dephasing parameters | Local isotropic Lindblad sensitivity model with exact-system rates unavailable | No semiconductor relaxation parameter enters active chemistry. |
| Electron transfer | Mapped from state counts or parity | Integrated kD Tr(PD rho) and kQ Tr(PQ rho) loss fluxes | Reaction yield is a time-integrated kinetic quantity. |
| Encounter escape | Absent or not separated from circuit depth | Independent first-order escape channel kescape Tr(rho) | Reaction, escape, and unresolved survival close the probability balance. |
| Primary superoxide | Assigned from singlet/triplet or basis outcomes | YD + YQ from integrated electron-transfer flux | One primary radical equivalent is assigned per reacted encounter. |
| Hydrogen peroxide | Mapped directly from a spin population | Separate HO2/O2-minus speciation and dismutation with two radicals per H2O2 | Spin populations are not treated as H2O2. |
| Rate units | Bulk and encounter rates could be conflated | Bulk M^-1 s^-1 constants remain separate from encounter s^-1 rates | No unsupported dimensional conversion is made. |
| Quantum-circuit role | Circuit outcomes interpreted as chemistry | Three-qubit coherent embedding consistency check | The circuit does not validate open-system chemistry or quantum advantage. |
| Classical/circuit boundary | Classical and circuit workloads could be presented as interchangeable | Classical model includes reaction, escape, relaxation, and yields; circuit executes only one supplied statevector and dense coherent unitary | Only the coherent statevector workload is used for the bounded NumPy/Qiskit comparison. |
| Endpoint interpretation | Finite-time results risked being read as final chemical yields | Finite-time reaction and escape yields are reported with unresolved survival | No asymptotic-yield claim is made without endpoint convergence. |
| Interpretation | Risk of direct biological or cardiotoxicity extrapolation | Conditional dimensionless sensitivity analysis with evidence-gated inputs | Cellular ROS and clinical outcomes remain outside scope. |
