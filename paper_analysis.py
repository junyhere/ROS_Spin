"""Generate the complete reproducible ROS_Spin paper-results bundle.

All encounter inputs used by this command are bounded, illustrative sensitivity
coordinates. They are not measured ranges, priors, or biological predictions.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import statistics
import subprocess
import sys
import time
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/ros_spin_matplotlib")

import matplotlib
import numpy as np

from plotting import (
    plot_benchmark_trajectories,
    plot_bulk_rates,
    plot_circuit_validation,
    plot_controls,
    plot_downstream,
    plot_initial_state_comparison,
    plot_mixing_escape_heatmaps,
    plot_mixing_relaxation_heatmaps,
    plot_model_overview,
    plot_reaction_selectivity,
    plot_solver_validation,
    save_figure,
)
from sensitivity_analysis import (
    INITIAL_STATES,
    LIMITATIONS,
    encounter_result_row,
    run_mixing_escape_sweep,
    run_mixing_relaxation_sweep,
    run_selectivity_sweep,
    run_sweep,
)
from table_rendering import render_table_from_csv
from spin_chemistry import (
    EncounterParameters,
    P_DOUBLET,
    P_QUARTET,
    REACTION_SQ_O2_ET,
    downstream_ros_species_resolved,
    execute_three_qubit_unitary_circuit,
    hamiltonian,
    independent_circuit_capability,
    propagate_encounter_independent,
    propagate_encounter_reference,
    validate_authority_bundle,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = ROOT / "configs" / "doxorubicin_parameters.json"
DEFAULT_PROVENANCE = ROOT / "configs" / "parameter_provenance.csv"
SEED = 1729
TOLERANCES = {
    "independent_atol": 1e-12,
    "independent_rtol": 1e-10,
    "density": 2e-9,
    "observable": 3e-9,
    "balance": 1e-8,
    "circuit_statevector": 1e-12,
    "circuit_observable": 1e-12,
    "circuit_leakage": 1e-24,
    "circuit_basis_ordering": 1e-15,
}
SCIENTIFIC_LIMITATIONS = (
    "Illustrative computational sensitivity study only. Encounter preparation, "
    "kD, kQ, escape, duration, relaxation, exchange, dipolar coupling, geometry, "
    "and encounter O2 tensors are unavailable. No per-encounter yield is converted "
    "to concentration or biological flux. No measured coherence, D/Q selectivity, "
    "quantum advantage, entanglement, magnetic control, or quantitative biological "
    "ROS prediction is claimed."
)

TABLE_PRESENTATION = {
    "table01_original_versus_corrected_model": {
        "number": 1,
        "title": "Original versus corrected model",
        "note": "The comparison identifies the scientific effect of each model correction.",
        "panels": [[
            "model_element", "original_implementation", "corrected_implementation",
            "scientific_consequence",
        ]],
    },
    "table02_parameter_provenance": {
        "number": 2,
        "title": "Parameter provenance",
        "note": "Each parameter retains its stated conditions, evidence status, and permitted use.",
        "panels": [
            ["symbol", "definition", "value_or_range", "unit", "status"],
            ["symbol", "species", "charge_state", "environment", "temperature", "pH"],
            [
                "symbol", "method", "source", "uncertainty", "limitations",
                "permitted_use", "code_location",
            ],
        ],
        "panel_titles": [
            "Parameter identity", "Chemical conditions", "Evidence and permitted use",
        ],
    },
    "table03_unavailable_parameters_and_consequences": {
        "number": 3,
        "title": "Unavailable parameters and consequences",
        "note": "Unavailable inputs remain sensitivity coordinates and limit the conclusions that can be drawn.",
        "panels": [
            [
                "missing_parameter", "why_needed", "status",
                "sensitivity_coordinate_allowed", "range_status",
            ],
            [
                "missing_parameter", "allowed_treatment", "conclusion_prohibited",
                "source_of_status",
            ],
        ],
        "panel_titles": ["Evidence status", "Allowed treatment and limitation"],
    },
    "table04_baseline_results": {
        "number": 4,
        "title": "Finite-time baseline results",
        "note": (
            "The five initial-state benchmarks use the declared baseline sensitivity "
            "scenario at t = 8/kref. Unresolved survival is reported explicitly; these "
            "finite-time values are not asymptotic yields."
        ),
        "panels": [[
            "initial_state_definition", "p_doublet_initial",
            "doublet_reaction_yield_per_encounter",
            "quartet_reaction_yield_per_encounter",
            "primary_superoxide_yield_per_encounter", "escape_yield_per_encounter",
            "survival_probability", "probability_balance",
        ]],
    },
    "table05_controls_and_extrema": {
        "number": 5,
        "title": "Controls and selected grid extrema",
        "note": "The extrema describe the sampled computational grid only.",
        "panels": [
            [
                "summary_type", "source_grid", "scenario_id", "control_label",
                "initial_state_definition",
            ],
            [
                "summary_type", "mixing_over_reference",
                "radical_relaxation_over_reference", "escape_over_reference",
                "kq_over_kd", "primary_superoxide_yield_per_encounter",
                "escape_yield_per_encounter", "unresolved_probability",
            ],
        ],
        "panel_titles": ["Scenario identity", "Coordinates and results"],
    },
    "table06_independent_solver_validation": {
        "number": 6,
        "title": "Independent solver validation",
        "note": (
            "Worst discrepancies compare the matrix exponential with the adaptive "
            "Dormand-Prince solver. The normalized metric divides by the fixed unit "
            "scale for bounded dimensionless quantities, never by a near-zero result."
        ),
        "panels": [[
            "observable", "worst_absolute_error",
            "worst_unit_scale_normalized_error", "normalization_scale",
            "near_zero_policy", "declared_absolute_tolerance", "passes",
            "reference_solver", "independent_solver", "scope",
        ]],
    },
    "table07_circuit_validation": {
        "number": 7,
        "title": "Circuit validation",
        "note": "Circuit checks cover coherent statevector execution and six state to eight state embedding only.",
        "panels": [
            [
                "metric", "validation_layer", "value", "tolerance", "passes",
                "qiskit_execution_status", "status",
            ],
            [
                "metric", "seed", "qubits", "basis_order", "physical_indices",
                "duration_s", "input_state_origin", "state_preparation",
                "implementation_type",
            ],
            [
                "metric", "measurement", "post_processing", "scope", "reason",
                "does_not_validate",
            ],
        ],
        "panel_titles": [
            "Validation metrics", "Executed coherent calculation",
            "Observables and exclusions",
        ],
    },
    "table08_experimental_validation_evidence": {
        "number": 8,
        "title": "Experimental validation evidence",
        "note": "The records distinguish exact chemical systems and direct validation eligibility.",
        "panels": [
            [
                "dataset_id", "exact_chemical_or_biological_system", "observable",
                "conditions", "reported_result", "unit",
            ],
            [
                "dataset_id", "method", "source", "usable_for_direct_validation",
                "limitations",
            ],
        ],
        "panel_titles": ["System and reported result", "Method and validation use"],
    },
    "table09_requirements_traceability": {
        "number": 9,
        "title": "Requirements and critique traceability",
        "note": (
            "The critique checklist is traced item by item. Except where the critique "
            "itself attributes a request to Reviewer 1, original reviewer correspondence "
            "was unavailable; withdrawn requests are not treated as reviewer acceptance."
        ),
        "panels": [
            [
                "requirement_id", "criticism", "request_attribution", "disposition",
            ],
            [
                "requirement_id", "repository_evidence", "remaining_limitation",
            ],
            [
                "requirement_id", "exact_file", "exact_code_object",
                "figure_or_table", "test", "evidence_status",
            ],
        ],
        "panel_titles": [
            "Criticism and disposition", "Evidence and remaining limitation",
            "Implementation and verification",
        ],
    },
    "table10_runtime_memory_benchmark": {
        "number": 10,
        "title": "Bounded runtime and memory benchmark",
        "note": (
            "Timings are environment-specific medians for the retained six-state model, "
            "not a scaling or speed-advantage claim. Process peak RSS is a whole-process "
            "high-water mark; theoretical array bytes are lower-bound storage estimates."
        ),
        "panels": [
            [
                "workload_group", "method", "phase", "status", "repetitions",
                "warmup_repetitions", "median_time_ms", "minimum_time_ms",
                "maximum_time_ms",
            ],
            [
                "workload_group", "method", "matched_problem",
                "timing_boundary", "measured_process_peak_rss_mib",
                "measured_peak_rss_limitation", "theoretical_core_array_bytes",
                "theoretical_storage_scope",
            ],
            [
                "workload_group", "method", "hardware", "software",
                "transpilation_status", "hardware_execution_status",
                "sampling_status", "interpretation",
            ],
        ],
        "panel_titles": [
            "Timing", "Memory measurement and estimate", "Execution context",
        ],
    },
}

COLUMN_LABELS = {
    "pH": "pH",
    "p_doublet_initial": "Initial doublet probability",
    "kq_over_kd": "kQ divided by kD",
    "value_or_range": "Value or range",
    "charge_state": "Charge state",
    "code_location": "Code location",
    "source_of_status": "Status source",
    "exact_chemical_or_biological_system": "Exact chemical or biological system",
    "figure_or_table": "Figure or table",
    "research_topic_or_output": "Research topic or output",
    "qiskit_execution_status": "Qiskit execution status",
    "worst_absolute_error": "Worst absolute error (dimensionless)",
    "worst_unit_scale_normalized_error": "Worst unit-scale-normalized error",
    "normalization_scale": "Fixed normalization scale",
    "declared_absolute_tolerance": "Absolute tolerance (dimensionless)",
    "p_doublet_initial": "Initial doublet probability",
    "doublet_reaction_yield_per_encounter": "Doublet reaction yield (per encounter)",
    "quartet_reaction_yield_per_encounter": "Quartet reaction yield (per encounter)",
    "primary_superoxide_yield_per_encounter": "Primary superoxide (radical eq. per encounter)",
    "escape_yield_per_encounter": "Escape yield (per encounter)",
    "survival_probability": "Unresolved survival probability",
    "probability_balance": "Reaction + escape + survival",
    "duration_s": "Duration (s; illustrative conversion)",
    "median_time_ms": "Median time (ms)",
    "minimum_time_ms": "Minimum time (ms)",
    "maximum_time_ms": "Maximum time (ms)",
    "measured_process_peak_rss_mib": "Process peak RSS (MiB)",
    "theoretical_core_array_bytes": "Core-array estimate (bytes)",
}


def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"not JSON serializable: {type(value).__name__}")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _identity(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def _git(command: list[str]) -> str:
    return subprocess.run(
        ["git", *command], cwd=ROOT, check=True, text=True, capture_output=True
    ).stdout.strip()


def _source_metadata(config: Path, provenance: Path) -> dict[str, Any]:
    source_files = [
        ROOT / "README.md",
        ROOT / "ROS.py",
        ROOT / "spin_chemistry.py",
        ROOT / "sensitivity_analysis.py",
        ROOT / "plotting.py",
        ROOT / "table_rendering.py",
        Path(__file__),
        ROOT / "requirements.txt",
        ROOT / ".github" / "workflows" / "tests.yml",
        ROOT / "tests" / "test_spin_chemistry.py",
        ROOT / "tests" / "test_paper_pipeline.py",
        config,
        provenance,
        ROOT / "references" / "doxorubicin_parameter_review.md",
    ]
    digest = hashlib.sha256()
    for path in source_files:
        digest.update(_identity(path).encode("utf-8"))
        digest.update(path.read_bytes())
    dirty = _git(["status", "--porcelain"]).splitlines()
    try:
        import qiskit
        qiskit_version = qiskit.__version__
    except Exception as error:
        qiskit_version = f"unavailable: {error}"
    try:
        import scipy
        scipy_version = scipy.__version__
    except Exception as error:
        scipy_version = f"unavailable: {error}"
    return {
        "git_commit": _git(["rev-parse", "HEAD"]),
        "git_branch": _git(["rev-parse", "--abbrev-ref", "HEAD"]),
        "git_remote_origin": _git(["remote", "get-url", "origin"]),
        "dirty_tree": bool(dirty),
        "dirty_entry_count_at_start": len(dirty),
        "source_fingerprint_sha256": digest.hexdigest(),
        "source_files": [_identity(path) for path in source_files],
        "config_path": _identity(config),
        "config_sha256": _sha256(config),
        "provenance_path": _identity(provenance),
        "provenance_sha256": _sha256(provenance),
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "operating_system": platform.system(),
        "architecture": platform.machine(),
        "numpy_version": np.__version__,
        "matplotlib_version": matplotlib.__version__,
        "scipy_version": scipy_version,
        "qiskit_version": qiskit_version,
        "source_hashes_sha256": {
            _identity(path): _sha256(path) for path in source_files
        },
    }


def _fieldnames(rows: list[dict]) -> list[str]:
    names: list[str] = []
    for row in rows:
        for key in row:
            if key not in names:
                names.append(key)
    return names


def _cell(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True, default=_jsonable)
    if isinstance(value, np.generic):
        return value.item()
    return value


def write_csv(rows: list[dict], path: Path) -> Path:
    if not rows:
        raise ValueError(f"cannot write empty source table: {path}")
    fields = _fieldnames(rows)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _cell(row.get(field)) for field in fields})
    return path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"source table is empty: {path}")
    return rows


def write_markdown(
    rows: list[dict], path: Path, fields: list[str] | None = None
) -> Path:
    if not rows:
        raise ValueError(f"cannot write empty table: {path}")
    fields = fields or _fieldnames(rows)

    def clean(value: Any) -> str:
        return _display_cell(value).replace("|", "\\|").replace("\n", " ")

    lines = [
        "| " + " | ".join(fields) + " |",
        "| " + " | ".join("---" for _ in fields) + " |",
    ]
    lines.extend(
        "| " + " | ".join(clean(row.get(field, "")) for field in fields) + " |"
        for row in rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _display_cell(value: Any) -> str:
    """Return a compact human-readable value without JSON or CSV punctuation."""
    if value is None:
        return ""
    if isinstance(value, (bool, np.bool_)):
        return "Yes" if bool(value) else "No"
    if isinstance(value, (float, np.floating)):
        number = float(value)
        if not np.isfinite(number):
            return "n.a."
        if number == 0:
            return "0"
        if abs(number) >= 1e4 or abs(number) < 1e-3:
            return f"{number:.3e}"
        return f"{number:.6g}"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, dict):
        return "; ".join(
            f"{_column_label(str(key))}: {_display_cell(item)}"
            for key, item in value.items()
        )
    if isinstance(value, (list, tuple)):
        return "; ".join(_display_cell(item) for item in value)
    return str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def _column_label(field: str) -> str:
    if field in COLUMN_LABELS:
        return COLUMN_LABELS[field]
    words = field.replace("_s^-1", " rate per second").replace("_m", " molar")
    return words.replace("_", " ").capitalize()


def _presentation_fields(spec: dict[str, Any]) -> list[str]:
    fields: list[str] = []
    for panel in spec["panels"]:
        for field in panel:
            if field not in fields:
                fields.append(field)
    return fields


def _augment_rows(rows: list[dict], metadata: dict[str, Any], config: Path, provenance: Path) -> None:
    for row in rows:
        row.update({
            "git_commit": metadata["git_commit"],
            "source_fingerprint_sha256": metadata["source_fingerprint_sha256"],
            "config_path": _identity(config),
            "config_sha256": metadata["config_sha256"],
            "provenance_path": _identity(provenance),
            "provenance_sha256": metadata["provenance_sha256"],
        })


def _prepare_output(output: Path, overwrite: bool) -> dict[str, Path]:
    if output.exists() and any(output.iterdir()):
        if not overwrite:
            raise FileExistsError(
                f"output directory is not empty: {output}; use --overwrite or choose a new directory"
            )
        shutil.rmtree(output)
    directories = {
        name: output / name
        for name in ("data", "figures", "tables", "captions", "metadata")
    }
    output.mkdir(parents=True, exist_ok=True)
    for directory in directories.values():
        directory.mkdir()
    return directories


def _model_overview_rows() -> list[dict]:
    return [
        {"order": 1, "stage": "upstream", "label": "Doxorubicin\nquinone", "scientific_role": "starting redox species", "evidence_status": "identity established", "limitation": "no spin preparation implied"},
        {"order": 2, "stage": "upstream", "label": "Semiquinone\nformation", "scientific_role": "classical upstream preparation", "evidence_status": "mechanism context", "limitation": "formation fraction unavailable"},
        {"order": 3, "stage": "encounter", "label": "Semiquinone/O2\nencounter", "scientific_role": "association boundary", "evidence_status": "framework only", "limitation": "association, geometry, duration unavailable"},
        {"order": 4, "stage": "encounter", "label": "Doublet/quartet\ndynamics", "scientific_role": "six-state density matrix", "evidence_status": "exact spin algebra", "limitation": "coherence and prepared D/Q state not established"},
        {"order": 5, "stage": "encounter", "label": "Reaction\nor escape", "scientific_role": "spin-selective loss and escape", "evidence_status": "sensitivity variables", "limitation": "kD, kQ and escape unavailable"},
        {"order": 6, "stage": "downstream", "label": "Primary\nsuperoxide", "scientific_role": "one radical equivalent per reacted encounter", "evidence_status": "stoichiometric mapping", "limitation": "per-encounter yield is not concentration"},
        {"order": 7, "stage": "downstream", "label": "Spontaneous/SOD\ndismutation", "scientific_role": "two radicals per H2O2", "evidence_status": "condition-specific aqueous constants", "limitation": "fixed pH and constant SOD"},
        {"order": 8, "stage": "downstream", "label": "H2O2\nformation/loss", "scientific_role": "separate downstream product and sink", "evidence_status": "bounded kinetic model", "limitation": "loss is illustrative and compartment-specific"},
    ]


def _model_comparison_rows() -> list[dict]:
    """Structured source for Table 1; rendered values come from its saved CSV."""
    return [
        {"model_element": "Physical species", "original_implementation": "Generic two-spin or semiconductor shuttling model", "corrected_implementation": "Doxorubicin/adriamycin semiquinone with S=1/2 and ground-state O2 with S=1", "scientific_consequence": "The chemistry is tied to the stated redox pair."},
        {"model_element": "Species multiplicity and scope", "original_implementation": "Species and multiplicities were not consistently separated", "corrected_implementation": "Closed-shell AQ/H2O2: 0 unpaired; semiquinone/superoxide: 1, doublet; ground-state O2: 2, triplet; singlet oxygen excluded", "scientific_consequence": "The retained model cannot be read as a singlet-oxygen sensitization mechanism."},
        {"model_element": "Spin dimensions", "original_implementation": "Two qubits; dimension 4", "corrected_implementation": "Spin-1/2 x spin-1 electronic space; dimension 6", "scientific_consequence": "All physical electronic states are retained."},
        {"model_element": "Manifold classification", "original_implementation": "Singlet/triplet or basis-parity labels", "corrected_implementation": "Doublet (dimension 2) and quartet (dimension 4) projectors", "scientific_consequence": "Projectors follow angular-momentum addition."},
        {"model_element": "Initial states", "original_implementation": "Bell states or computational-basis preparations", "corrected_implementation": "I6/6, PD/2, PQ/4, and bounded pD mixtures", "scientific_consequence": "All cases are computational benchmarks; no chemical preparation is asserted."},
        {"model_element": "Hamiltonian", "original_implementation": "Repeated CZ or identity gates used as chemical dynamics", "corrected_implementation": "Zeeman, exchange, dipolar, O2 ZFS, and local-field terms in the six-state space", "scientific_consequence": "Unavailable interactions remain disabled or sensitivity-only."},
        {"model_element": "Relaxation", "original_implementation": "Semiconductor dephasing parameters", "corrected_implementation": "Local isotropic Lindblad sensitivity model with exact-system rates unavailable", "scientific_consequence": "No semiconductor relaxation parameter enters active chemistry."},
        {"model_element": "Electron transfer", "original_implementation": "Mapped from state counts or parity", "corrected_implementation": "Integrated kD Tr(PD rho) and kQ Tr(PQ rho) loss fluxes", "scientific_consequence": "Reaction yield is a time-integrated kinetic quantity."},
        {"model_element": "Encounter escape", "original_implementation": "Absent or not separated from circuit depth", "corrected_implementation": "Independent first-order escape channel kescape Tr(rho)", "scientific_consequence": "Reaction, escape, and unresolved survival close the probability balance."},
        {"model_element": "Primary superoxide", "original_implementation": "Assigned from singlet/triplet or basis outcomes", "corrected_implementation": "YD + YQ from integrated electron-transfer flux", "scientific_consequence": "One primary radical equivalent is assigned per reacted encounter."},
        {"model_element": "Hydrogen peroxide", "original_implementation": "Mapped directly from a spin population", "corrected_implementation": "Separate HO2/O2-minus speciation and dismutation with two radicals per H2O2", "scientific_consequence": "Spin populations are not treated as H2O2."},
        {"model_element": "Rate units", "original_implementation": "Bulk and encounter rates could be conflated", "corrected_implementation": "Bulk M^-1 s^-1 constants remain separate from encounter s^-1 rates", "scientific_consequence": "No unsupported dimensional conversion is made."},
        {"model_element": "Quantum-circuit role", "original_implementation": "Circuit outcomes interpreted as chemistry", "corrected_implementation": "Three-qubit coherent embedding consistency check", "scientific_consequence": "The circuit does not validate open-system chemistry or quantum advantage."},
        {"model_element": "Classical/circuit boundary", "original_implementation": "Classical and circuit workloads could be presented as interchangeable", "corrected_implementation": "Classical model includes reaction, escape, relaxation, and yields; circuit executes only one supplied statevector and dense coherent unitary", "scientific_consequence": "Only the coherent statevector workload is used for the bounded NumPy/Qiskit comparison."},
        {"model_element": "Endpoint interpretation", "original_implementation": "Finite-time results risked being read as final chemical yields", "corrected_implementation": "Finite-time reaction and escape yields are reported with unresolved survival", "scientific_consequence": "No asymptotic-yield claim is made without endpoint convergence."},
        {"model_element": "Interpretation", "original_implementation": "Risk of direct biological or cardiotoxicity extrapolation", "corrected_implementation": "Conditional dimensionless sensitivity analysis with evidence-gated inputs", "scientific_consequence": "Cellular ROS and clinical outcomes remain outside scope."},
    ]


def _benchmark_data(reference_rate: float, samples: int) -> tuple[list[dict], list[dict]]:
    trajectory_rows: list[dict] = []
    final_rows: list[dict] = []
    parameters = EncounterParameters(
        local_field_proxy_rad_s=(reference_rate, 0.0, 0.0),
        radical_relaxation_s=0.1 * reference_rate,
        oxygen_relaxation_s=0.1 * reference_rate,
        k_doublet_s=reference_rate,
        k_quartet_s=0.1 * reference_rate,
        k_escape_s=reference_rate,
    )
    duration_s = 8 / reference_rate
    for state, p_doublet, label in INITIAL_STATES:
        result = propagate_encounter_reference(parameters, duration_s, state, samples, p_doublet)
        for index, time_s in enumerate(result["time_s"]):
            trajectory_rows.append({
                "scenario_id": "illustrative_benchmark",
                "initial_state": state,
                "initial_state_definition": label,
                "p_doublet_requested": "" if p_doublet is None else p_doublet,
                "p_doublet_initial": result["p_doublet_initial"],
                "time_s": float(time_s),
                "normalized_time": float(time_s * reference_rate),
                "surviving_doublet_population": float(result["p_doublet"][index]),
                "surviving_quartet_population": float(result["p_quartet"][index]),
                "total_survival": float(result["survival"][index]),
                "cumulative_doublet_reaction_yield": float(result["cumulative_doublet_reaction_yield"][index]),
                "cumulative_quartet_reaction_yield": float(result["cumulative_quartet_reaction_yield"][index]),
                "cumulative_escape_yield": float(result["cumulative_escape_yield"][index]),
                "cumulative_primary_superoxide_yield": float(result["cumulative_primary_superoxide_yield"][index]),
                "reference_rate_s^-1": reference_rate,
                "duration_s": duration_s,
                "mixing_over_reference": 1.0,
                "local_mixing_proxy_rad_s": reference_rate,
                "radical_relaxation_over_reference": 0.1,
                "radical_relaxation_s^-1": 0.1 * reference_rate,
                "oxygen_relaxation_over_reference": 0.1,
                "oxygen_relaxation_s^-1": 0.1 * reference_rate,
                "escape_over_reference": 1.0,
                "k_escape_s^-1": reference_rate,
                "kd_over_reference": 1.0,
                "k_doublet_s^-1": reference_rate,
                "kq_over_kd": 0.1,
                "k_quartet_s^-1": 0.1 * reference_rate,
                "solver": result["numerical_method"],
                "solver_samples": samples,
                "non_predictive": True,
                "output_class": "dimensionless populations/per-encounter cumulative yields",
                "units": "dimensionless except time_s and reference_rate_s^-1",
                "limitations": LIMITATIONS,
            })
        final_rows.append(encounter_result_row(
            scenario_id="illustrative_benchmark", reference_rate_s=reference_rate,
            mixing_over_reference=1.0, radical_relaxation_over_reference=0.1,
            oxygen_relaxation_over_reference=0.1, escape_over_reference=1.0,
            kq_over_kd=0.1, initial_state=state, p_doublet=p_doublet, samples=2,
        ))
    return trajectory_rows, final_rows


def _controls(reference_rate: float) -> list[dict]:
    definitions = [
        (
            "baseline_sensitivity", "Baseline sensitivity",
            EncounterParameters(
                local_field_proxy_rad_s=(reference_rate, 0, 0),
                radical_relaxation_s=.1 * reference_rate,
                oxygen_relaxation_s=.1 * reference_rate,
                k_doublet_s=reference_rate, k_quartet_s=.1 * reference_rate,
                k_escape_s=reference_rate,
            ),
            "Illustrative reference case; no parameter is an exact-system estimate",
            "unpolarized", None,
        ),
        (
            "zero_mixing", "Zero mixing",
            EncounterParameters(
                radical_relaxation_s=.1 * reference_rate,
                oxygen_relaxation_s=.1 * reference_rate,
                k_doublet_s=reference_rate, k_quartet_s=.1 * reference_rate,
                k_escape_s=reference_rate,
            ),
            "No coherent D/Q mixing term", "unpolarized", None,
        ),
        (
            "spin_independent_reaction", "Spin-independent reaction",
            EncounterParameters(
                local_field_proxy_rad_s=(reference_rate, 0, 0),
                radical_relaxation_s=.1 * reference_rate,
                oxygen_relaxation_s=.1 * reference_rate,
                k_doublet_s=reference_rate, k_quartet_s=reference_rate,
                k_escape_s=reference_rate,
            ),
            "kQ=kD null; mixing cannot change total reaction probability",
            "unpolarized", None,
        ),
        (
            "fast_relaxation", "Fast relaxation",
            EncounterParameters(
                local_field_proxy_rad_s=(reference_rate, 0, 0),
                radical_relaxation_s=100 * reference_rate,
                oxygen_relaxation_s=100 * reference_rate,
                k_doublet_s=reference_rate, k_quartet_s=.1 * reference_rate,
                k_escape_s=reference_rate,
            ),
            "Illustrative local isotropic depolarization", "unpolarized", None,
        ),
        (
            "rapid_escape", "Rapid escape",
            EncounterParameters(
                local_field_proxy_rad_s=(reference_rate, 0, 0),
                radical_relaxation_s=.1 * reference_rate,
                oxygen_relaxation_s=.1 * reference_rate,
                k_doublet_s=reference_rate, k_quartet_s=.1 * reference_rate,
                k_escape_s=100 * reference_rate,
            ),
            "Escape dominates before most modeled electron transfer", "unpolarized", None,
        ),
        (
            "doublet_benchmark", "Doublet benchmark",
            EncounterParameters(
                local_field_proxy_rad_s=(reference_rate, 0, 0),
                radical_relaxation_s=.1 * reference_rate,
                oxygen_relaxation_s=.1 * reference_rate,
                k_doublet_s=reference_rate, k_quartet_s=.1 * reference_rate,
                k_escape_s=reference_rate,
            ),
            "Computational limiting benchmark; not a demonstrated prepared state",
            "doublet", None,
        ),
        (
            "quartet_benchmark", "Quartet benchmark",
            EncounterParameters(
                local_field_proxy_rad_s=(reference_rate, 0, 0),
                radical_relaxation_s=.1 * reference_rate,
                oxygen_relaxation_s=.1 * reference_rate,
                k_doublet_s=reference_rate, k_quartet_s=.1 * reference_rate,
                k_escape_s=reference_rate,
            ),
            "Computational limiting benchmark; not a demonstrated prepared state",
            "quartet", None,
        ),
    ]
    rows = []
    for scenario_id, label, parameters, meaning, state, p_doublet in definitions:
        result = propagate_encounter_reference(
            parameters, 8/reference_rate, state, 2, p_doublet
        )
        matrix = hamiltonian(parameters)
        commutator = float(
            np.linalg.norm(matrix @ P_DOUBLET - P_DOUBLET @ matrix) / reference_rate
        )
        rows.append({
            "scenario_id": scenario_id, "control_label": label,
            "initial_state": state, "p_doublet_initial": result["p_doublet_initial"],
            "initial_state_definition": {
                "unpolarized": "unpolarized I6/6",
                "doublet": "doublet-manifold mixture PD/2",
                "quartet": "quartet-manifold mixture PQ/4",
            }[state],
            "reference_rate_s^-1": reference_rate, "duration_over_reference": 8.0,
            "duration_s": 8/reference_rate, "field_t": parameters.field_t,
            "g_radical": parameters.g_radical, "g_oxygen": parameters.g_oxygen,
            "exchange_rad_s": parameters.exchange_rad_s,
            "local_mixing_proxy_rad_s": parameters.local_field_proxy_rad_s,
            "mixing_over_reference": float(
                np.linalg.norm(parameters.local_field_proxy_rad_s) / reference_rate
            ),
            "radical_relaxation_s^-1": parameters.radical_relaxation_s,
            "radical_relaxation_over_reference": (
                parameters.radical_relaxation_s / reference_rate
            ),
            "oxygen_relaxation_s^-1": parameters.oxygen_relaxation_s,
            "k_doublet_s^-1": parameters.k_doublet_s,
            "k_quartet_s^-1": parameters.k_quartet_s,
            "k_escape_s^-1": parameters.k_escape_s,
            "escape_over_reference": parameters.k_escape_s / reference_rate,
            "kq_over_kd": (
                parameters.k_quartet_s / parameters.k_doublet_s
                if parameters.k_doublet_s else None
            ),
            "dq_commutator_frobenius": commutator,
            "dq_commutator_normalization": "divided by k_ref",
            "can_coherently_mix_dq": commutator > 1e-12,
            "control_meaning": meaning,
            "final_doublet_population": float(result["p_doublet"][-1]),
            "final_quartet_population": float(result["p_quartet"][-1]),
            "doublet_reaction_yield_per_encounter": result["doublet_reaction_yield"],
            "quartet_reaction_yield_per_encounter": result["quartet_reaction_yield"],
            "primary_superoxide_yield_per_encounter": result["primary_superoxide_yield"],
            "escape_yield_per_encounter": result["escape_yield"],
            "unresolved_probability": result["unresolved_probability"],
            "probability_balance": result["probability_balance"],
            "probability_balance_error": result["probability_balance_error"],
            "solver": result["numerical_method"], "solver_samples": 2,
            "non_predictive": True,
            "output_class": "dimensionless control probabilities/per-encounter yields",
            "units": "dimensionless except explicit rate/frequency fields",
            "limitations": LIMITATIONS,
        })
    return rows


def _solver_validation(reference_rate: float, samples: int) -> tuple[list[dict], list[dict]]:
    parameters = EncounterParameters(
        field_t=(0, 0, 2e-6), g_radical=2.0035, g_oxygen=2.0023,
        exchange_rad_s=0.4 * reference_rate,
        dipolar_rad_s=0.2 * reference_rate, dipolar_axis=(1, 1, 2),
        oxygen_zfs_d_rad_s=0.15 * reference_rate,
        oxygen_zfs_e_rad_s=0.07 * reference_rate,
        local_field_proxy_rad_s=(1.7 * reference_rate, 0.3 * reference_rate, 0),
        radical_relaxation_s=0.3 * reference_rate,
        oxygen_relaxation_s=0.2 * reference_rate,
        k_doublet_s=reference_rate,
        k_quartet_s=0.23 * reference_rate,
        k_escape_s=0.7 * reference_rate,
    )
    duration = 6/reference_rate
    ref = propagate_encounter_reference(parameters, duration, "mixture", samples, 0.75)
    ind = propagate_encounter_independent(
        parameters, duration, "mixture", samples, 0.75,
        atol=TOLERANCES["independent_atol"],
        rtol=TOLERANCES["independent_rtol"],
    )
    rows = []
    for index, time_s in enumerate(ref["time_s"]):
        ref_balance = (
            ref["cumulative_primary_superoxide_yield"][index]
            + ref["cumulative_escape_yield"][index] + ref["survival"][index]
        )
        ind_balance = (
            ind["cumulative_primary_superoxide_yield"][index]
            + ind["cumulative_escape_yield"][index] + ind["survival"][index]
        )
        rows.append({
            "scenario_id": "independent_solver_validation",
            "normalized_time": float(time_s * reference_rate),
            "time_s": float(time_s),
            "max_density_matrix_abs_error": float(np.max(np.abs(ref["density_matrices"][index] - ind["density_matrices"][index]))),
            "doublet_population_abs_difference": float(abs(ref["p_doublet"][index] - ind["p_doublet"][index])),
            "quartet_population_abs_difference": float(abs(ref["p_quartet"][index] - ind["p_quartet"][index])),
            "survival_abs_difference": float(abs(ref["survival"][index] - ind["survival"][index])),
            "doublet_reaction_yield_abs_difference": float(abs(ref["cumulative_doublet_reaction_yield"][index] - ind["cumulative_doublet_reaction_yield"][index])),
            "quartet_reaction_yield_abs_difference": float(abs(ref["cumulative_quartet_reaction_yield"][index] - ind["cumulative_quartet_reaction_yield"][index])),
            "primary_superoxide_yield_abs_difference": float(abs(ref["cumulative_primary_superoxide_yield"][index] - ind["cumulative_primary_superoxide_yield"][index])),
            "escape_yield_abs_difference": float(abs(ref["cumulative_escape_yield"][index] - ind["cumulative_escape_yield"][index])),
            "reference_probability_accounting_error": float(abs(ref_balance - 1)),
            "independent_probability_accounting_error": float(abs(ind_balance - 1)),
            "reference_solver": ref["numerical_method"],
            "independent_solver": ind["numerical_method"],
            "independent_atol": TOLERANCES["independent_atol"],
            "independent_rtol": TOLERANCES["independent_rtol"],
            "independent_accepted_steps": ind["accepted_steps"],
            "reference_rate_s^-1": reference_rate,
            "initial_state_definition": "mixture(pD=0.75)",
            "mixing_over_reference": 1.7,
            "local_mixing_proxy_rad_s": 1.7 * reference_rate,
            "radical_relaxation_over_reference": 0.3,
            "radical_relaxation_s^-1": 0.3 * reference_rate,
            "oxygen_relaxation_over_reference": 0.2,
            "oxygen_relaxation_s^-1": 0.2 * reference_rate,
            "escape_over_reference": 0.7,
            "k_escape_s^-1": 0.7 * reference_rate,
            "kd_over_reference": 1.0,
            "k_doublet_s^-1": reference_rate,
            "kq_over_kd": 0.23,
            "k_quartet_s^-1": 0.23 * reference_rate,
            "non_predictive": True,
            "output_class": "absolute numerical error",
            "units": "dimensionless except time_s and reference_rate_s^-1",
            "limitations": "Numerical agreement only; not validation of physical encounter inputs.",
        })
    metric_specs = {
        "density_matrix": ("max_density_matrix_abs_error", TOLERANCES["density"]),
        "doublet_population": ("doublet_population_abs_difference", TOLERANCES["observable"]),
        "quartet_population": ("quartet_population_abs_difference", TOLERANCES["observable"]),
        "survival": ("survival_abs_difference", TOLERANCES["observable"]),
        "doublet_reaction_yield": ("doublet_reaction_yield_abs_difference", TOLERANCES["observable"]),
        "quartet_reaction_yield": ("quartet_reaction_yield_abs_difference", TOLERANCES["observable"]),
        "primary_superoxide_yield": ("primary_superoxide_yield_abs_difference", TOLERANCES["observable"]),
        "escape_yield": ("escape_yield_abs_difference", TOLERANCES["observable"]),
        "reference_probability_accounting": ("reference_probability_accounting_error", TOLERANCES["balance"]),
        "independent_probability_accounting": ("independent_probability_accounting_error", TOLERANCES["balance"]),
    }
    summary = []
    for observable, (field, tolerance) in metric_specs.items():
        absolute = max(float(row[field]) for row in rows)
        # Every reported quantity is dimensionless and bounded on a unit scale.
        # This deliberately avoids a data-derived denominator near zero.
        scale = 1.0
        summary.append({
            "observable": observable,
            "worst_absolute_error": absolute,
            "worst_unit_scale_normalized_error": absolute / scale,
            "normalization_scale": scale,
            "near_zero_policy": "fixed unit scale; no data-derived denominator",
            "declared_absolute_tolerance": tolerance,
            "passes": absolute <= tolerance,
            "reference_solver": ref["numerical_method"],
            "independent_solver": ind["numerical_method"],
            "scope": "numerical implementation agreement, not physical validation",
        })
    return rows, summary


def _process_peak_rss_mib() -> float | None:
    """Return whole-process peak RSS in MiB when the platform exposes it."""
    try:
        import resource

        peak = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    except (ImportError, OSError, ValueError):
        return None
    # macOS reports bytes; Linux and the GitHub Actions runner report KiB.
    peak_bytes = peak if sys.platform == "darwin" else peak * 1024.0
    return peak_bytes / (1024.0**2)


def _timed_operation(operation, repetitions: int = 7, warmups: int = 2) -> dict[str, Any]:
    for _ in range(warmups):
        operation()
    elapsed_ms = []
    for _ in range(repetitions):
        start = time.perf_counter_ns()
        operation()
        elapsed_ms.append((time.perf_counter_ns() - start) / 1e6)
    return {
        "repetitions": repetitions,
        "warmup_repetitions": warmups,
        "median_time_ms": statistics.median(elapsed_ms),
        "minimum_time_ms": min(elapsed_ms),
        "maximum_time_ms": max(elapsed_ms),
        "measured_process_peak_rss_mib": _process_peak_rss_mib(),
    }


def _runtime_memory_benchmark(
    reference_rate: float, seed: int, execute_circuit: bool, metadata: dict[str, Any]
) -> list[dict]:
    """Benchmark only matched retained six-state calculations, never scaling."""
    hardware = (
        f"{platform.system()} {platform.machine()}; logical CPUs={os.cpu_count()}; "
        f"{platform.platform()}"
    )
    software = (
        f"Python {platform.python_version()}; NumPy {np.__version__}; "
        f"Qiskit {metadata['qiskit_version']}"
    )
    common = {
        "hardware": hardware,
        "software": software,
        "measured_peak_rss_limitation": (
            "whole Python-process high-water mark; not incremental or attributable "
            "to one operation"
        ),
        "transpilation_status": "not applicable unless stated",
        "hardware_execution_status": "not used",
        "sampling_status": "not used; exact statevectors and expectations",
        "interpretation": (
            "bounded reproducibility benchmark only; no speed, efficiency, scaling, "
            "or quantum-advantage inference"
        ),
    }

    coherent_parameters = EncounterParameters(
        field_t=(0, 0, 1e-4), g_radical=2.0035, g_oxygen=2.0023,
        exchange_rad_s=2 * reference_rate, dipolar_rad_s=0.1 * reference_rate,
        local_field_proxy_rad_s=(0.3 * reference_rate, 0, 0),
    )
    duration = 1 / reference_rate
    energies, vectors = np.linalg.eigh(hamiltonian(coherent_parameters))
    unitary6 = (vectors * np.exp(-1j * energies * duration)) @ vectors.conj().T
    rng = np.random.default_rng(seed)
    state6 = rng.normal(size=6) + 1j * rng.normal(size=6)
    state6 /= np.linalg.norm(state6)
    physical = np.array([0, 1, 2, 4, 5, 6])
    unused = np.array([3, 7])
    unitary8 = np.eye(8, dtype=complex)
    unitary8[np.ix_(physical, physical)] = unitary6
    state8 = np.zeros(8, dtype=complex)
    state8[physical] = state6
    pd8 = np.zeros((8, 8), dtype=complex)
    pq8 = np.zeros((8, 8), dtype=complex)
    pd8[np.ix_(physical, physical)] = P_DOUBLET
    pq8[np.ix_(physical, physical)] = P_QUARTET

    rows: list[dict[str, Any]] = []

    def add_timed(
        workload_group: str,
        method: str,
        phase: str,
        matched_problem: str,
        timing_boundary: str,
        theoretical_bytes: int,
        theoretical_scope: str,
        operation,
        **overrides: Any,
    ) -> None:
        rows.append({
            "workload_group": workload_group,
            "method": method,
            "phase": phase,
            "status": "measured",
            "matched_problem": matched_problem,
            "timing_boundary": timing_boundary,
            "theoretical_core_array_bytes": theoretical_bytes,
            "theoretical_storage_scope": theoretical_scope,
            **common,
            **_timed_operation(operation),
            **overrides,
        })

    def add_not_applicable(method: str, phase: str, reason: str) -> None:
        rows.append({
            "workload_group": "matched coherent statevector",
            "method": method,
            "phase": phase,
            "status": "not applicable",
            "repetitions": 0,
            "warmup_repetitions": 0,
            "median_time_ms": None,
            "minimum_time_ms": None,
            "maximum_time_ms": None,
            "matched_problem": (
                "same supplied normalized six-state vector, dense U=exp(-iHt), "
                f"t=1/kref={duration:.6g} s"
            ),
            "timing_boundary": reason,
            "measured_process_peak_rss_mib": None,
            "theoretical_core_array_bytes": 0,
            "theoretical_storage_scope": "not applicable",
            **common,
        })

    coherent_problem = (
        "same supplied normalized six-state vector, dense U=exp(-iHt), "
        f"t=1/kref={duration:.6g} s"
    )
    add_timed(
        "matched coherent statevector",
        "NumPy dense six-state propagation",
        "execution",
        coherent_problem,
        "precomputed 6x6 unitary multiplied by supplied 6-element statevector",
        unitary6.nbytes + 2 * state6.nbytes,
        "one 6x6 complex unitary plus input and output statevectors",
        lambda: unitary6 @ state6,
        transpilation_status="not applicable",
    )

    capability = independent_circuit_capability()
    if execute_circuit and capability["available"]:
        from qiskit import QuantumCircuit
        from qiskit.circuit.library import UnitaryGate
        from qiskit.quantum_info import Statevector

        def construct_circuit():
            circuit = QuantumCircuit(3)
            circuit.append(UnitaryGate(unitary8, label="exp(-iH6t) embedded"), [0, 1, 2])
            return circuit

        add_timed(
            "matched coherent statevector", "Qiskit Statevector dense unitary",
            "circuit construction", coherent_problem,
            "construct QuantumCircuit and append precomputed 8x8 UnitaryGate",
            unitary8.nbytes,
            "embedded 8x8 complex unitary only; excludes Qiskit object overhead",
            construct_circuit,
        )
        circuit = construct_circuit()
        add_timed(
            "matched coherent statevector", "Qiskit Statevector dense unitary",
            "execution", coherent_problem,
            "Statevector(input).evolve(preconstructed circuit); no transpilation or shots",
            unitary8.nbytes + 2 * state8.nbytes,
            "embedded unitary plus input and output 8-element statevectors",
            lambda: Statevector(state8).evolve(circuit).data,
        )
        simulated = Statevector(state8).evolve(circuit).data
        add_timed(
            "matched coherent statevector", "Qiskit Statevector dense unitary",
            "post-processing", coherent_problem,
            "compute D/Q expectations and unused-state leakage from exact statevector",
            simulated.nbytes + pd8.nbytes + pq8.nbytes,
            "output statevector plus embedded D/Q projector arrays",
            lambda: (
                np.vdot(simulated, pd8 @ simulated),
                np.vdot(simulated, pq8 @ simulated),
                np.sum(np.abs(simulated[unused]) ** 2),
            ),
        )
        add_not_applicable(
            "Qiskit Statevector dense unitary", "transpilation",
            "not applicable: Statevector evolves the inserted dense unitary directly",
        )
        add_not_applicable(
            "Qiskit Statevector dense unitary", "sampling",
            "not applicable: no finite-shot measurement is executed",
        )
        add_not_applicable(
            "quantum hardware", "hardware execution",
            "not applicable: no quantum hardware is used",
        )
    else:
        reason = (
            "not requested by pipeline"
            if not execute_circuit else str(capability["reason"])
        )
        for phase in (
            "circuit construction", "execution", "post-processing",
            "transpilation", "sampling", "hardware execution",
        ):
            add_not_applicable("Qiskit Statevector dense unitary", phase, reason)

    open_parameters = EncounterParameters(
        field_t=(0, 0, 2e-6), g_radical=2.0035, g_oxygen=2.0023,
        exchange_rad_s=0.4 * reference_rate,
        dipolar_rad_s=0.2 * reference_rate, dipolar_axis=(1, 1, 2),
        oxygen_zfs_d_rad_s=0.15 * reference_rate,
        oxygen_zfs_e_rad_s=0.07 * reference_rate,
        local_field_proxy_rad_s=(1.7 * reference_rate, 0.3 * reference_rate, 0),
        radical_relaxation_s=0.3 * reference_rate,
        oxygen_relaxation_s=0.2 * reference_rate,
        k_doublet_s=reference_rate,
        k_quartet_s=0.23 * reference_rate,
        k_escape_s=0.7 * reference_rate,
    )
    open_duration = 4 / reference_rate
    open_samples = 31
    open_problem = (
        "same 39-component open-system IVP; mixture pD=0.75; matched H, "
        f"reaction, relaxation, escape, {open_samples} output times, t=4/kref"
    )
    complex_bytes = np.dtype(complex).itemsize
    add_timed(
        "matched classical open system",
        "constant-generator matrix exponential",
        "execution",
        open_problem,
        "complete solver call including generator construction, propagation, and outputs",
        (39 * 39 + open_samples * 39) * complex_bytes,
        "39x39 generator plus 39-component stored trajectory; lower bound",
        lambda: propagate_encounter_reference(
            open_parameters, open_duration, "mixture", open_samples, 0.75
        ),
        transpilation_status="not applicable",
    )
    add_timed(
        "matched classical open system",
        "adaptive Dormand-Prince RK5(4)",
        "execution",
        open_problem,
        "complete solver call including adaptive integration and output assembly",
        (open_samples * 39 + 8 * 39) * complex_bytes,
        "stored 39-component trajectory plus eight stage/state vectors; lower bound",
        lambda: propagate_encounter_independent(
            open_parameters, open_duration, "mixture", open_samples, 0.75,
            atol=TOLERANCES["independent_atol"],
            rtol=TOLERANCES["independent_rtol"],
        ),
        transpilation_status="not applicable",
    )
    return rows


def _downstream_data(authority: dict, samples: int) -> list[dict]:
    evidence = authority["reference_evidence"]["downstream"]
    pka = evidence["HO2_pKa"]["value"]
    k_hh = evidence["k_HO2_HO2"]["value"]
    k_ha = evidence["k_HO2_O2minus"]["value"]
    k_sod = evidence["k_CuZnSOD_overall"]["value"]
    superoxide0 = 1e-5
    duration = 0.1
    definitions = [
        ("spontaneous_only", "Spontaneous only", None, 0.0, None),
        ("sod_mediated", "SOD-mediated", k_sod, 1e-8, None),
        ("sod_plus_h2o2_loss", "SOD + H2O2 loss", k_sod, 1e-8, 20.0),
    ]
    rows = []
    for scenario, label, sod_rate, sod_m, loss in definitions:
        result = downstream_ros_species_resolved(
            superoxide0, duration, 7.4, pka, k_hh, k_ha,
            k_sod_m_inv_s=sod_rate, sod_m=sod_m,
            k_h2o2_loss_s=loss, samples=samples,
        )
        for index, time_s in enumerate(result["time_s"]):
            rows.append({
                "scenario_id": scenario, "scenario_label": label,
                "time_s": float(time_s), "initial_radical_pulse_m": superoxide0,
                "radical_pool_m": float(result["radical_pool_m"][index]),
                "ho2_m": float(result["ho2_m"][index]),
                "o2minus_m": float(result["o2minus_m"][index]),
                "hydrogen_peroxide_m": float(result["hydrogen_peroxide_m"][index]),
                "accumulated_hydrogen_peroxide_loss_m": float(result["accumulated_hydrogen_peroxide_loss_m"][index]),
                "radical_equivalent_balance_m": float(result["radical_equivalent_balance_m"][index]),
                "radical_equivalent_balance_error_m": float(abs(result["radical_equivalent_balance_m"][index] - superoxide0)),
                "pH": 7.4, "pKa_HO2": pka,
                "k_HO2_HO2_m^-1_s^-1": k_hh,
                "k_HO2_O2minus_m^-1_s^-1": k_ha,
                "k_SOD_m^-1_s^-1": "" if sod_rate is None else sod_rate,
                "SOD_m": sod_m,
                "k_H2O2_loss_s^-1": "" if loss is None else loss,
                "solver": result["numerical_method"], "solver_samples": samples,
                "assumptions": result["approximations"], "non_predictive": True,
                "output_class": "illustrative single-pulse molar concentration",
                "units": "M except time_s and rate fields",
                "limitations": "Fixed pH, rapid equilibrium, constant SOD and rates, single pulse; not continuously driven cellular ROS.",
            })
    return rows


def _bulk_rows(authority: dict, provenance_path: Path) -> list[dict]:
    with provenance_path.open(newline="", encoding="utf-8") as handle:
        provenance_by_location = {
            row["code_location"]: row for row in csv.DictReader(handle)
        }
    rows = []
    for parameter_id, record in authority["parameter_records"].items():
        if record.get("reaction_id") != REACTION_SQ_O2_ET:
            continue
        code_location = f"parameter_records.{parameter_id}"
        provenance = provenance_by_location.get(code_location)
        if provenance is None:
            raise ValueError(f"bulk parameter lacks provenance row: {code_location}")
        if provenance.get("unit") != record.get("unit"):
            raise ValueError(f"bulk parameter/provenance unit mismatch: {parameter_id}")
        conditions = record.get("conditions", {})
        rows.append({
            "parameter_id": parameter_id, "reaction_id": record["reaction_id"],
            "provenance_symbol": provenance["symbol"],
            "provenance_code_location": provenance["code_location"],
            "provenance_status": provenance["status"],
            "value_m^-1_s^-1": record["value"],
            "uncertainty_m^-1_s^-1": record.get("uncertainty") or 0,
            "pH": record.get("pH") if record.get("pH") is not None else "unknown",
            "temperature": record.get("temperature") or "unknown",
            "environment": record.get("environment") or "unknown",
            "species": record.get("species") or "unknown",
            "protonation": conditions.get("protonation") or "unknown",
            "oxygen_conditions": conditions.get("oxygen_conditions") or "unknown",
            "method": record.get("method") or "unknown",
            "source": provenance.get("source") or record.get("source"),
            "allowed_use": ";".join(record.get("allowed_uses", [])),
            "condition_matched_allowed": "condition_matched_prediction" in record.get("allowed_uses", []),
            "limitations": record.get("limitations"),
            "aggregation_policy": "separate condition-specific record; do not fit a universal rate",
        })
    return rows


def _provenance_rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        source = list(csv.DictReader(handle))
    fields = [
        "symbol", "definition", "value_or_range", "unit", "species", "charge_state",
        "environment", "temperature", "pH", "method", "uncertainty", "source",
        "status", "limitations", "use_class", "code_location",
    ]
    return [
        {field: row.get(field, "") for field in fields}
        | {"permitted_use": row.get("use_class", "")}
        for row in source
    ]


def _unavailable_rows(authority: dict) -> list[dict]:
    conclusions = {
        "hyperfine": "Prohibits an exact-system hyperfine Hamiltonian or mixing rate.",
        "association": "Prohibits a quantitative encounter-formation flux.",
        "k_doublet": "Prohibits a quantitative D-channel reaction yield.",
        "k_quartet": "Prohibits a quantitative Q-channel reaction yield.",
        "ratio": "Prohibits assigning a preferred manifold or physical selectivity range.",
        "escape": "Prohibits a quantitative reaction-versus-escape branching ratio.",
        "duration": "Prohibits a physical encounter time axis and endpoint.",
        "preparation": "Prohibits choosing a chemically prepared initial D/Q state.",
        "t1": "Prohibits an encounter-specific relaxation rate.",
        "t2": "Prohibits quantitative dephasing or coherence-lifetime claims.",
        "exchange": "Prohibits a quantitative exchange Hamiltonian.",
        "dipolar": "Prohibits a quantitative dipolar Hamiltonian and orientation average.",
        "geometry": "Prohibits a physical association or encounter model.",
        "o2": "Prohibits transfer of gas, matrix, or protein O2 tensors into solution.",
        "o2_relaxation": "Prohibits assigning a solution-encounter O2 relaxation rate.",
        "loss": "Prohibits compartment-level H2O2 prediction.",
    }
    rows = []
    for item in authority["unavailable_exact_system_parameters"]:
        text = item.lower()
        # Classify specific phrases before short substrings: for example,
        # "duration" contains "ratio" and "H2O2 loss" contains "O2".
        if "hyperfine" in text:
            key = "hyperfine"
        elif "association" in text:
            key = "association"
        elif "physical k_quartet_s/k_doublet_s" in text:
            key = "ratio"
        elif "k_doublet" in text:
            key = "k_doublet"
        elif "k_quartet" in text:
            key = "k_quartet"
        elif "escape" in text:
            key = "escape"
        elif "duration" in text:
            key = "duration"
        elif "h2o2 loss" in text or "compartment-specific" in text:
            key = "loss"
        elif "semiquinone t1" in text:
            key = "t1"
        elif "semiquinone t2" in text or "semiquinone tm" in text:
            key = "t2"
        elif "exchange" in text:
            key = "exchange"
        elif "dipolar" in text:
            key = "dipolar"
        elif "o2 relaxation" in text:
            key = "o2_relaxation"
        elif "o2" in text or "zfs" in text:
            key = "o2"
        else:
            key = "geometry"
        why_needed = {
            "hyperfine": "Defines nuclear-spin-driven electron-spin mixing in the exact radical.",
            "association": "Connects bulk semiquinone and oxygen concentrations to encounter formation.",
            "k_doublet": "Sets doublet-resolved electron-transfer flux.",
            "k_quartet": "Sets quartet-resolved electron-transfer flux.",
            "ratio": "Determines reaction selectivity between the two manifolds.",
            "escape": "Sets competition between reaction and encounter separation.",
            "duration": "Sets a physical observation endpoint or residence time.",
            "preparation": "Defines the encounter birth density matrix.",
            "t1": "Defines longitudinal relaxation within the encounter.",
            "t2": "Defines coherence decay under a physically matched model.",
            "exchange": "Sets the isotropic electron-electron coupling magnitude.",
            "dipolar": "Sets anisotropic coupling and requires geometry/orientation.",
            "geometry": "Defines separation, orientation, association, and residence physics.",
            "o2": "Defines encounter-specific oxygen Zeeman/ZFS evolution.",
            "o2_relaxation": "Defines spin relaxation of O2 within the solution encounter.",
            "loss": "Defines compartment-specific removal of downstream H2O2.",
        }
        rows.append({
            "missing_parameter": item,
            "why_needed": why_needed[key],
            "status": "unavailable",
            "sensitivity_coordinate_allowed": "yes, with explicit opt-in",
            "range_status": "illustrative only; no defensible physical range",
            "conclusion_prohibited": conclusions[key],
            "allowed_treatment": "bounded illustrative sensitivity coordinate only",
            "source_of_status": "configuration authority and targeted evidence review",
        })
    explicit_additions = [
        (
            "encounter association rate",
            "association",
            "bounded illustrative sensitivity coordinate only",
        ),
        (
            "encounter geometry separated from dipolar coupling",
            "geometry",
            "illustrative geometry scenarios only; none used in the baseline",
        ),
    ]
    for name, key, treatment in explicit_additions:
        if not any(row["missing_parameter"] == name for row in rows):
            rows.append({
                "missing_parameter": name,
                "why_needed": why_needed[key],
                "status": "unavailable",
                "sensitivity_coordinate_allowed": "yes, with explicit opt-in",
                "range_status": "illustrative only; no defensible physical range",
                "conclusion_prohibited": conclusions[key],
                "allowed_treatment": treatment,
                "source_of_status": "references/doxorubicin_parameter_review.md",
            })
    rows.append({
        "missing_parameter": "chemically prepared initial D/Q population",
        "why_needed": "Defines the encounter birth density matrix.",
        "status": "unavailable",
        "sensitivity_coordinate_allowed": "yes, as benchmark mixtures only",
        "range_status": "illustrative only; no measured preparation distribution",
        "conclusion_prohibited": conclusions["preparation"],
        "allowed_treatment": "unpolarized and manifold-mixture benchmarks only",
        "source_of_status": "references/doxorubicin_parameter_review.md",
    })
    return rows


def _validation_dataset_rows(authority: dict) -> list[dict]:
    reserved = {
        "dataset_id", "enabled", "source", "limitations", "suitability", "method",
        "unit", "dose", "environment", "preparation", "temperature", "pH", "oxygen",
    }
    rows = []
    for dataset_id, record in authority["validation_datasets"].items():
        values = {key: value for key, value in record.items() if key not in reserved}
        exact_system = (
            record.get("preparation")
            or record.get("environment")
            or dataset_id.replace("_", " ")
        )
        identifier = dataset_id.lower()
        if "h2o2" in identifier:
            observable = "hydrogen peroxide (H2O2)"
        elif "mitosox" in identifier:
            observable = "relative MitoSOX fluorescence (oxidant proxy)"
        elif "superoxide" in identifier or "negative_control" in identifier:
            observable = "superoxide proxy rate"
        else:
            observable = "reported ROS-related observable"
        conditions = "; ".join(
            f"{label}={value}" for label, value in (
                ("dose", record.get("dose")),
                ("environment", record.get("environment") or record.get("preparation")),
                ("temperature", record.get("temperature")),
                ("pH", record.get("pH")),
                ("oxygen", record.get("oxygen")),
            ) if value not in (None, "")
        )
        rows.append({
            "dataset_id": dataset_id,
            "exact_chemical_or_biological_system": exact_system,
            "observable": observable,
            "conditions": conditions or "not fully reported",
            "reported_result": values,
            "reported_values_and_uncertainties": values,
            "unit": record.get("unit"), "dose": record.get("dose"),
            "environment": record.get("environment") or record.get("preparation") or "not fully reported",
            "temperature": record.get("temperature") or "unknown",
            "pH": record.get("pH") if record.get("pH") is not None else "unknown",
            "oxygen": record.get("oxygen") or "unknown", "method": record.get("method"),
            "source": record.get("source"),
            "suitability": record.get("suitability", "context/validation only; exact condition match required"),
            "usable_for_direct_validation": "no" if not record.get("enabled") else "condition-specific only",
            "enabled_for_fit": record.get("enabled"),
            "limitations": record.get("limitations", "Conditions and normalization differ; do not pool or fit."),
        })
    return rows


def _traceability_rows() -> list[dict]:
    critique_author = (
        "Critique-author recommendation; original reviewer correspondence unavailable"
    )
    rows = [
        ("A1", "Correct the physical spin system", critique_author,
         "corrected with evidence", "Complete S=1/2 x 1 six-state space and exact D/Q projectors.",
         "No chemically measured encounter preparation is available.",
         "spin_chemistry.py", "P_DOUBLET; P_QUARTET; initial_density", "Table 1", "SpinAlgebraTests; InitialStateAndHamiltonianTests"),
        ("A2", "Construct the anthracycline reaction network", critique_author,
         "corrected with evidence", "Upstream semiquinone, primary superoxide, HO2/O2-minus, H2O2, and loss are separate stages.",
         "Formal semiquinone protonation remains condition-dependent; the downstream calculation is a fixed-pH pulse.",
         "spin_chemistry.py; README.md", "primary_superoxide_formation; downstream_ros_species_resolved", "Figures 1 and 9; Table 1", "StagedChemistryTests"),
        ("A3", "Encode a real spin Hamiltonian", critique_author,
         "partially addressed", "Hermitian six-state Zeeman, exchange, dipolar, O2 ZFS, and local electronic-field terms are implemented and unit documented.",
         "Exact doxorubicin hyperfine, geometry, J, dipolar, and encounter O2 tensors are unavailable and disabled or sensitivity-only.",
         "spin_chemistry.py; configs/doxorubicin_parameters.json", "hamiltonian; EncounterParameters", "Tables 1-3", "InitialStateAndHamiltonianTests; AuthorityAndModeTests"),
        ("A4", "Replace CZ repetition with Hamiltonian propagation", critique_author,
         "corrected with evidence", "Classical propagation uses exp(L dt); the circuit inserts one dense Hamiltonian-derived unitary.",
         "No product-formula circuit exists; Trotter terminology and convergence claims are withdrawn.",
         "spin_chemistry.py", "propagate_encounter_reference; execute_three_qubit_unitary_circuit", "Figure 11; Table 7", "CoherentValidationTests; EncounterDynamicsTests"),
        ("A5", "Correct the singlet measurement", critique_author,
         "superseded by corrected model", "The invalid singlet/triplet observable is replaced by exact doublet/quartet projectors and numerical expectations.",
         "No finite-shot or measurement-gate experiment is executed.",
         "spin_chemistry.py", "P_DOUBLET; P_QUARTET", "Figures 2, 3, and 11; Table 7", "SpinAlgebraTests; CoherentValidationTests"),
        ("A6", "Implement spin-selective reaction dynamics", critique_author,
         "corrected with evidence", "Separate kD, kQ, and escape loss channels are integrated with survival closure.",
         "D/Q selectivity is a conditional modeling assumption; exact rates are unavailable.",
         "spin_chemistry.py", "propagate_encounter_reference; propagate_encounter_independent", "Figures 2-7; Tables 4-6", "EncounterDynamicsTests; IndependentEncounterValidationTests"),
        ("A7", "Add a matched classical comparison", critique_author,
         "partially addressed", "Independent matrix-exponential and adaptive solvers match for the same full open-system problem; coherent Qiskit matches the same dense unitary.",
         "The coherent circuit does not implement reaction, relaxation, or escape, so no open-system circuit validation is claimed.",
         "spin_chemistry.py; paper_analysis.py", "_solver_validation; _circuit_rows", "Figures 8 and 11; Tables 6 and 7", "IndependentEncounterValidationTests; CoherentValidationTests"),
        ("A8", "Add computation-time and memory benchmarks",
         "Explicitly attributed by the critique to Reviewer 1; original correspondence unavailable",
         "corrected with evidence", "A bounded table separately times matched coherent NumPy/Qiskit work and matched classical open-system solvers, with phase and memory scope.",
         "Process RSS is a whole-process high-water mark; no quantum hardware, transpilation, sampling, scaling, or performance advantage is evaluated.",
         "paper_analysis.py", "_runtime_memory_benchmark", "Table 10", "PaperPipelineTests"),
        ("A9", "Test scaling with hyperfine-coupled nuclei", critique_author,
         "withdrawn with justification", "Scaling and efficiency claims are removed; no artificial nuclear-spin model is introduced.",
         "Exact doxorubicin hyperfine tensors are unavailable. Withdrawal is not reviewer acceptance.",
         "README.md; metadata/run_manifest.json", "SCIENTIFIC_LIMITATIONS", "Table 9", "DocumentationTraceabilityTests; PaperPipelineTests"),
        ("A10", "Replace the T1/T2 treatment", critique_author,
         "superseded by corrected model", "The old amplitude/phase partition is retired; the retained generator is explicitly phenomenological isotropic Lindblad depolarization.",
         "It is not a measured T1, T2, or Tphi model; exact-system relaxation values are unavailable.",
         "spin_chemistry.py; README.md", "_lindblad; propagate_encounter_reference", "Figures 5, 7, and 8; Tables 2-3", "EncounterDynamicsTests; IndependentEncounterValidationTests"),
        ("A11", "Replace or rescale seconds-long parameters", critique_author,
         "corrected with evidence", "All encounter axes use t*kref and rates/kref; seconds are explicitly illustrative conversions from kref.",
         "The reference rate is not a measured encounter clock.",
         "README.md; sensitivity_analysis.py; paper_analysis.py", "encounter_result_row", "Figures 2, 4-8; Tables 4-6", "PaperPipelineTests"),
        ("A12", "Rebuild the magnetic-field model", critique_author,
         "partially addressed", "Static Zeeman inputs in tesla are distinct from a named local electronic-field proxy in rad/s.",
         "No stochastic magnetic-noise distribution, spectrum, correlation time, or physical field-effect prediction is implemented.",
         "spin_chemistry.py; README.md", "hamiltonian; EncounterParameters", "Tables 1-3", "InitialStateAndHamiltonianTests; AuthorityAndModeTests"),
        ("A13", "Correct the legacy beta definition", critique_author,
         "superseded by corrected model", "The legacy beta variable and its conclusions are absent from active calculations and outputs.",
         "No old beta boundary or trend is retained; this scope reduction is not reviewer acceptance.",
         "paper_analysis.py; sensitivity_analysis.py", "run_sweep", "Table 1", "PaperPipelineTests"),
        ("A14", "Resolve implementation defects in public code", critique_author,
         "superseded by corrected model", "CZ depth, parity-as-singlet, one-kick, semiconductor calibration, and protocol-weighting code paths are retired.",
         "The current model does not reproduce or reinterpret those obsolete outputs.",
         "spin_chemistry.py; tests/test_paper_pipeline.py", "active six-state API", "Table 1", "PaperPipelineTests.test_active_files_do_not_reference_removed_two_qubit_api"),
        ("A15", "Regenerate all results", critique_author,
         "corrected with evidence", "Figures 1-11, Tables 1-10, captions, CSV sources, checksums, and manifests are regenerated from the corrected model.",
         "Removed amplitude/phase, protocol, one-kick, beta, and Trotter analyses are not rerun because their claims were retired.",
         "paper_analysis.py", "main", "Figures 1-11; Tables 1-10", "PaperPipelineTests"),
        ("B1", "Add a chemical reaction schematic", critique_author,
         "partially addressed", "The existing author-created Figure 1, its caption, README equations, and Table 1 document species, charges, unpaired-electron counts, multiplicities, and singlet-oxygen exclusion.",
         "No new detailed chemical drawing was added under the restricted-imagery scope.",
         "README.md; paper_analysis.py", "CAPTIONS; _model_overview_rows", "Figure 1; Table 1", "PaperPipelineTests"),
        ("B2", "Add the reviewer-requested model-comparison schematic",
         "Called reviewer-requested by the critique; original correspondence unavailable",
         "partially addressed", "Figure 1 plus Table 1 and Table 7 separate the full classical open-system model from the coherent-only circuit.",
         "No new side-by-side schematic was added; the visual request was reduced to text/table evidence and is not accepted by assumption.",
         "README.md; paper_analysis.py", "_model_comparison_rows; _circuit_rows", "Figure 1; Tables 1 and 7", "PaperPipelineTests"),
        ("B3", "Add a complete circuit diagram", critique_author,
         "partially addressed", "Table 7 and the Figure 11 caption document basis indices, supplied statevector, dense unitary, duration, and numerical D/Q post-processing.",
         "No preparation, measurement, noise, reaction, delay, or Trotter gates are executed; no new diagram was added.",
         "spin_chemistry.py; paper_analysis.py", "execute_three_qubit_unitary_circuit; _circuit_rows", "Figure 11; Table 7", "CoherentValidationTests"),
        ("B4", "Add a parameter-provenance table", critique_author,
         "corrected with evidence", "Complete machine-readable parameter provenance is retained with curated readable panels.",
         "Many exact encounter parameters remain unavailable; no surrogate is activated as evidence.",
         "configs/parameter_provenance.csv; paper_analysis.py", "_provenance_rows", "Table 2", "AuthorityAndModeTests; PaperPipelineTests"),
        ("B5", "Add validation and performance figures", critique_author,
         "partially addressed", "Numerical solver validation and coherent circuit consistency remain Figures 8 and 11; bounded runtime/memory is Table 10.",
         "Trotter, nuclear-spin scaling, and quantum-performance figures are withdrawn with their claims.",
         "paper_analysis.py; plotting.py", "_solver_validation; _runtime_memory_benchmark", "Figures 8 and 11; Tables 6, 7, and 10", "PaperPipelineTests"),
        ("B6", "Rebuild or remove the existing figures", critique_author,
         "corrected with evidence", "Figures 1-11 are corrected-model outputs with one CSV and caption each; Figure 1 alone is an author-created conceptual map.",
         "All quantitative encounter coordinates remain illustrative and non-predictive.",
         "paper_analysis.py; plotting.py", "main; CAPTIONS", "Figures 1-11", "PaperPipelineTests.test_every_rendered_figure_is_nonblank_and_has_source_csv"),
        ("B7", "Rebuild or delete the tables", critique_author,
         "corrected with evidence", "Tables 1-10 are rebuilt as complete CSV sources plus readable PDF and numbered high-resolution PNG pages.",
         "Presentation panels round values and select columns, while the linked CSV retains full precision and provenance.",
         "paper_analysis.py; table_rendering.py", "TABLE_PRESENTATION; render_table_from_csv", "Tables 1-10", "PaperPipelineTests"),
    ]
    keys = (
        "requirement_id", "criticism", "request_attribution", "disposition",
        "repository_evidence", "remaining_limitation", "exact_file",
        "exact_code_object", "figure_or_table", "test",
    )
    return [
        {**dict(zip(keys, row)), "evidence_status": row[3]}
        for row in rows
    ]


def _circuit_rows(
    reference_rate: float, execute: bool, seed: int
) -> tuple[list[dict], dict[str, Any]]:
    capability = independent_circuit_capability()
    parameters = EncounterParameters(
        field_t=(0, 0, 1e-4), g_radical=2.0035, g_oxygen=2.0023,
        exchange_rad_s=2 * reference_rate, dipolar_rad_s=.1 * reference_rate,
        local_field_proxy_rad_s=(.3 * reference_rate, 0, 0),
    )
    duration = 1 / reference_rate
    physical = np.array([0, 1, 2, 4, 5, 6])
    unused = np.array([3, 7])
    energies, vectors = np.linalg.eigh(hamiltonian(parameters))
    unitary6 = (vectors * np.exp(-1j * energies * duration)) @ vectors.conj().T
    unitary8 = np.eye(8, dtype=complex)
    unitary8[np.ix_(physical, physical)] = unitary6
    rng = np.random.default_rng(seed)
    state6 = rng.normal(size=6) + 1j * rng.normal(size=6)
    state6 /= np.linalg.norm(state6)
    state8 = np.zeros(8, dtype=complex)
    state8[physical] = state6
    expected6 = unitary6 @ state6
    embedded8 = unitary8 @ state8
    pd8 = np.zeros((8, 8), dtype=complex)
    pq8 = np.zeros((8, 8), dtype=complex)
    pd8[np.ix_(physical, physical)] = P_DOUBLET
    pq8[np.ix_(physical, physical)] = P_QUARTET
    expected8 = np.zeros(8, dtype=complex)
    expected8[physical] = expected6

    numpy_metrics = [
        ("numpy_embedding_statevector_error", float(np.max(np.abs(embedded8[physical] - expected6))), TOLERANCES["circuit_statevector"]),
        ("numpy_embedding_doublet_observable_error", float(abs(np.vdot(embedded8, pd8 @ embedded8) - np.vdot(expected8, pd8 @ expected8))), TOLERANCES["circuit_observable"]),
        ("numpy_embedding_quartet_observable_error", float(abs(np.vdot(embedded8, pq8 @ embedded8) - np.vdot(expected8, pq8 @ expected8))), TOLERANCES["circuit_observable"]),
        ("numpy_embedding_leakage_probability", float(np.sum(np.abs(embedded8[unused]) ** 2)), TOLERANCES["circuit_leakage"]),
    ]
    qiskit_status = "not requested"
    reason = "--execute-circuit was not requested"
    implementation = "NumPy six-state unitary embedded in an 8-state array"
    qiskit_metrics: list[tuple[str, float, float]] = []
    qiskit_result: dict[str, Any] | None = None
    if execute and capability["available"]:
        qiskit_result = execute_three_qubit_unitary_circuit(
            parameters, duration, seed=seed
        )
        qiskit_status = "executed"
        reason = ""
        implementation = qiskit_result["gate_implementation"]
        qiskit_metrics = [
            ("qiskit_statevector_error", qiskit_result["max_statevector_error"], TOLERANCES["circuit_statevector"]),
            ("qiskit_doublet_observable_error", qiskit_result["doublet_observable_error"], TOLERANCES["circuit_observable"]),
            ("qiskit_quartet_observable_error", qiskit_result["quartet_observable_error"], TOLERANCES["circuit_observable"]),
            ("qiskit_leakage_probability", qiskit_result["physical_subspace_leakage_probability"], TOLERANCES["circuit_leakage"]),
            ("qiskit_basis_ordering_error", qiskit_result["max_basis_ordering_error"], TOLERANCES["circuit_basis_ordering"]),
        ]
    elif execute:
        qiskit_status = "unavailable"
        reason = str(capability["reason"])

    status = "qiskit_executed" if qiskit_status == "executed" else "numpy_embedding_only"
    scope = "closed coherent simulator and encoding consistency only"
    rows = [{
        "metric": metric,
        "validation_layer": "Qiskit" if metric.startswith("qiskit_") else "NumPy embedding",
        "value": value,
        "tolerance": tolerance,
        "passes": value <= tolerance,
        "status": status,
        "qiskit_available": bool(capability["available"]),
        "qiskit_execution_status": qiskit_status,
        "implementation_type": implementation,
        "scope": scope,
        "reason": reason,
        "seed": seed,
        "qubits": 3,
        "duration_over_reference": 1.0,
        "duration_s": duration,
        "input_state_origin": (
            "normalized pseudorandom six-element complex statevector supplied "
            f"directly with seed {seed}"
        ),
        "state_preparation": "none; Statevector input is supplied numerically",
        "measurement": "none; no measurement gates or finite-shot sampling",
        "post_processing": (
            "numerical statevector difference, D/Q projector expectations, basis "
            "ordering, and unused-state leakage"
        ),
        "basis_order": ["000", "001", "010", "100", "101", "110"],
        "physical_indices": physical.tolist(),
        "unused_indices": unused.tolist(),
        "does_not_validate": [
            "independent Hamiltonian construction", "spin-selective reaction",
            "escape", "Lindblad relaxation", "upstream or downstream chemistry",
            "hardware performance", "quantum advantage",
        ],
    } for metric, value, tolerance in numpy_metrics + qiskit_metrics]
    record = {
        "status": status,
        "qiskit_available": bool(capability["available"]),
        "qiskit_execution_status": qiskit_status,
        "reason": reason,
        "implementation_type": implementation,
        "scope": scope,
        "all_metrics_pass": all(row["passes"] for row in rows),
    }
    if qiskit_result is not None:
        record["qiskit_version_scope"] = qiskit_result["scope"]
    return rows, record


CAPTIONS = {
    "figure01_model_overview": "Figure 1. Corrected model overview. Formal electron-transfer bookkeeping is AQ + e- -> AQ radical anion, followed by AQ radical anion + ground-state triplet O2 -> AQ + O2 radical anion; the modeled downstream pool uses HO2/O2 radical-anion equilibrium and two radical equivalents per H2O2. Closed-shell AQ has zero unpaired electrons, semiquinone and superoxide each have one (doublet), and ground-state O2 has two (triplet); singlet oxygen is excluded. The semiquinone protonation/formal charge is condition-dependent and is not assigned as a measured active-model input. A semiquinone-triplet-O2 encounter enters the six-state doublet/quartet stage, followed by competing conditional reaction or escape. Conceptual workflow created by the authors; not a simulation output. The diagram is non-quantitative and does not assert coherent preparation or a measured encounter.",
    "figure02_benchmark_trajectories": "Figure 2. Benchmark encounter trajectories. Surviving doublet and quartet populations, total survival, cumulative doublet and quartet reaction yields, primary-superoxide yield, and escape yield are shown against normalized time for five benchmark mixtures. At every sampled time, YD + YQ + Yescape + Psurvival = 1 within numerical tolerance. The settings kD/kref=1, kQ/kD=0.1, kescape/kref=1, omega_local/kref=1, and gammaR/kref=gammaO/kref=0.1 are illustrative sensitivity coordinates, not measured parameters.",
    "figure03_initial_state_comparison": "Figure 3. Initial-state comparison. Final reaction, escape, and unresolved survival outcomes are compared for the unpolarized encounter, doublet-manifold benchmark, quartet-manifold benchmark, and pD=0.25 and pD=0.75 mixtures. The doublet and quartet cases are computational limiting benchmarks rather than demonstrated chemically prepared states. Primary superoxide is the sum of integrated D and Q reaction yields, not a final spin population.",
    "figure04_mixing_escape_heatmaps": "Figure 4. Mixing-versus-escape heatmaps. Primary-superoxide yield per encounter is calculated over normalized local D/Q mixing and escape rates for kQ/kD=0, 0.1, 1, and 10. The kQ/kD=1 panel is the spin-independent null condition. Zero and logarithmic nonzero coordinates are illustrative; no displayed range is experimentally established.",
    "figure05_mixing_relaxation_heatmaps": "Figure 5. Mixing-versus-relaxation heatmaps. Primary-superoxide yield is calculated for an unpolarized benchmark over normalized local mixing and equal radical/O2 relaxation rates. The fixed settings are kD/kref=1, kescape/kref=1, and duration=8/kref; panels show kQ/kD=0, 0.1, 1, and 10, with equality marked as the spin-independent null. All coordinates are illustrative and non-predictive.",
    "figure06_reaction_selectivity": "Figure 6. Reaction-selectivity sensitivity. Primary-superoxide yield is plotted against kQ/kD for unpolarized, doublet, quartet, pD=0.25, and pD=0.75 benchmarks across three mixing and escape settings. The vertical line marks kQ/kD=1 as the spin-independent null. The kQ/kD axis is a sensitivity coordinate, not a measured or chemically defensible probability distribution or physical range.",
    "figure07_controls": "Figure 7. Controls and limiting cases. Baseline sensitivity, zero mixing, kQ=kD spin-independent reaction, fast relaxation, rapid escape, doublet benchmark, and quartet benchmark cases show primary reaction, escape, and unresolved survival. Stacked outcomes make probability accounting visible. All rates and initial-state restrictions are computational controls, not established encounter parameters or preparations.",
    "figure08_solver_validation": "Figure 8. Independent numerical validation. A separately constructed adaptive Dormand-Prince 5(4) solver is compared with the constant-generator matrix-exponential reference for full density matrices, doublet and quartet populations, survival, reaction yields, escape yield, and probability balance. Dashed lines show declared absolute tolerances. The symmetric-log display includes exact zero without replacing it by an artificial logarithmic floor. Agreement validates numerical implementation only, not physical encounter inputs or the coherent circuit.",
    "figure09_downstream_kinetics": "Figure 9. Downstream ROS kinetics. Total superoxide-family radical pool, HO2, O2-minus, accumulated H2O2, and accumulated H2O2 loss are shown for spontaneous dismutation, SOD-mediated dismutation, and SOD plus H2O2 loss. The calculation starts from an illustrative 10 micromolar single radical pulse at fixed pH 7.4, assumes rapid acid-base equilibrium and constant rates/SOD, and consumes two radical equivalents per H2O2. No spin population is mapped directly to H2O2, and the trajectories are not cellular predictions.",
    "figure10_bulk_rate_comparison": "Figure 10. Literature bulk-rate comparison. Actual numerical values stored in the configuration JSON and matched provenance CSV are plotted with reported uncertainty, pH, temperature, environment, species/protonation, method, and source identifiers. Missing conditions remain visibly marked and records are not combined into a universal rate. These M^-1 s^-1 bulk total constants are not encounter-level kD or kQ in s^-1.",
    "figure11_circuit_validation": "Figure 11. Three-qubit coherent-embedding validation. Physical six-state basis indices [0,1,2,4,5,6] are embedded in an eight-state register; indices [3,7] are unused. A normalized pseudorandom six-component statevector (seed 1729) is supplied numerically, not prepared by gates. One dense 8x8 UnitaryGate contains U=exp(-iH t) for t=1/kref; no Trotter, delay, noise, reaction, measurement, or sampling gates are executed. Statevector, numerical D/Q projector expectations, basis ordering, and unused-state leakage are compared with the direct six-state result. This validates coherent simulator and encoding consistency only; it does not independently validate H, the open-system model, chemistry, hardware performance, or quantum advantage.",
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--mode", choices=("quick", "full"), help="analysis mode")
    mode.add_argument("--quick", action="store_true", help="small grids for CI and rapid review")
    mode.add_argument("--full", action="store_true", help="publication grid (default)")
    parser.add_argument("--reference-rate-s", type=float, default=1e6)
    parser.add_argument("--samples", type=int, help="trajectory samples; mode-dependent default")
    parser.add_argument("--grid-size", type=int, help="points per normalized axis (4-41)")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--provenance", type=Path, default=DEFAULT_PROVENANCE)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results" / "paper")
    parser.add_argument("--formats", default="png,pdf,svg")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--execute-circuit", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> None:
    started = time.perf_counter()
    args = build_parser().parse_args(argv)
    mode = args.mode or ("quick" if args.quick else "full")
    samples = args.samples if args.samples is not None else (61 if mode == "quick" else 241)
    grid_size = args.grid_size if args.grid_size is not None else (5 if mode == "quick" else 17)
    if not np.isfinite(args.reference_rate_s) or args.reference_rate_s <= 0:
        raise SystemExit("--reference-rate-s must be finite and positive")
    if isinstance(samples, bool) or not 11 <= samples <= 2001:
        raise SystemExit("--samples must be an integer from 11 through 2001")
    if isinstance(grid_size, bool) or not 4 <= grid_size <= 41:
        raise SystemExit("--grid-size must be an integer from 4 through 41")
    if not 72 <= args.dpi <= 1200:
        raise SystemExit("--dpi must be from 72 through 1200")
    if isinstance(args.seed, bool):
        raise SystemExit("--seed must be an integer")
    formats = tuple(dict.fromkeys(
        part.strip().lower() for part in args.formats.split(",") if part.strip()
    ))
    if not formats or set(formats) - {"png", "pdf", "svg"}:
        raise SystemExit("--formats must be a comma list drawn from png,pdf,svg")

    config, provenance = args.config.resolve(), args.provenance.resolve()
    authority = validate_authority_bundle(config, provenance)
    metadata = _source_metadata(config, provenance)
    output_dir = args.output_dir.resolve()
    paths = _prepare_output(output_dir, args.overwrite)

    data_sets: dict[str, list[dict]] = {}
    data_sets["figure01_model_overview"] = _model_overview_rows()
    (
        data_sets["figure02_benchmark_trajectories"],
        data_sets["figure03_initial_state_comparison"],
    ) = _benchmark_data(args.reference_rate_s, samples)
    data_sets["figure04_mixing_escape_heatmaps"] = run_mixing_escape_sweep(
        args.reference_rate_s, grid_size
    )
    data_sets["figure05_mixing_relaxation_heatmaps"] = run_mixing_relaxation_sweep(
        args.reference_rate_s, grid_size
    )
    data_sets["figure06_reaction_selectivity"] = run_selectivity_sweep(
        args.reference_rate_s, grid_size
    )
    data_sets["figure07_controls"] = _controls(args.reference_rate_s)
    (
        data_sets["figure08_solver_validation"], solver_summary,
    ) = _solver_validation(args.reference_rate_s, samples)
    data_sets["figure09_downstream_kinetics"] = _downstream_data(authority, samples)
    data_sets["figure10_bulk_rate_comparison"] = _bulk_rows(authority, provenance)
    (
        data_sets["figure11_circuit_validation"], circuit_record,
    ) = _circuit_rows(args.reference_rate_s, args.execute_circuit, args.seed)
    for rows in data_sets.values():
        _augment_rows(rows, metadata, config, provenance)

    source_paths: dict[str, Path] = {}
    for stem, rows in data_sets.items():
        source_paths[stem] = write_csv(rows, paths["data"] / f"{stem}.csv")

    # Plot only values read back from the saved source CSVs.
    plotters = {
        "figure01_model_overview": plot_model_overview,
        "figure02_benchmark_trajectories": plot_benchmark_trajectories,
        "figure03_initial_state_comparison": plot_initial_state_comparison,
        "figure04_mixing_escape_heatmaps": plot_mixing_escape_heatmaps,
        "figure05_mixing_relaxation_heatmaps": plot_mixing_relaxation_heatmaps,
        "figure06_reaction_selectivity": plot_reaction_selectivity,
        "figure07_controls": plot_controls,
        "figure08_solver_validation": lambda rows: plot_solver_validation(rows, TOLERANCES),
        "figure09_downstream_kinetics": plot_downstream,
        "figure10_bulk_rate_comparison": plot_bulk_rates,
    }
    plotters["figure11_circuit_validation"] = plot_circuit_validation
    generated_figures: list[Path] = []
    for stem, plotter in plotters.items():
        rows_from_disk = read_csv(source_paths[stem])
        generated_figures.extend(save_figure(
            plotter(rows_from_disk), paths["figures"] / stem, formats, args.dpi
        ))

    provenance_rows = _provenance_rows(provenance)
    unavailable_rows = _unavailable_rows(authority)
    baseline_rows = [dict(row) for row in data_sets["figure03_initial_state_comparison"]]
    extrema = []
    sensitivity_sources = (
        "figure04_mixing_escape_heatmaps",
        "figure05_mixing_relaxation_heatmaps",
        "figure06_reaction_selectivity",
    )
    for source in sensitivity_sources:
        rows = data_sets[source]
        for kind, selection in (
            ("minimum", min(rows, key=lambda row: row["primary_superoxide_yield_per_encounter"])),
            ("maximum", max(rows, key=lambda row: row["primary_superoxide_yield_per_encounter"])),
        ):
            extrema.append({
                "summary_type": f"selected_grid_{kind}",
                "source_grid": source,
                **selection,
            })
    control_summary = [
        {"summary_type": "control", "source_grid": "figure07_controls", **row}
        for row in data_sets["figure07_controls"]
    ]
    for row in baseline_rows:
        row["survival_probability"] = row["unresolved_probability"]
    table_sets = {
        "table01_original_versus_corrected_model": _model_comparison_rows(),
        "table02_parameter_provenance": provenance_rows,
        "table03_unavailable_parameters_and_consequences": unavailable_rows,
        "table04_baseline_results": baseline_rows,
        "table05_controls_and_extrema": control_summary + extrema,
        "table06_independent_solver_validation": solver_summary,
        "table07_circuit_validation": data_sets["figure11_circuit_validation"],
        "table08_experimental_validation_evidence": _validation_dataset_rows(authority),
        "table09_requirements_traceability": _traceability_rows(),
        "table10_runtime_memory_benchmark": _runtime_memory_benchmark(
            args.reference_rate_s, args.seed, args.execute_circuit, metadata
        ),
    }
    for rows in table_sets.values():
        _augment_rows(rows, metadata, config, provenance)
    table_linkage: dict[str, dict[str, Any]] = {}
    for stem, rows in table_sets.items():
        fields = _presentation_fields(TABLE_PRESENTATION[stem])
        csv_path = write_csv(rows, paths["tables"] / f"{stem}.csv")
        markdown_path = write_markdown(
            rows, paths["tables"] / f"{stem}.md", fields
        )
        rendered = render_table_from_csv(
            csv_path,
            paths["tables"] / f"{stem}.pdf",
            TABLE_PRESENTATION[stem],
            COLUMN_LABELS,
            requested_dpi=args.dpi,
        )
        table_linkage[stem] = {
            **rendered,
            "source_csv": f"tables/{csv_path.name}",
            "markdown": f"tables/{markdown_path.name}",
            "pdf": f"tables/{rendered['pdf']}",
            "png_previews": [
                f"tables/{name}" for name in rendered["png_previews"]
            ],
        }

    ofat_rows = run_sweep(args.reference_rate_s, samples)
    _augment_rows(ofat_rows, metadata, config, provenance)
    write_csv(ofat_rows, paths["data"] / "one_factor_at_a_time_sensitivity.csv")

    caption_lines = [
        "# Generated figure captions", "",
        "All encounter coordinates are normalized by the declared kref and are illustrative computational bounds, not measured ranges, confidence intervals, or priors.", "",
    ]
    for stem, caption in CAPTIONS.items():
        if stem == "figure11_circuit_validation" and circuit_record.get("qiskit_execution_status") != "executed":
            caption += (
                f" Qiskit execution status: {circuit_record.get('qiskit_execution_status')}. "
                f"Reason: {circuit_record.get('reason')}. The plotted NumPy embedding comparison remains available."
            )
        caption_lines.extend([f"## {stem}", "", caption, ""])
        (paths["captions"] / f"{stem}_caption.txt").write_text(
            caption + "\n", encoding="utf-8"
        )
    (paths["metadata"] / "figure_captions.md").write_text(
        "\n".join(caption_lines), encoding="utf-8"
    )

    baseline_yields = [
        row["primary_superoxide_yield_per_encounter"] for row in baseline_rows
    ]
    grid_yields = [
        row["primary_superoxide_yield_per_encounter"]
        for source in sensitivity_sources for row in data_sets[source]
    ]
    worst_solver = max(row["worst_absolute_error"] for row in solver_summary)
    downstream_final = {
        scenario: [
            row for row in data_sets["figure09_downstream_kinetics"]
            if row["scenario_id"] == scenario
        ][-1]
        for scenario in (
            "spontaneous_only", "sod_mediated", "sod_plus_h2o2_loss"
        )
    }
    balance_sources = (
        "figure03_initial_state_comparison", "figure04_mixing_escape_heatmaps",
        "figure05_mixing_relaxation_heatmaps", "figure06_reaction_selectivity",
        "figure07_controls",
    )
    result_summary = {
        "scope": "non-predictive sensitivity study",
        "baseline_primary_superoxide_yield_range_per_encounter": [
            min(baseline_yields), max(baseline_yields)
        ],
        "baseline_by_initial_state": {
            row["initial_state_definition"]: row["primary_superoxide_yield_per_encounter"]
            for row in baseline_rows
        },
        "selected_grid_primary_superoxide_yield_range_per_encounter": [
            min(grid_yields), max(grid_yields)
        ],
        "worst_independent_solver_absolute_error": worst_solver,
        "maximum_encounter_probability_balance_error": max(
            row.get("probability_balance_error", 0)
            for source in balance_sources for row in data_sets[source]
        ),
        "maximum_downstream_radical_equivalent_balance_error_m": max(
            row["radical_equivalent_balance_error_m"]
            for row in data_sets["figure09_downstream_kinetics"]
        ),
        "downstream_final_h2o2_m": {
            key: value["hydrogen_peroxide_m"]
            for key, value in downstream_final.items()
        },
        "circuit_status": circuit_record.get("status"),
        "circuit_scope": circuit_record.get("scope"),
        "limitations": SCIENTIFIC_LIMITATIONS,
    }
    (paths["metadata"] / "paper_results_summary.json").write_text(
        json.dumps(result_summary, indent=2, default=_jsonable) + "\n",
        encoding="utf-8",
    )
    summary_lines = [
        "# Paper-results summary", "",
        "This is a non-predictive sensitivity study. Its normalized encounter settings do not describe a measured doxorubicin encounter.", "",
        f"Across the five baseline initial-state benchmarks, primary-superoxide yield per encounter ranged from {min(baseline_yields):.6g} to {max(baseline_yields):.6g}.",
        f"Across the selected grids, the range was {min(grid_yields):.6g} to {max(grid_yields):.6g}; these extrema are properties of this grid only.",
        f"The worst matrix-exponential versus Dormand-Prince absolute discrepancy was {worst_solver:.3e}.",
        f"The maximum encounter probability-balance error was {result_summary['maximum_encounter_probability_balance_error']:.3e}; the maximum downstream radical-equivalent balance error was {result_summary['maximum_downstream_radical_equivalent_balance_error_m']:.3e} M.",
        f"Circuit status: {circuit_record.get('status')}. Scope: {circuit_record.get('scope')}.", "",
        "The computations show how assumed selectivity, mixing, relaxation, and escape interact in the declared equations. They do not establish measured D/Q selectivity, encounter coherence, magnetic control, or biological ROS concentrations.", "",
    ]
    (paths["metadata"] / "paper_results_summary.md").write_text(
        "\n".join(summary_lines), encoding="utf-8"
    )

    command_args = argv if argv is not None else sys.argv[1:]
    invocation = " ".join([
        shlex.quote(sys.executable), shlex.quote(str(Path(__file__).resolve())),
        *[shlex.quote(argument) for argument in command_args],
    ])
    timestamp = datetime.now(timezone.utc)
    figure_linkage = {
        stem: {
            "source_data": str(source_paths[stem].relative_to(output_dir)),
            "caption": f"captions/{stem}_caption.txt",
            "rendered_formats": [
                f"figures/{stem}.{extension}" for extension in formats
            ],
        }
        for stem in data_sets
    }
    manifest = {
        "schema_version": "2.1",
        "analysis_mode": mode,
        "output_label": "NON-PREDICTIVE ROS_SPIN PAPER SENSITIVITY STUDY",
        "command_line_invocation": invocation,
        "timestamp_iso8601": timestamp.isoformat(),
        "timezone": "UTC",
        "random_seeds": {"circuit_state": args.seed, "plotting": "deterministic/no randomness"},
        "reference_rate_s^-1": args.reference_rate_s,
        "formats": list(formats),
        "dpi": args.dpi,
        "grid_definition": {
            "grid_size_request": grid_size,
            "nonzero_mixing_relaxation_escape_bounds": [1e-2, 1e2],
            "kq_over_kd_nonzero_bounds": [1e-2, 1e1],
            "includes_exact_zero": True,
            "includes_exact_kq_over_kd_one": True,
            "bounds_status": "illustrative computational bounds only",
        },
        "sample_counts": {
            "trajectory": samples, "final_sweep": 2, "ofat": samples,
            "runtime_benchmark_repetitions": 7,
            "runtime_benchmark_warmups": 2,
        },
        "numerical_tolerances": TOLERANCES,
        "environment": metadata,
        "source_state_acceptance": {
            "generation_commit_is_clean_source_identity": not metadata["dirty_tree"],
            "results_commit_may_differ_from_generation_commit": True,
            "acceptance_basis": (
                "recorded source commit plus source/configuration/provenance hashes; "
                "results are deposited in a later results-only commit"
            ),
        },
        "qiskit_execution": circuit_record,
        "skipped_features": [] if circuit_record.get("qiskit_execution_status") == "executed" else [{
            "feature": "Qiskit statevector execution", "reason": circuit_record.get("reason")
        }],
        "assumptions": [
            "encounter rates and frequencies are normalized illustrative coordinates",
            "kD/kref=1 defines scale but is not a measured kD",
            "initial D/Q states are computational benchmarks",
            "downstream kinetics use fixed pH, rapid acid-base equilibrium, constant SOD and rate constants, and one radical pulse",
            "two radical equivalents form one H2O2",
            "bulk M^-1 s^-1 constants are not converted to encounter s^-1 rates",
            "finite-time encounter yields retain unresolved survival and are not asymptotic",
            "local Lindblad coefficients define phenomenological isotropic depolarization, not measured T1/T2",
            "the local electronic-field term is static and is not stochastic magnetic noise or nuclear hyperfine dynamics",
        ],
        "scientific_limitations": SCIENTIFIC_LIMITATIONS,
        "figure_source_and_render_linkage": figure_linkage,
        "table_source_and_render_linkage": table_linkage,
        "table_files": sorted({
            item
            for linkage in table_linkage.values()
            for item in (
                linkage["source_csv"], linkage["markdown"], linkage["pdf"],
                *linkage["png_previews"],
            )
        }),
        "runtime_seconds": time.perf_counter() - started,
        "generated_file_hashes_sha256": {},
        "checksum_scope_note": "Self-referential metadata files are inventoried but excluded from their own embedded hash maps; output_checksums.csv hashes the final run manifest and all preceding artifacts.",
    }

    readme_path = output_dir / "README.md"
    known_files = sorted(
        path.relative_to(output_dir) for path in output_dir.rglob("*") if path.is_file()
    )
    known_files.extend([
        Path("README.md"), Path("metadata/run_manifest.json"),
        Path("metadata/output_checksums.csv"), Path("metadata/SHA256SUMS"),
    ])
    known_files = sorted(set(known_files))
    readme_lines = [
        "# ROS_Spin paper-results bundle", "",
        f"Generated by `{invocation}`.", "",
        "This directory is a complete non-predictive sensitivity-study run. Use `--overwrite` to replace it deliberately; otherwise choose a new output directory.", "",
        "## Artifact inventory", "",
        "Every figure reads its values from the same-named CSV in `data/`. Figure 11 always includes the NumPy six-to-eight-state embedding comparison and also includes Qiskit statevector metrics when execution was requested and available. Table 7 records the exact status.", "",
        "Tables 1-10 each have a complete machine-readable CSV, a multipage typeset PDF, and numbered high-resolution PNG page previews. The rendered panels use documented presentation columns and rounding; the linked CSV retains every row, every column, full numerical precision, and provenance fields. `Not available` is never rendered as numeric zero.", "",
    ]
    for relative in known_files:
        if str(relative).startswith("data/"):
            meaning = "machine-readable figure or sensitivity source data"
        elif str(relative).startswith("figures/"):
            meaning = "publication figure rendered from its source CSV"
        elif str(relative).startswith("tables/"):
            meaning = "publication-ready or machine-readable evidence/results table"
        else:
            meaning = "run metadata, captions, or computed summary"
        readme_lines.append(
            f"- `{relative}` — {meaning}; apply its recorded limitations."
        )
    readme_lines.extend([
        "", "## Exact reproduction", "",
        "From the repository root, run:", "",
        "```bash",
        ".venv/bin/python paper_analysis.py --mode full --output-dir results/paper --execute-circuit --overwrite",
        "```", "",
        "The quick verification run is:", "",
        "```bash",
        ".venv/bin/python paper_analysis.py --mode quick --output-dir results/paper-quick --execute-circuit --overwrite",
        "```", "",
        "", "## Scientific meaning and limitations", "", SCIENTIFIC_LIMITATIONS, "",
        "The downstream figure assumes fixed pH, rapid HO2/O2-minus equilibrium, constant SOD, constant rates, and one radical pulse. Bulk constants remain separate condition-specific records and are not encounter rates.", "",
        "All encounter yields are finite-time quantities with unresolved survival reported explicitly; none is an asserted asymptotic yield. Any seconds shown come from the illustrative reference-rate conversion and are not measured encounter times.", "",
        "The local-spin GKSL term is phenomenological isotropic depolarization, not measured T1/T2. The local electronic-field proxy is static, not stochastic magnetic noise or explicit nuclear hyperfine dynamics.", "",
        "Table 9 records every critique item A1-A15 and B1-B7. Original reviewer correspondence was unavailable; withdrawn or scope-reduced requests are not represented as reviewer acceptance.", "",
        "The manifest records the clean generation commit separately from the later results commit; source, configuration, and provenance hashes are the acceptance basis.", "",
    ])
    readme_path.write_text("\n".join(readme_lines), encoding="utf-8")

    files_before_manifest = sorted(
        path for path in output_dir.rglob("*") if path.is_file()
    )
    manifest["generated_files"] = sorted(
        {str(path.relative_to(output_dir)) for path in files_before_manifest}
        | {
            "metadata/run_manifest.json",
            "metadata/output_checksums.csv",
            "metadata/SHA256SUMS",
        }
    )
    manifest["generated_file_hashes_sha256"] = {
        str(path.relative_to(output_dir)): _sha256(path)
        for path in files_before_manifest
    }
    manifest_path = paths["metadata"] / "run_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, default=_jsonable) + "\n", encoding="utf-8"
    )
    checksum_rows = []
    for path in sorted(path for path in output_dir.rglob("*") if path.is_file()):
        relative = str(path.relative_to(output_dir))
        checksum_rows.append({
            "path": relative,
            "sha256": _sha256(path),
            "bytes": path.stat().st_size,
            "artifact_kind": relative.split("/", 1)[0] if "/" in relative else "run_root",
        })
    write_csv(checksum_rows, paths["metadata"] / "output_checksums.csv")
    hash_lines = [
        f"{_sha256(path)}  {path.relative_to(output_dir)}"
        for path in sorted(path for path in output_dir.rglob("*") if path.is_file())
    ]
    (paths["metadata"] / "SHA256SUMS").write_text(
        "\n".join(hash_lines) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": "complete", "mode": mode, "output_dir": str(output_dir),
        "runtime_seconds": manifest["runtime_seconds"],
        "data_tables": len(data_sets) + 1,
        "figure_files": len(generated_figures),
        "publication_tables": len(table_sets),
        "circuit_status": circuit_record.get("status"),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
