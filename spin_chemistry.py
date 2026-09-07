"""Evidence-gated doxorubicin semiquinone--triplet-oxygen model.

The reference encounter implementation propagates the complete classical
six-state density matrix for spin-1/2 x spin-1. Semiquinone formation and
encounter association are deliberately separate classical stages. A numeric
zero used for a disabled interaction is only an algebraic placeholder and is
never interpreted as evidence that the physical interaction vanishes.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import numpy as np


EVIDENCE_BACKED = "evidence_backed"
SENSITIVITY = "sensitivity"
ACTIVATABLE_STATUSES = {"measured", "calculated", "fitted", "approved_surrogate"}
PARAMETER_STATUSES = ACTIVATABLE_STATUSES | {
    "assumed", "illustrative", "surrogate", "unavailable"
}
PARAMETER_FIELDS = {
    "value", "unit", "status", "enabled", "species", "environment",
    "temperature", "pH", "source", "limitations",
}


class ModelPolicyError(ValueError):
    """Raised when an execution request violates the evidence policy."""


def spin_matrices(spin: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dimension = int(round(2 * spin + 1))
    magnetic_numbers = np.arange(spin, -spin - 1, -1)
    raising = np.zeros((dimension, dimension), complex)
    for column in range(1, dimension):
        m = magnetic_numbers[column]
        raising[column - 1, column] = np.sqrt(spin * (spin + 1) - m * (m + 1))
    lowering = raising.T
    return (
        (raising + lowering) / 2,
        (raising - lowering) / (2j),
        np.diag(magnetic_numbers),
    )


I2, I3 = np.eye(2), np.eye(3)
R = tuple(np.kron(operator, I3) for operator in spin_matrices(0.5))
O = tuple(np.kron(I2, operator) for operator in spin_matrices(1.0))
I6 = np.eye(6, dtype=complex)
R_DOT_O = sum(radical @ oxygen for radical, oxygen in zip(R, O))
P_DOUBLET = (0.5 * I6 - R_DOT_O) / 1.5
P_QUARTET = (R_DOT_O + I6) / 1.5


def clebsch_gordan_states() -> dict[str, np.ndarray]:
    """Return coupled |J,M> states in the |ms,mS> product-basis ordering."""
    basis = np.eye(6, dtype=complex)
    return {
        "Q,+3/2": basis[0],
        "Q,+1/2": np.sqrt(2 / 3) * basis[1] + np.sqrt(1 / 3) * basis[3],
        "Q,-1/2": np.sqrt(1 / 3) * basis[2] + np.sqrt(2 / 3) * basis[4],
        "Q,-3/2": basis[5],
        "D,+1/2": np.sqrt(1 / 3) * basis[1] - np.sqrt(2 / 3) * basis[3],
        "D,-1/2": np.sqrt(2 / 3) * basis[2] - np.sqrt(1 / 3) * basis[4],
    }


@dataclass(frozen=True)
class EncounterParameters:
    """Encounter Hamiltonian (rad s^-1) and loss rates (s^-1).

    Zero defaults make no physical assertion; provenance and enabled/disabled
    state live in the authority bundle, not in this numerical value object.
    """

    field_t: tuple[float, float, float] = (0.0, 0.0, 0.0)
    g_radical: float = 0.0
    g_oxygen: float = 0.0
    exchange_rad_s: float = 0.0
    dipolar_rad_s: float = 0.0
    dipolar_axis: tuple[float, float, float] = (0.0, 0.0, 1.0)
    oxygen_zfs_d_rad_s: float = 0.0
    oxygen_zfs_e_rad_s: float = 0.0
    effective_hyperfine_rad_s: tuple[float, float, float] = (0.0, 0.0, 0.0)
    radical_relaxation_s: float = 0.0
    oxygen_relaxation_s: float = 0.0
    k_doublet_s: float = 0.0
    k_quartet_s: float = 0.0
    k_escape_s: float = 0.0


@dataclass(frozen=True)
class ExecutionContext:
    mode: str
    parameters: EncounterParameters
    duration_s: float | None
    authority: dict[str, Any]
    sensitivity_scenario: str | None = None
    non_predictive: bool = False


def _unit(vector: tuple[float, float, float]) -> np.ndarray:
    value = np.asarray(vector, float)
    norm = np.linalg.norm(value)
    if norm == 0:
        raise ValueError("axis must be nonzero")
    return value / norm


def hamiltonian(parameters: EncounterParameters) -> np.ndarray:
    """Return H/hbar for the six-state electronic encounter space."""
    beta_e_over_hbar = 8.79410005e10  # rad s^-1 T^-1
    field = np.asarray(parameters.field_t, float)
    result = beta_e_over_hbar * sum(
        field[index]
        * (parameters.g_radical * R[index] + parameters.g_oxygen * O[index])
        for index in range(3)
    )
    result += parameters.exchange_rad_s * R_DOT_O
    if parameters.dipolar_rad_s:
        axis = _unit(parameters.dipolar_axis)
        radical_axis = sum(axis[index] * R[index] for index in range(3))
        oxygen_axis = sum(axis[index] * O[index] for index in range(3))
        result += parameters.dipolar_rad_s * (R_DOT_O - 3 * radical_axis @ oxygen_axis)
    ox, oy, oz = O
    result += parameters.oxygen_zfs_d_rad_s * (oz @ oz - (2 / 3) * I6)
    result += parameters.oxygen_zfs_e_rad_s * (ox @ ox - oy @ oy)
    result += sum(
        parameters.effective_hyperfine_rad_s[index] * R[index]
        for index in range(3)
    )
    return np.asarray(result, complex)


def initial_density(scenario: str) -> np.ndarray:
    """Return an explicitly benchmark-only initial electronic density matrix."""
    if scenario == "unpolarized":
        return I6 / 6
    if scenario == "doublet":
        return P_DOUBLET / 2
    if scenario == "quartet":
        return P_QUARTET / 4
    raise ValueError("scenario must be unpolarized, doublet, or quartet")


def _lindblad(density: np.ndarray, operator: np.ndarray) -> np.ndarray:
    adjoint = operator.conj().T
    return operator @ density @ adjoint - 0.5 * (
        adjoint @ operator @ density + density @ adjoint @ operator
    )


def _rk4(rhs, initial, times, max_frequency=0.0):
    """Dependency-free RK4 with frequency-aware internal substeps."""
    values = [np.asarray(initial, complex)]
    value = values[0].copy()
    for start, stop in zip(times[:-1], times[1:]):
        interval = stop - start
        substeps = max(1, int(np.ceil(interval * max_frequency / 0.05)))
        dt = interval / substeps
        time = start
        for _ in range(substeps):
            k1 = rhs(time, value)
            k2 = rhs(time + dt / 2, value + dt * k1 / 2)
            k3 = rhs(time + dt / 2, value + dt * k2 / 2)
            k4 = rhs(time + dt, value + dt * k3)
            value = value + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            time += dt
        values.append(value.copy())
    return np.asarray(values)


def _matrix_exponential(matrix: np.ndarray) -> np.ndarray:
    """Scaling-and-squaring Padé(13) exponential without a SciPy dependency."""
    coefficients = (
        64764752532480000.0, 32382376266240000.0, 7771770303897600.0,
        1187353796428800.0, 129060195264000.0, 10559470521600.0,
        670442572800.0, 33522128640.0, 1323241920.0, 40840800.0,
        960960.0, 16380.0, 182.0, 1.0,
    )
    norm = np.linalg.norm(matrix, 1)
    squarings = max(0, int(np.ceil(np.log2(norm / 5.371920351148152)))) if norm else 0
    scaled = matrix / (2**squarings)
    identity = np.eye(matrix.shape[0], dtype=complex)
    a2 = scaled @ scaled
    a4 = a2 @ a2
    a6 = a4 @ a2
    odd = scaled @ (
        a6 @ (coefficients[13] * a6 + coefficients[11] * a4
              + coefficients[9] * a2)
        + coefficients[7] * a6 + coefficients[5] * a4
        + coefficients[3] * a2 + coefficients[1] * identity
    )
    even = (
        a6 @ (coefficients[12] * a6 + coefficients[10] * a4
              + coefficients[8] * a2)
        + coefficients[6] * a6 + coefficients[4] * a4
        + coefficients[2] * a2 + coefficients[0] * identity
    )
    result = np.linalg.solve(even - odd, even + odd)
    for _ in range(squarings):
        result = result @ result
    return result


def _exact_linear_trajectory(rhs, initial: np.ndarray, times: np.ndarray) -> np.ndarray:
    """Propagate a constant linear ODE by one exact matrix exponential per step."""
    dimension = initial.size
    basis = np.eye(dimension, dtype=complex)
    generator = np.column_stack([rhs(0.0, basis[:, index])
                                 for index in range(dimension)])
    step = _matrix_exponential(generator * (times[1] - times[0]))
    trajectory = [np.asarray(initial, complex)]
    for _ in times[1:]:
        trajectory.append(step @ trajectory[-1])
    return np.asarray(trajectory)


def propagate_encounter_reference(
    parameters: EncounterParameters,
    duration_s: float,
    initial_state: str = "unpolarized",
    samples: int = 1001,
) -> dict[str, Any]:
    """Propagate the full classical six-state density matrix by matrix exponential.

    d rho/dt = -i[H,rho] - {K + k_escape I,rho}/2 + L_relax(rho),
    K = k_D P_D + k_Q P_Q. Reaction and escape yields are time integrals
    of their respective fluxes. This is not an evidence-backed encounter
    prediction until every required encounter input is supported. The constant
    39-component linear system (36 density elements plus three accumulated
    yields) is advanced with a scaling-and-squaring Padé matrix exponential.
    """
    if duration_s <= 0 or samples < 2:
        raise ValueError("duration_s must be positive and samples >= 2")
    rates = (
        parameters.radical_relaxation_s, parameters.oxygen_relaxation_s,
        parameters.k_doublet_s, parameters.k_quartet_s, parameters.k_escape_s,
    )
    if any(rate < 0 for rate in rates):
        raise ValueError("rates must be nonnegative")

    coherent_hamiltonian = hamiltonian(parameters)
    loss = (
        parameters.k_doublet_s * P_DOUBLET
        + parameters.k_quartet_s * P_QUARTET
        + parameters.k_escape_s * I6
    )
    relaxation_operators = []
    if parameters.radical_relaxation_s:
        relaxation_operators.extend(
            np.sqrt(parameters.radical_relaxation_s) * operator for operator in R
        )
    if parameters.oxygen_relaxation_s:
        relaxation_operators.extend(
            np.sqrt(parameters.oxygen_relaxation_s) * operator for operator in O
        )

    def rhs(_time, state):
        density = state[:36].reshape(6, 6)
        derivative = -1j * (
            coherent_hamiltonian @ density - density @ coherent_hamiltonian
        ) - 0.5 * (loss @ density + density @ loss)
        for operator in relaxation_operators:
            derivative += _lindblad(density, operator)
        p_doublet_now = np.trace(P_DOUBLET @ density)
        p_quartet_now = np.trace(P_QUARTET @ density)
        survival_now = np.trace(density)
        return np.concatenate(
            (
                derivative.reshape(-1),
                np.array(
                    (
                        parameters.k_doublet_s * p_doublet_now,
                        parameters.k_quartet_s * p_quartet_now,
                        parameters.k_escape_s * survival_now,
                    ),
                    complex,
                ),
            )
        )

    times = np.linspace(0, duration_s, samples)
    initial = np.concatenate(
        (initial_density(initial_state).reshape(-1), np.zeros(3, complex))
    )
    trajectory = _exact_linear_trajectory(rhs, initial, times)
    densities = trajectory[:, :36].reshape(-1, 6, 6)
    p_doublet = np.real(np.einsum("ij,tji->t", P_DOUBLET, densities))
    p_quartet = np.real(np.einsum("ij,tji->t", P_QUARTET, densities))
    survival = np.real(np.trace(densities, axis1=1, axis2=2))
    doublet_yield, quartet_yield, escape_yield = (
        float(value) for value in np.real(trajectory[-1, 36:39])
    )
    return {
        "time_s": times,
        "density_matrices": densities,
        "p_doublet": p_doublet,
        "p_quartet": p_quartet,
        "survival": survival,
        "doublet_reaction_yield": doublet_yield,
        "quartet_reaction_yield": quartet_yield,
        "primary_superoxide_yield": doublet_yield + quartet_yield,
        "superoxide_yield": doublet_yield + quartet_yield,
        "escape_yield": escape_yield,
        "unresolved_probability": float(survival[-1]),
    }


# Backward-compatible name; the complete six-state implementation is reference.
propagate_encounter = propagate_encounter_reference


def form_semiquinone(quinone_m: float, formation_fraction: float) -> float:
    """Classical upstream preparation stage, kept outside the density matrix."""
    if quinone_m < 0 or not 0 <= formation_fraction <= 1:
        raise ValueError("quinone_m must be nonnegative and formation_fraction in [0,1]")
    return quinone_m * formation_fraction


def associate_encounters(
    semiquinone_m: float, oxygen_m: float, association_fraction: float
) -> float:
    """Classical encounter-association stage, limited by both reactants."""
    if min(semiquinone_m, oxygen_m) < 0 or not 0 <= association_fraction <= 1:
        raise ValueError("concentrations must be nonnegative and association_fraction in [0,1]")
    return min(semiquinone_m, oxygen_m) * association_fraction


def primary_superoxide_formation(encounter_m: float, reaction_yield: float) -> float:
    """Map reacted encounters to primary superoxide at 1:1 stoichiometry."""
    if encounter_m < 0 or not 0 <= reaction_yield <= 1 + 1e-9:
        raise ValueError("encounter_m must be nonnegative and reaction_yield in [0,1]")
    return encounter_m * min(reaction_yield, 1.0)


def bulk_superoxide_formation_rate(
    semiquinone_m: float, oxygen_m: float, k_bulk_m_inv_s: float
) -> float:
    """Measured bulk law r=k_bulk[SQ][O2], never an encounter k_D or k_Q."""
    if min(semiquinone_m, oxygen_m, k_bulk_m_inv_s) < 0:
        raise ValueError("bulk concentrations and rate constant must be nonnegative")
    return k_bulk_m_inv_s * semiquinone_m * oxygen_m


def acid_base_fractions(pH: float, pKa_ho2: float) -> tuple[float, float]:
    fraction_ho2 = 1 / (1 + 10 ** (pH - pKa_ho2))
    return fraction_ho2, 1 - fraction_ho2


def spontaneous_dismutation_rate(
    radical_pool_m: float,
    pH: float,
    pKa_ho2: float,
    k_ho2_ho2_m_inv_s: float,
    k_ho2_o2minus_m_inv_s: float,
    k_o2minus_o2minus_m_inv_s: float | None = None,
) -> float:
    """Species-resolved spontaneous H2O2 event rate in M s^-1."""
    fraction_ho2, fraction_o2minus = acid_base_fractions(pH, pKa_ho2)
    ho2 = fraction_ho2 * radical_pool_m
    o2minus = fraction_o2minus * radical_pool_m
    rate = k_ho2_ho2_m_inv_s * ho2**2 + k_ho2_o2minus_m_inv_s * ho2 * o2minus
    if k_o2minus_o2minus_m_inv_s is not None:
        rate += k_o2minus_o2minus_m_inv_s * o2minus**2
    return rate


def sod_dismutation_loss_rate(
    radical_pool_m: float,
    sod_m: float,
    k_sod_m_inv_s: float,
    pH: float,
    pKa_ho2: float,
) -> float:
    """SOD-mediated O2-minus disappearance rate in M s^-1."""
    _, fraction_o2minus = acid_base_fractions(pH, pKa_ho2)
    return k_sod_m_inv_s * sod_m * fraction_o2minus * radical_pool_m


def h2o2_loss_rate(h2o2_m: float, k_h2o2_loss_s: float | None) -> float:
    """Compartmental first-order H2O2 loss; None means omitted, not zero evidence."""
    return 0.0 if k_h2o2_loss_s is None else k_h2o2_loss_s * h2o2_m


def downstream_ros_species_resolved(
    superoxide0_m: float,
    duration_s: float,
    pH: float,
    pKa_ho2: float,
    k_ho2_ho2_m_inv_s: float,
    k_ho2_o2minus_m_inv_s: float,
    k_sod_m_inv_s: float | None = None,
    sod_m: float = 0.0,
    k_o2minus_o2minus_m_inv_s: float | None = None,
    k_h2o2_loss_s: float | None = None,
    samples: int = 201,
) -> dict[str, Any]:
    """Propagate spontaneous, SOD, and H2O2-loss stages without conflation."""
    numeric = [
        superoxide0_m, duration_s, pKa_ho2, k_ho2_ho2_m_inv_s,
        k_ho2_o2minus_m_inv_s, sod_m,
    ]
    numeric.extend(
        value for value in (k_sod_m_inv_s, k_o2minus_o2minus_m_inv_s,
                            k_h2o2_loss_s) if value is not None
    )
    if any(value < 0 for value in numeric) or samples < 2:
        raise ValueError("concentrations, durations, rates must be nonnegative; samples >= 2")

    def rhs(_time, values):
        radical_pool, h2o2 = np.real(values)
        spontaneous = spontaneous_dismutation_rate(
            radical_pool, pH, pKa_ho2, k_ho2_ho2_m_inv_s,
            k_ho2_o2minus_m_inv_s, k_o2minus_o2minus_m_inv_s,
        )
        sod_loss = 0.0
        if k_sod_m_inv_s is not None:
            sod_loss = sod_dismutation_loss_rate(
                radical_pool, sod_m, k_sod_m_inv_s, pH, pKa_ho2
            )
        return np.array(
            (-2 * spontaneous - sod_loss,
             spontaneous + 0.5 * sod_loss - h2o2_loss_rate(h2o2, k_h2o2_loss_s))
        )

    times = np.linspace(0, duration_s, samples)
    max_rate = 2 * (k_ho2_ho2_m_inv_s + k_ho2_o2minus_m_inv_s) * superoxide0_m
    max_rate += (k_sod_m_inv_s or 0.0) * sod_m + (k_h2o2_loss_s or 0.0)
    values = np.real(_rk4(rhs, np.array((superoxide0_m, 0.0)), times, max_rate))
    return {
        "time_s": times,
        "superoxide_pool_m": values[:, 0],
        "superoxide_m": values[:, 0],
        "hydrogen_peroxide_m": values[:, 1],
        "omitted_interactions": [
            name
            for name, value in (
                ("direct_O2minus_O2minus", k_o2minus_o2minus_m_inv_s),
                ("SOD_mediated_dismutation", k_sod_m_inv_s),
                ("H2O2_loss", k_h2o2_loss_s),
            )
            if value is None
        ],
    }


def downstream_ros(
    superoxide0_m, duration_s, k_spont_m_inv_s, k_sod_m_inv_s, sod_m,
    k_h2o2_loss_s=0.0, samples=201,
):
    """Legacy effective-rate sensitivity model; not species-resolved evidence."""
    def rhs(_time, values):
        superoxide, h2o2 = np.real(values)
        spontaneous = k_spont_m_inv_s * superoxide**2
        enzyme_loss = k_sod_m_inv_s * sod_m * superoxide
        return np.array(
            (-2 * spontaneous - enzyme_loss,
             spontaneous + 0.5 * enzyme_loss - k_h2o2_loss_s * h2o2)
        )

    times = np.linspace(0, duration_s, samples)
    scale = (
        2 * k_spont_m_inv_s * superoxide0_m
        + k_sod_m_inv_s * sod_m + k_h2o2_loss_s
    )
    values = np.real(_rk4(rhs, np.array((superoxide0_m, 0.0)), times, scale))
    return {
        "time_s": times,
        "superoxide_m": values[:, 0],
        "hydrogen_peroxide_m": values[:, 1],
    }


def unitary_embedding_consistency_check(
    parameters: EncounterParameters, duration_s: float
) -> float:
    """Algebraically check a 6x6 coherent unitary embedded in 8x8.

    This is not an independently implemented circuit. It covers only closed,
    coherent Hamiltonian evolution--never Lindblad relaxation, spin-selective
    reaction, escape, semiquinone/association stages, or downstream kinetics.
    """
    energies, vectors = np.linalg.eigh(hamiltonian(parameters))
    unitary6 = (vectors * np.exp(-1j * energies * duration_s)) @ vectors.conj().T
    unitary8 = np.eye(8, dtype=complex)
    physical = np.array([0, 1, 2, 4, 5, 6])
    unitary8[np.ix_(physical, physical)] = unitary6
    rng = np.random.default_rng(1729)
    state6 = rng.normal(size=6) + 1j * rng.normal(size=6)
    state6 /= np.linalg.norm(state6)
    state8 = np.zeros(8, complex)
    state8[physical] = state6
    return float(np.max(np.abs((unitary8 @ state8)[physical] - unitary6 @ state6)))


def independent_circuit_capability() -> dict[str, str | bool]:
    """Report whether the optional local Qiskit runtime can actually execute."""
    try:
        import qiskit  # noqa: F401
    except Exception as error:  # Import may fail at native-extension load time.
        return {"available": False, "reason": f"qiskit import failed: {error}"}
    return {
        "available": True,
        "reason": "qiskit imports; no independent circuit implementation is claimed",
    }


def _validate_parameter(name: str, parameter: Any) -> None:
    if not isinstance(parameter, dict) or not PARAMETER_FIELDS.issubset(parameter):
        present = set(parameter) if isinstance(parameter, dict) else set()
        raise ValueError(
            f"parameter {name!r} is missing schema fields: "
            f"{sorted(PARAMETER_FIELDS - present)}"
        )
    if parameter["status"] not in PARAMETER_STATUSES:
        raise ValueError(f"parameter {name!r} has invalid status {parameter['status']!r}")
    if not isinstance(parameter["enabled"], bool):
        raise ValueError(f"parameter {name!r} enabled must be boolean")
    if parameter["enabled"] and parameter["status"] not in ACTIVATABLE_STATUSES:
        raise ValueError(
            f"parameter {name!r} cannot be active with status {parameter['status']!r}"
        )
    if parameter["enabled"] and parameter["value"] is None:
        raise ValueError(f"enabled parameter {name!r} must have a value")


def parameter_value(parameter, disabled_value=None):
    """Return enabled evidence value or a caller-supplied algebraic placeholder."""
    _validate_parameter("value", parameter)
    return parameter["value"] if parameter["enabled"] else disabled_value


def resolve_parameter_block(block, disabled_values=None):
    disabled_values = disabled_values or {}
    resolved = {}
    for name, parameter in block.items():
        _validate_parameter(name, parameter)
        resolved[name] = parameter_value(parameter, disabled_values.get(name))
    return resolved


def validate_authority_bundle(json_path, provenance_path) -> dict[str, Any]:
    """Validate configuration, provenance ledger, and validation-data schemas."""
    raw = json.loads(Path(json_path).read_text())
    if raw.get("schema_version") != "2.1":
        raise ValueError("unsupported parameter schema_version")
    if raw.get("quantitative_prediction_supported") is not False:
        raise ValueError("quantitative_prediction_supported must be false")
    active = raw["active_model"]
    _validate_parameter("duration_s", active["duration_s"])
    for block_name in ("encounter", "downstream"):
        for name, parameter in active[block_name].items():
            _validate_parameter(f"{block_name}.{name}", parameter)
            if (
                not parameter["enabled"]
                and parameter["status"] == "unavailable"
                and parameter["value"] == 0
                and "not" not in parameter["limitations"].lower()
            ):
                raise ValueError(f"disabled zero {name!r} lacks a non-evidence warning")

    required_columns = {
        "topic", "symbol", "definition", "value_or_range", "unit", "species",
        "environment", "method", "source", "status", "limitations", "use_class",
        "code_location",
    }
    with Path(provenance_path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not required_columns.issubset(reader.fieldnames or []):
            raise ValueError("provenance CSV is missing required columns")
        rows = list(reader)
    symbols = [row["symbol"] for row in rows]
    if len(symbols) != len(set(symbols)):
        raise ValueError("provenance symbols must be unique")
    for row in rows:
        root = row["code_location"].split(".", 1)[0]
        if root not in raw:
            raise ValueError(f"provenance path has unknown root: {row['code_location']}")

    g_row = next(row for row in rows if row["symbol"] == "g_SQ")
    configured_g = active["encounter"]["g_radical"]
    if float(g_row["value_or_range"]) != configured_g["value"]:
        raise ValueError("active g_radical disagrees with provenance g_SQ")
    if "10.1016/0003-9861(91)90023-c" not in configured_g["source"].lower():
        raise ValueError("active g_radical source disagrees with provenance")

    for name in ("k_doublet_s", "k_quartet_s", "k_escape_s"):
        if active["encounter"][name]["enabled"]:
            raise ValueError(f"unavailable encounter rate {name!r} must remain disabled")
    for name, parameter in raw["measured_parameters"].items():
        if "oxygen" in name and "semiquinone" in name:
            if parameter["unit"] != "M^-1 s^-1":
                raise ValueError("bulk oxygen-reaction constants must remain bimolecular")

    required_validation = {"unit", "dose", "method", "source", "enabled"}
    for name, dataset in raw["validation_datasets"].items():
        missing = required_validation - set(dataset)
        if missing:
            raise ValueError(f"validation dataset {name!r} missing {sorted(missing)}")
        if dataset["enabled"] is not False:
            raise ValueError(f"validation dataset {name!r} must not drive predictions")
        if not any(key in dataset for key in ("uncertainty", "uncertainty_SE")):
            raise ValueError(f"validation dataset {name!r} lacks uncertainty metadata")
    return raw


def load_execution_context(
    json_path,
    mode: str = EVIDENCE_BACKED,
    sensitivity_scenario: str | None = None,
    allow_sensitivity: bool = False,
    provenance_path=None,
) -> ExecutionContext:
    """Resolve one of two policy-enforced execution modes."""
    raw = (
        validate_authority_bundle(json_path, provenance_path)
        if provenance_path is not None
        else json.loads(Path(json_path).read_text())
    )
    active = raw["active_model"]
    defaults = EncounterParameters()
    values = {
        name: getattr(defaults, name) for name in defaults.__dataclass_fields__
    }
    for name, parameter in active["encounter"].items():
        _validate_parameter(name, parameter)
        if parameter["enabled"]:
            values[name] = parameter["value"]

    duration = active["duration_s"]["value"] if active["duration_s"]["enabled"] else None
    if mode == EVIDENCE_BACKED:
        if sensitivity_scenario is not None or allow_sensitivity:
            raise ModelPolicyError("sensitivity options are invalid in evidence_backed mode")
    elif mode == SENSITIVITY:
        if not allow_sensitivity:
            raise ModelPolicyError("sensitivity mode requires explicit allow_sensitivity opt-in")
        if not sensitivity_scenario:
            raise ModelPolicyError("sensitivity mode requires a named sensitivity scenario")
        scenarios = raw.get("sensitivity_scenarios", {})
        if sensitivity_scenario not in scenarios:
            raise ModelPolicyError(
                f"unknown sensitivity scenario {sensitivity_scenario!r}; "
                f"choose one of {sorted(scenarios)}"
            )
        scenario = scenarios[sensitivity_scenario]
        if scenario.get("non_predictive") is not True:
            raise ModelPolicyError("sensitivity scenario must be marked non_predictive")
        duration = scenario["duration_s"]
        for name in scenario.get("use_config_values", []):
            if name not in active["encounter"]:
                raise ValueError(f"scenario references unknown encounter parameter {name!r}")
            values[name] = active["encounter"][name]["value"]
        for name, value in scenario.get("encounter_overrides", {}).items():
            if name not in values:
                raise ValueError(f"scenario overrides unknown encounter parameter {name!r}")
            values[name] = value
    else:
        raise ModelPolicyError("mode must be evidence_backed or sensitivity")

    for vector_name in ("field_t", "dipolar_axis", "effective_hyperfine_rad_s"):
        values[vector_name] = tuple(values[vector_name])
    return ExecutionContext(
        mode=mode,
        parameters=EncounterParameters(**values),
        duration_s=duration,
        authority=raw,
        sensitivity_scenario=sensitivity_scenario,
        non_predictive=mode == SENSITIVITY,
    )


def parameters_from_json(
    path,
    mode: str = EVIDENCE_BACKED,
    sensitivity_scenario: str | None = None,
    allow_sensitivity: bool = False,
):
    """Compatibility loader returning (parameters, active_model)."""
    context = load_execution_context(
        path, mode=mode, sensitivity_scenario=sensitivity_scenario,
        allow_sensitivity=allow_sensitivity,
    )
    return context.parameters, context.authority["active_model"]


def evidence_backed_encounter_prediction(*_args, **_kwargs):
    """Refuse encounter yields until k_D, k_Q, escape, and preparation exist."""
    raise ModelPolicyError(
        "evidence_backed mode refuses encounter-level yield predictions: "
        "k_D, k_Q, k_escape, encounter duration/association, and prepared D/Q state "
        "are unavailable"
    )


def illustrative_run_authorized(_active_model, cli_opt_in=False):
    """Deprecated compatibility policy: illustrative execution always needs opt-in."""
    return bool(cli_opt_in)
