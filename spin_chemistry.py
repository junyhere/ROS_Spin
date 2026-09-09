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
REACTION_SQ_O2_ET = "rxn.bulk.doxorubicin_semiquinone_plus_triplet_o2_et"
REACTION_SOD_O2MINUS = "rxn.downstream.cuznsod_plus_o2minus"
USE_LITERATURE_ARITHMETIC = "literature_arithmetic"
USE_CONDITION_MATCHED = "condition_matched_prediction"
USE_SENSITIVITY_ONLY = "sensitivity_only"
ACTIVATABLE_STATUSES = {"measured", "calculated", "fitted", "approved_surrogate"}
PARAMETER_STATUSES = ACTIVATABLE_STATUSES | {
    "assumed", "illustrative", "surrogate", "unavailable"
}
PARAMETER_FIELDS = {
    "value", "unit", "status", "enabled", "species", "environment",
    "temperature", "pH", "source", "limitations",
}
CONDITION_FIELDS = {
    "species", "protonation", "pH", "temperature", "solvent_buffer",
    "oxygen_conditions",
}


class ModelPolicyError(ValueError):
    """Raised when an execution request violates the evidence policy."""


def _finite_scalar(name: str, value: Any, *, nonnegative: bool = False) -> float:
    """Return a finite real scalar, rejecting booleans and array coercions."""
    if isinstance(value, (bool, np.bool_, str, bytes)) or not np.isscalar(value):
        raise ValueError(f"{name} must be a finite real scalar")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"{name} must be a finite real scalar") from error
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite")
    if nonnegative and result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _finite_vector(name: str, value: Any, length: int = 3) -> np.ndarray:
    try:
        result = np.asarray(value, dtype=float)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"{name} must contain {length} finite real values") from error
    if result.shape != (length,) or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain {length} finite real values")
    return result


def _sample_count(samples: Any) -> int:
    if isinstance(samples, (bool, np.bool_, float, np.floating)) or not isinstance(
        samples, (int, np.integer)
    ):
        raise ValueError("samples must be an integer >= 2")
    if int(samples) < 2:
        raise ValueError("samples must be an integer >= 2")
    return int(samples)


def _validate_encounter_parameters(parameters: "EncounterParameters") -> None:
    if not isinstance(parameters, EncounterParameters):
        raise ValueError("parameters must be EncounterParameters")
    for name in ("field_t", "dipolar_axis", "local_field_proxy_rad_s"):
        vector = _finite_vector(name, getattr(parameters, name))
        if name == "dipolar_axis" and parameters.dipolar_rad_s and not np.linalg.norm(vector):
            raise ValueError("dipolar_axis must be nonzero when dipolar coupling is active")
    for name in (
        "g_radical", "g_oxygen", "exchange_rad_s", "dipolar_rad_s",
        "oxygen_zfs_d_rad_s", "oxygen_zfs_e_rad_s",
    ):
        _finite_scalar(name, getattr(parameters, name))
    for name in (
        "radical_relaxation_s", "oxygen_relaxation_s", "k_doublet_s",
        "k_quartet_s", "k_escape_s",
    ):
        _finite_scalar(name, getattr(parameters, name), nonnegative=True)


def spin_matrices(spin: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    spin = _finite_scalar("spin", spin, nonnegative=True)
    if not np.isclose(2 * spin, round(2 * spin)):
        raise ValueError("spin must be a nonnegative integer or half-integer")
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
    local_field_proxy_rad_s: tuple[float, float, float] = (0.0, 0.0, 0.0)
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
    resolved_inputs: dict[str, Any] | None = None


def _unit(vector: tuple[float, float, float]) -> np.ndarray:
    value = _finite_vector("axis", vector)
    norm = np.linalg.norm(value)
    if norm == 0:
        raise ValueError("axis must be nonzero")
    return value / norm


def hamiltonian(parameters: EncounterParameters) -> np.ndarray:
    """Return H/hbar for the six-state electronic encounter space."""
    _validate_encounter_parameters(parameters)
    beta_e_over_hbar = 8.79410005e10  # rad s^-1 T^-1
    field = _finite_vector("field_t", parameters.field_t)
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
        parameters.local_field_proxy_rad_s[index] * R[index]
        for index in range(3)
    )
    return np.asarray(result, complex)


def initial_density(scenario: str, p_doublet: float | None = None) -> np.ndarray:
    """Return a benchmark manifold mixture, never a prepared pure wavefunction.

    ``doublet`` and ``quartet`` are maximally mixed states restricted to their
    respective manifolds. ``mixture`` requires 0 <= ``p_doublet`` <= 1.
    Chemically prepared states remain unsupported because no preparation
    evidence is available for this encounter.
    """
    if scenario == "unpolarized":
        if p_doublet is not None:
            raise ValueError("p_doublet is only valid for scenario='mixture'")
        return I6 / 6
    if scenario == "doublet":
        if p_doublet is not None:
            raise ValueError("p_doublet is only valid for scenario='mixture'")
        return P_DOUBLET / 2
    if scenario == "quartet":
        if p_doublet is not None:
            raise ValueError("p_doublet is only valid for scenario='mixture'")
        return P_QUARTET / 4
    if scenario == "mixture":
        if p_doublet is None:
            raise ValueError("scenario='mixture' requires p_doublet")
        p_doublet = _finite_scalar("p_doublet", p_doublet)
        if not 0 <= p_doublet <= 1:
            raise ValueError("p_doublet must be in [0, 1]")
        return p_doublet * P_DOUBLET / 2 + (1 - p_doublet) * P_QUARTET / 4
    raise ValueError("scenario must be unpolarized, doublet, quartet, or mixture")


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
    p_doublet: float | None = None,
) -> dict[str, Any]:
    """Propagate the full classical six-state density matrix by matrix exponential.

    d rho/dt = -i[H,rho] - {K + k_escape I,rho}/2 + L_relax(rho),
    K = k_D P_D + k_Q P_Q. Reaction and escape yields are time integrals
    of their respective fluxes. This is not an evidence-backed encounter
    prediction until every required encounter input is supported. The constant
    39-component linear system (36 density elements plus three accumulated
    yields) is advanced with a scaling-and-squaring Padé matrix exponential.
    """
    _validate_encounter_parameters(parameters)
    duration_s = _finite_scalar("duration_s", duration_s)
    if duration_s <= 0:
        raise ValueError("duration_s must be positive")
    samples = _sample_count(samples)

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
        (initial_density(initial_state, p_doublet).reshape(-1), np.zeros(3, complex))
    )
    trajectory = _exact_linear_trajectory(rhs, initial, times)
    densities = trajectory[:, :36].reshape(-1, 6, 6)
    p_doublet_population = np.real(np.einsum("ij,tji->t", P_DOUBLET, densities))
    p_quartet = np.real(np.einsum("ij,tji->t", P_QUARTET, densities))
    survival = np.real(np.trace(densities, axis1=1, axis2=2))
    doublet_yield, quartet_yield, escape_yield = (
        float(value) for value in np.real(trajectory[-1, 36:39])
    )
    cumulative_yields = np.real(trajectory[:, 36:39])
    probability_balance = doublet_yield + quartet_yield + escape_yield + float(
        survival[-1]
    )
    if not all(
        np.all(np.isfinite(value))
        for value in (densities, p_doublet_population, p_quartet, survival, trajectory[:, 36:39])
    ):
        raise FloatingPointError("reference encounter result is nonfinite")
    if min(doublet_yield, quartet_yield, escape_yield, float(survival.min())) < -1e-10:
        raise FloatingPointError("reference encounter result has substantial negative probability")
    return {
        "time_s": times,
        "density_matrices": densities,
        "p_doublet": p_doublet_population,
        "p_quartet": p_quartet,
        "survival": survival,
        "doublet_reaction_yield": doublet_yield,
        "quartet_reaction_yield": quartet_yield,
        "cumulative_doublet_reaction_yield": cumulative_yields[:, 0],
        "cumulative_quartet_reaction_yield": cumulative_yields[:, 1],
        "cumulative_primary_superoxide_yield": (
            cumulative_yields[:, 0] + cumulative_yields[:, 1]
        ),
        "cumulative_escape_yield": cumulative_yields[:, 2],
        "primary_superoxide_yield": doublet_yield + quartet_yield,
        "superoxide_yield": doublet_yield + quartet_yield,
        "escape_yield": escape_yield,
        "unresolved_probability": float(survival[-1]),
        "probability_balance": probability_balance,
        "probability_balance_error": abs(probability_balance - 1.0),
        "initial_state": initial_state,
        "p_doublet_initial": (
            float(np.trace(P_DOUBLET @ initial_density(initial_state, p_doublet)).real)
        ),
        "numerical_method": "constant-generator matrix exponential Pade(13)",
        "duration_s": duration_s,
        "samples": samples,
    }


def _dormand_prince_trajectory(
    rhs,
    initial: np.ndarray,
    times: np.ndarray,
    *,
    atol: float,
    rtol: float,
    max_steps: int,
) -> tuple[np.ndarray, int]:
    """Adaptive Dormand--Prince 5(4), implemented independently of the reference."""
    atol = _finite_scalar("atol", atol)
    rtol = _finite_scalar("rtol", rtol)
    if atol <= 0 or rtol <= 0:
        raise ValueError("atol and rtol must be positive")
    if isinstance(max_steps, bool) or not isinstance(max_steps, (int, np.integer)):
        raise ValueError("max_steps must be a positive integer")
    if max_steps <= 0:
        raise ValueError("max_steps must be a positive integer")

    value = np.asarray(initial, complex).copy()
    output = [value.copy()]
    time = float(times[0])
    step = max((float(times[-1]) - time) / 200, np.finfo(float).eps)
    accepted = 0
    attempted = 0
    for target in times[1:]:
        target = float(target)
        while time < target:
            attempted += 1
            if attempted > max_steps:
                raise RuntimeError("independent Dormand-Prince solver exceeded max_steps")
            step = min(step, target - time)
            k1 = rhs(time, value)
            k2 = rhs(time + step / 5, value + step * k1 / 5)
            k3 = rhs(
                time + 3 * step / 10,
                value + step * (3 * k1 / 40 + 9 * k2 / 40),
            )
            k4 = rhs(
                time + 4 * step / 5,
                value + step * (44 * k1 / 45 - 56 * k2 / 15 + 32 * k3 / 9),
            )
            k5 = rhs(
                time + 8 * step / 9,
                value + step * (
                    19372 * k1 / 6561 - 25360 * k2 / 2187
                    + 64448 * k3 / 6561 - 212 * k4 / 729
                ),
            )
            k6 = rhs(
                time + step,
                value + step * (
                    9017 * k1 / 3168 - 355 * k2 / 33
                    + 46732 * k3 / 5247 + 49 * k4 / 176
                    - 5103 * k5 / 18656
                ),
            )
            fifth = value + step * (
                35 * k1 / 384 + 500 * k3 / 1113 + 125 * k4 / 192
                - 2187 * k5 / 6784 + 11 * k6 / 84
            )
            k7 = rhs(time + step, fifth)
            fourth = value + step * (
                5179 * k1 / 57600 + 7571 * k3 / 16695 + 393 * k4 / 640
                - 92097 * k5 / 339200 + 187 * k6 / 2100 + k7 / 40
            )
            scale = atol + rtol * np.maximum(np.abs(value), np.abs(fifth))
            error = float(np.max(np.abs(fifth - fourth) / scale))
            if not np.isfinite(error):
                raise FloatingPointError("independent solver produced a nonfinite error")
            if error <= 1:
                value = fifth
                time += step
                accepted += 1
                if not np.all(np.isfinite(value)):
                    raise FloatingPointError("independent solver produced nonfinite values")
            factor = 5.0 if error == 0 else min(5.0, max(0.1, 0.9 * error ** -0.2))
            step *= factor
            if step <= np.finfo(float).eps * max(1.0, abs(time)):
                raise RuntimeError("independent solver step size underflow")
        output.append(value.copy())
    return np.asarray(output), accepted


def propagate_encounter_independent(
    parameters: EncounterParameters,
    duration_s: float,
    initial_state: str = "unpolarized",
    samples: int = 101,
    p_doublet: float | None = None,
    *,
    atol: float = 1e-11,
    rtol: float = 1e-9,
    max_steps: int = 200000,
) -> dict[str, Any]:
    """Independently integrate the encounter ODE with adaptive RK5(4).

    This routine constructs the density-matrix derivative directly and never
    forms, exponentiates, or calls the 39x39 generator used by the reference.
    """
    _validate_encounter_parameters(parameters)
    duration_s = _finite_scalar("duration_s", duration_s)
    if duration_s <= 0:
        raise ValueError("duration_s must be positive")
    samples = _sample_count(samples)
    initial = initial_density(initial_state, p_doublet)
    coherent = hamiltonian(parameters)
    reaction = parameters.k_doublet_s * P_DOUBLET + parameters.k_quartet_s * P_QUARTET
    escape = parameters.k_escape_s

    radical_ops = tuple(R) if parameters.radical_relaxation_s else ()
    oxygen_ops = tuple(O) if parameters.oxygen_relaxation_s else ()

    def independent_rhs(_time, state):
        density = state[:36].reshape((6, 6))
        derivative = -1j * (coherent @ density - density @ coherent)
        derivative -= 0.5 * (reaction @ density + density @ reaction)
        derivative -= escape * density
        for operator in radical_ops:
            derivative += parameters.radical_relaxation_s * (
                operator @ density @ operator
                - 0.5 * (operator @ operator @ density + density @ operator @ operator)
            )
        for operator in oxygen_ops:
            derivative += parameters.oxygen_relaxation_s * (
                operator @ density @ operator
                - 0.5 * (operator @ operator @ density + density @ operator @ operator)
            )
        flux_d = parameters.k_doublet_s * np.trace(P_DOUBLET @ density)
        flux_q = parameters.k_quartet_s * np.trace(P_QUARTET @ density)
        flux_e = escape * np.trace(density)
        return np.concatenate((derivative.ravel(), [flux_d, flux_q, flux_e]))

    times = np.linspace(0.0, duration_s, samples)
    trajectory, accepted_steps = _dormand_prince_trajectory(
        independent_rhs,
        np.concatenate((initial.ravel(), np.zeros(3, complex))),
        times,
        atol=atol,
        rtol=rtol,
        max_steps=max_steps,
    )
    densities = trajectory[:, :36].reshape((-1, 6, 6))
    p_doublet_values = np.einsum("ij,tji->t", P_DOUBLET, densities).real
    p_quartet_values = np.einsum("ij,tji->t", P_QUARTET, densities).real
    survival = np.trace(densities, axis1=1, axis2=2).real
    yields = trajectory[-1, 36:39].real
    cumulative_yields = trajectory[:, 36:39].real
    if not all(np.all(np.isfinite(item)) for item in (densities, yields, survival)):
        raise FloatingPointError("independent encounter result is nonfinite")
    if min(float(yields.min()), float(survival.min())) < -1e-8:
        raise FloatingPointError("independent encounter result has substantial negative probability")
    probability_balance = float(yields.sum() + survival[-1])
    return {
        "time_s": times,
        "density_matrices": densities,
        "p_doublet": p_doublet_values,
        "p_quartet": p_quartet_values,
        "survival": survival,
        "doublet_reaction_yield": float(yields[0]),
        "quartet_reaction_yield": float(yields[1]),
        "cumulative_doublet_reaction_yield": cumulative_yields[:, 0],
        "cumulative_quartet_reaction_yield": cumulative_yields[:, 1],
        "cumulative_primary_superoxide_yield": (
            cumulative_yields[:, 0] + cumulative_yields[:, 1]
        ),
        "cumulative_escape_yield": cumulative_yields[:, 2],
        "primary_superoxide_yield": float(yields[0] + yields[1]),
        "superoxide_yield": float(yields[0] + yields[1]),
        "escape_yield": float(yields[2]),
        "unresolved_probability": float(survival[-1]),
        "probability_balance": probability_balance,
        "probability_balance_error": abs(probability_balance - 1.0),
        "initial_state": initial_state,
        "p_doublet_initial": float(np.trace(P_DOUBLET @ initial).real),
        "numerical_method": "independent adaptive Dormand-Prince RK5(4)",
        "accepted_steps": accepted_steps,
        "atol": atol,
        "rtol": rtol,
        "duration_s": duration_s,
        "samples": samples,
    }


# Backward-compatible name; the complete six-state implementation is reference.
propagate_encounter = propagate_encounter_reference


def form_semiquinone(quinone_m: float, formation_fraction: float) -> float:
    """Classical upstream preparation stage, kept outside the density matrix."""
    quinone_m = _finite_scalar("quinone_m", quinone_m, nonnegative=True)
    formation_fraction = _finite_scalar("formation_fraction", formation_fraction)
    if not 0 <= formation_fraction <= 1:
        raise ValueError("quinone_m must be nonnegative and formation_fraction in [0,1]")
    return quinone_m * formation_fraction


def associate_encounters(
    semiquinone_m: float, oxygen_m: float, association_fraction: float
) -> float:
    """Classical encounter-association stage, limited by both reactants."""
    semiquinone_m = _finite_scalar("semiquinone_m", semiquinone_m, nonnegative=True)
    oxygen_m = _finite_scalar("oxygen_m", oxygen_m, nonnegative=True)
    association_fraction = _finite_scalar("association_fraction", association_fraction)
    if not 0 <= association_fraction <= 1:
        raise ValueError("concentrations must be nonnegative and association_fraction in [0,1]")
    return min(semiquinone_m, oxygen_m) * association_fraction


def primary_superoxide_formation(encounter_m: float, reaction_yield: float) -> float:
    """Map reacted encounters to primary superoxide at 1:1 stoichiometry."""
    encounter_m = _finite_scalar("encounter_m", encounter_m, nonnegative=True)
    reaction_yield = _finite_scalar("reaction_yield", reaction_yield)
    if not 0 <= reaction_yield <= 1:
        raise ValueError("encounter_m must be nonnegative and reaction_yield in [0,1]")
    return encounter_m * reaction_yield


def bulk_superoxide_formation_rate(
    semiquinone_m: float, oxygen_m: float, k_bulk_m_inv_s: float
) -> float:
    """Arithmetic for r=k_bulk[SQ][O2]; evidence policy is enforced separately."""
    semiquinone_m = _finite_scalar("semiquinone_m", semiquinone_m, nonnegative=True)
    oxygen_m = _finite_scalar("oxygen_m", oxygen_m, nonnegative=True)
    k_bulk_m_inv_s = _finite_scalar(
        "k_bulk_m_inv_s", k_bulk_m_inv_s, nonnegative=True
    )
    return k_bulk_m_inv_s * semiquinone_m * oxygen_m


def _condition_match(record: dict[str, Any], profile: dict[str, Any]) -> list[str]:
    """Return exact missing/mismatched condition fields without invented tolerances."""
    if not isinstance(profile, dict):
        raise ModelPolicyError("condition-matched use requires an explicit condition profile")
    missing = CONDITION_FIELDS - set(profile)
    if missing:
        raise ModelPolicyError(
            f"condition profile is missing required fields: {sorted(missing)}"
        )
    mismatches = []
    expected = record.get("conditions", {})
    for field in sorted(CONDITION_FIELDS):
        if profile[field] is None:
            mismatches.append(f"{field}: requested condition is unknown")
        elif expected.get(field) is None:
            mismatches.append(f"{field}: literature condition is unknown")
        elif profile[field] != expected[field]:
            mismatches.append(
                f"{field}: requested {profile[field]!r} != literature {expected[field]!r}"
            )
    return mismatches


def evaluate_bulk_superoxide_rate(
    authority: dict[str, Any],
    parameter_id: str,
    semiquinone_m: float,
    oxygen_m: float,
    *,
    use_scope: str,
    condition_profile: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Apply reaction identity and evidence/condition policy to a bulk rate law."""
    records = authority.get("parameter_records", {})
    if parameter_id not in records:
        raise ModelPolicyError(f"unknown bulk parameter ID {parameter_id!r}")
    record = records[parameter_id]
    if record.get("parameter_id") != parameter_id:
        raise ModelPolicyError("parameter key and stable parameter_id disagree")
    if record.get("reaction_id") != REACTION_SQ_O2_ET:
        raise ModelPolicyError(
            f"parameter {parameter_id!r} is for reaction {record.get('reaction_id')!r}, "
            f"not {REACTION_SQ_O2_ET!r}"
        )
    if record.get("unit") != "M^-1 s^-1":
        raise ModelPolicyError("bulk SQ/O2 parameter must have unit M^-1 s^-1")
    if use_scope not in record.get("allowed_uses", []):
        raise ModelPolicyError(
            f"parameter {parameter_id!r} does not permit use scope {use_scope!r}"
        )
    if record.get("status") not in {"measured", "calculated", "fitted"}:
        raise ModelPolicyError("selected bulk rate is not literature-backed")
    mismatches: list[str] = []
    if use_scope == USE_CONDITION_MATCHED:
        mismatches = _condition_match(record, condition_profile)
        if mismatches:
            raise ModelPolicyError(
                "condition-matched calculation refused: " + "; ".join(mismatches)
            )
        interpretation = "condition-matched bulk-law reference calculation"
    elif use_scope == USE_LITERATURE_ARITHMETIC:
        interpretation = (
            "literature-backed arithmetic at user-supplied concentrations; "
            "not a prediction and sensitivity-only outside the reported conditions"
        )
    elif use_scope == USE_SENSITIVITY_ONLY:
        interpretation = "sensitivity-only bulk-law arithmetic"
    else:
        raise ModelPolicyError(f"unsupported bulk use scope {use_scope!r}")
    return {
        "rate_m_s": bulk_superoxide_formation_rate(
            semiquinone_m, oxygen_m, record["value"]
        ),
        "interpretation": interpretation,
        "predictive_status": (
            "condition_matched_bulk_reference_only"
            if use_scope == USE_CONDITION_MATCHED
            else "not_a_prediction"
        ),
        "use_scope": use_scope,
        "parameter": {
            key: record.get(key)
            for key in (
                "parameter_id", "reaction_id", "value", "unit", "source",
                "conditions", "uncertainty", "limitations",
            )
        },
        "inputs": {
            "semiquinone_m": semiquinone_m,
            "triplet_oxygen_m": oxygen_m,
            "condition_profile": condition_profile,
        },
    }


def acid_base_fractions(pH: float, pKa_ho2: float) -> tuple[float, float]:
    pH = _finite_scalar("pH", pH)
    pKa_ho2 = _finite_scalar("pKa_ho2", pKa_ho2)
    exponent = np.log(10.0) * (pH - pKa_ho2)
    if exponent >= 0:
        exp_negative = np.exp(-exponent)
        fraction_ho2 = exp_negative / (1 + exp_negative)
    else:
        fraction_ho2 = 1 / (1 + np.exp(exponent))
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
    radical_pool_m = _finite_scalar("radical_pool_m", radical_pool_m, nonnegative=True)
    for name, value in (
        ("k_ho2_ho2_m_inv_s", k_ho2_ho2_m_inv_s),
        ("k_ho2_o2minus_m_inv_s", k_ho2_o2minus_m_inv_s),
    ):
        _finite_scalar(name, value, nonnegative=True)
    if k_o2minus_o2minus_m_inv_s is not None:
        _finite_scalar(
            "k_o2minus_o2minus_m_inv_s",
            k_o2minus_o2minus_m_inv_s,
            nonnegative=True,
        )
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
    radical_pool_m = _finite_scalar("radical_pool_m", radical_pool_m, nonnegative=True)
    sod_m = _finite_scalar("sod_m", sod_m, nonnegative=True)
    k_sod_m_inv_s = _finite_scalar("k_sod_m_inv_s", k_sod_m_inv_s, nonnegative=True)
    _, fraction_o2minus = acid_base_fractions(pH, pKa_ho2)
    return k_sod_m_inv_s * sod_m * fraction_o2minus * radical_pool_m


def h2o2_loss_rate(h2o2_m: float, k_h2o2_loss_s: float | None) -> float:
    """Compartmental first-order H2O2 loss; None means omitted, not zero evidence."""
    h2o2_m = _finite_scalar("h2o2_m", h2o2_m, nonnegative=True)
    if k_h2o2_loss_s is None:
        return 0.0
    return _finite_scalar("k_h2o2_loss_s", k_h2o2_loss_s, nonnegative=True) * h2o2_m


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
    """Solve a constant-pH pulse decay with an analytic radical-pool solution.

    Rapid HO2/O2-minus acid-base equilibrium, fixed pH, constant SOD, constant
    elementary rate coefficients, and first-order H2O2 loss are assumed. The
    spontaneous terms are reaction-event rates, while the radical-pool loss is
    twice their sum. SOD's radical species-loss rate produces half as many H2O2
    events. This is a single-pulse decay, not a continuously driven ROS source.
    """
    superoxide0_m = _finite_scalar("superoxide0_m", superoxide0_m, nonnegative=True)
    duration_s = _finite_scalar("duration_s", duration_s, nonnegative=True)
    pH = _finite_scalar("pH", pH)
    pKa_ho2 = _finite_scalar("pKa_ho2", pKa_ho2)
    k_ho2_ho2_m_inv_s = _finite_scalar(
        "k_ho2_ho2_m_inv_s", k_ho2_ho2_m_inv_s, nonnegative=True
    )
    k_ho2_o2minus_m_inv_s = _finite_scalar(
        "k_ho2_o2minus_m_inv_s", k_ho2_o2minus_m_inv_s, nonnegative=True
    )
    sod_m = _finite_scalar("sod_m", sod_m, nonnegative=True)
    if k_sod_m_inv_s is not None:
        k_sod_m_inv_s = _finite_scalar(
            "k_sod_m_inv_s", k_sod_m_inv_s, nonnegative=True
        )
    if k_o2minus_o2minus_m_inv_s is not None:
        k_o2minus_o2minus_m_inv_s = _finite_scalar(
            "k_o2minus_o2minus_m_inv_s",
            k_o2minus_o2minus_m_inv_s,
            nonnegative=True,
        )
    if k_h2o2_loss_s is not None:
        k_h2o2_loss_s = _finite_scalar(
            "k_h2o2_loss_s", k_h2o2_loss_s, nonnegative=True
        )
    samples = _sample_count(samples)

    fraction_ho2, fraction_o2minus = acid_base_fractions(pH, pKa_ho2)
    event_coefficient = (
        k_ho2_ho2_m_inv_s * fraction_ho2**2
        + k_ho2_o2minus_m_inv_s * fraction_ho2 * fraction_o2minus
        + (k_o2minus_o2minus_m_inv_s or 0.0) * fraction_o2minus**2
    )
    quadratic_species_loss = 2 * event_coefficient
    linear_species_loss = (k_sod_m_inv_s or 0.0) * sod_m * fraction_o2minus

    times = np.linspace(0.0, duration_s, samples)

    def radical_pool(time):
        if superoxide0_m == 0:
            return np.zeros_like(np.asarray(time, dtype=float))
        time = np.asarray(time, dtype=float)
        if linear_species_loss == 0:
            return superoxide0_m / (1 + quadratic_species_loss * superoxide0_m * time)
        exp_term = np.exp(-linear_species_loss * time)
        denominator = linear_species_loss + (
            quadratic_species_loss * superoxide0_m * (1 - exp_term)
        )
        return superoxide0_m * linear_species_loss * exp_term / denominator

    radical_values = np.asarray(radical_pool(times), float)
    h2o2 = np.zeros(samples, float)
    loss_constant = k_h2o2_loss_s or 0.0

    # Exact integrating factor for H2O2 loss plus adaptive Simpson quadrature
    # for the analytic radical-consumption forcing F(t) = -0.5 dC/dt.
    def forcing(time):
        pool = float(radical_pool(time))
        return event_coefficient * pool**2 + 0.5 * linear_species_loss * pool

    def adaptive_simpson(function, left, right, tolerance, depth=20):
        middle = (left + right) / 2
        f_left, f_middle, f_right = function(left), function(middle), function(right)
        whole = (right - left) * (f_left + 4 * f_middle + f_right) / 6

        def refine(a, b, fa, fm, fb, estimate, tol, remaining):
            m = (a + b) / 2
            lm, rm = (a + m) / 2, (m + b) / 2
            flm, frm = function(lm), function(rm)
            left_est = (m - a) * (fa + 4 * flm + fm) / 6
            right_est = (b - m) * (fm + 4 * frm + fb) / 6
            correction = left_est + right_est - estimate
            if remaining == 0:
                if abs(correction) > 15 * tol:
                    raise RuntimeError("downstream H2O2 quadrature did not converge")
                return left_est + right_est + correction / 15
            if abs(correction) <= 15 * tol:
                return left_est + right_est + correction / 15
            return refine(a, m, fa, flm, fm, left_est, tol / 2, remaining - 1) + refine(
                m, b, fm, frm, fb, right_est, tol / 2, remaining - 1
            )

        return refine(
            left, right, f_left, f_middle, f_right, whole, tolerance, depth
        )

    if loss_constant == 0:
        h2o2 = 0.5 * (superoxide0_m - radical_values)
    else:
        for index in range(1, samples):
            left, right = float(times[index - 1]), float(times[index])
            interval = right - left
            decay = np.exp(-loss_constant * interval)
            integrand = lambda time, end=right: np.exp(
                -loss_constant * (end - time)
            ) * forcing(time)
            tolerance = max(1e-18, 1e-11 * max(superoxide0_m, 1e-12))
            h2o2[index] = decay * h2o2[index - 1] + adaptive_simpson(
                integrand, left, right, tolerance
            )

    accumulated_h2o2_loss = 0.5 * (superoxide0_m - radical_values) - h2o2
    tolerance = 5e-12 * max(superoxide0_m, 1e-12)
    for name, values in (
        ("radical pool", radical_values),
        ("H2O2", h2o2),
        ("accumulated H2O2 loss", accumulated_h2o2_loss),
    ):
        if not np.all(np.isfinite(values)):
            raise FloatingPointError(f"{name} trajectory is nonfinite")
        if float(np.min(values)) < -tolerance:
            raise FloatingPointError(f"{name} trajectory has a substantial negative value")
        values[np.abs(values) <= tolerance] = 0.0

    radical_balance = radical_values + 2 * h2o2 + 2 * accumulated_h2o2_loss
    balance_error = float(np.max(np.abs(radical_balance - superoxide0_m)))
    return {
        "time_s": times,
        "radical_pool_m": radical_values,
        "ho2_m": fraction_ho2 * radical_values,
        "o2minus_m": fraction_o2minus * radical_values,
        "hydrogen_peroxide_m": h2o2,
        "accumulated_hydrogen_peroxide_loss_m": accumulated_h2o2_loss,
        "radical_equivalent_balance_m": radical_balance,
        "max_radical_equivalent_balance_error_m": balance_error,
        "event_rate_convention": (
            "spontaneous terms are H2O2 reaction-event rates; radical-pool species "
            "loss is twice the event rate; SOD species loss yields 0.5 H2O2"
        ),
        "approximations": [
            "fixed pH",
            "rapid HO2/O2-minus acid-base equilibrium",
            "constant SOD concentration",
            "constant elementary rate coefficients",
            "single-pulse radical decay with no continuous source",
            "first-order H2O2 loss when enabled",
        ],
        "numerical_method": (
            "analytic Riccati radical-pool solution plus integrating-factor "
            "adaptive quadrature for H2O2"
        ),
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


def unitary_embedding_consistency_check(
    parameters: EncounterParameters, duration_s: float
) -> float:
    """Algebraically check a 6x6 coherent unitary embedded in 8x8.

    This is not an independently implemented circuit. It covers only closed,
    coherent Hamiltonian evolution--never Lindblad relaxation, spin-selective
    reaction, escape, semiquinone/association stages, or downstream kinetics.
    """
    duration_s = _finite_scalar("duration_s", duration_s, nonnegative=True)
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
        "reason": "qiskit imports and the executable circuit validator is available",
    }


def execute_three_qubit_unitary_circuit(
    parameters: EncounterParameters,
    duration_s: float,
    *,
    seed: int = 1729,
) -> dict[str, Any]:
    """Execute the coherent 6-state propagator as a dense 3-qubit unitary.

    Qiskit's little-endian statevector indices are checked explicitly. This is
    simulator/embedding consistency using one inserted dense unitary; it is not
    an independent Hamiltonian derivation and does not cover open-system terms.
    """
    _validate_encounter_parameters(parameters)
    duration_s = _finite_scalar("duration_s", duration_s, nonnegative=True)
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)):
        raise ValueError("seed must be an integer")
    try:
        from qiskit import QuantumCircuit
        from qiskit.circuit.library import UnitaryGate
        from qiskit.quantum_info import Statevector
    except Exception as error:
        raise RuntimeError(f"qiskit circuit execution unavailable: {error}") from error

    physical = np.array([0, 1, 2, 4, 5, 6])
    basis_order = ["000", "001", "010", "100", "101", "110"]
    ordering_error = 0.0
    for index, label in zip(physical, basis_order):
        prepared = Statevector.from_label(label)
        ordering_error = max(
            ordering_error,
            float(np.max(np.abs(prepared.data - np.eye(8, dtype=complex)[index]))),
        )

    energies, vectors = np.linalg.eigh(hamiltonian(parameters))
    unitary6 = (vectors * np.exp(-1j * energies * duration_s)) @ vectors.conj().T
    unitary8 = np.eye(8, dtype=complex)
    unitary8[np.ix_(physical, physical)] = unitary6

    rng = np.random.default_rng(int(seed))
    state6 = rng.normal(size=6) + 1j * rng.normal(size=6)
    state6 /= np.linalg.norm(state6)
    state8 = np.zeros(8, complex)
    state8[physical] = state6
    circuit = QuantumCircuit(3)
    circuit.append(UnitaryGate(unitary8, label="exp(-iH6t) embedded"), [0, 1, 2])
    simulated = Statevector(state8).evolve(circuit).data
    expected6 = unitary6 @ state6
    complement = np.array([3, 7])

    pd8 = np.zeros((8, 8), complex)
    pq8 = np.zeros((8, 8), complex)
    pd8[np.ix_(physical, physical)] = P_DOUBLET
    pq8[np.ix_(physical, physical)] = P_QUARTET
    expected8 = np.zeros(8, complex)
    expected8[physical] = expected6
    doublet_observable_error = abs(
        np.vdot(simulated, pd8 @ simulated) - np.vdot(expected8, pd8 @ expected8)
    )
    quartet_observable_error = abs(
        np.vdot(simulated, pq8 @ simulated) - np.vdot(expected8, pq8 @ expected8)
    )
    dq_error = max(doublet_observable_error, quartet_observable_error)
    return {
        "status": "executed",
        "qubits": 3,
        "basis_order": basis_order,
        "physical_indices": physical.tolist(),
        "gate_implementation": "one inserted dense 8x8 UnitaryGate; not decomposed gates",
        "max_basis_ordering_error": ordering_error,
        "max_statevector_error": float(
            np.max(np.abs(simulated[physical] - expected6))
        ),
        "doublet_observable_error": float(abs(doublet_observable_error)),
        "quartet_observable_error": float(abs(quartet_observable_error)),
        "physical_subspace_leakage_probability": float(
            np.sum(np.abs(simulated[complement]) ** 2)
        ),
        "max_doublet_quartet_observable_error": float(abs(dq_error)),
        "scope": "closed coherent unitary simulator/embedding consistency only",
        "does_not_validate": [
            "independent Hamiltonian construction",
            "spin-selective reaction",
            "escape",
            "Lindblad relaxation",
            "upstream or downstream chemistry",
        ],
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
    value = parameter["value"]
    if value is not None:
        if isinstance(value, list):
            _finite_vector(f"parameter {name!r} value", value)
        elif isinstance(value, (int, float)) and not isinstance(value, bool):
            _finite_scalar(f"parameter {name!r} value", value)


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


def _resolve_json_path(document: dict[str, Any], path: str) -> Any:
    """Resolve every component of a dot-separated provenance object path."""
    current: Any = document
    traversed = []
    for component in path.split("."):
        traversed.append(component)
        if not isinstance(current, dict) or component not in current:
            raise ValueError(
                f"broken provenance path {path!r} at {'.'.join(traversed)!r}"
            )
        current = current[component]
    return current


def _source_token(source: Any) -> str | None:
    if source is None:
        return None
    text = str(source).lower().replace("https://doi.org/", "doi:")
    return text.split("doi:", 1)[1].strip() if "doi:" in text else text.strip()


def validate_authority_bundle(json_path, provenance_path) -> dict[str, Any]:
    """Validate configuration, provenance ledger, and validation-data schemas."""
    raw = json.loads(Path(json_path).read_text())
    if raw.get("schema_version") != "2.2":
        raise ValueError("unsupported parameter schema_version")
    if raw.get("quantitative_prediction_supported") is not False:
        raise ValueError("quantitative_prediction_supported must be false")
    active = raw["active_model"]
    _validate_parameter("duration_s", active["duration_s"])
    for block_name in ("encounter", "downstream"):
        for name, parameter in active[block_name].items():
            _validate_parameter(f"{block_name}.{name}", parameter)
            if name.startswith("k_") or name.endswith("relaxation_s"):
                for field in ("parameter_id", "reaction_id", "allowed_uses"):
                    if field not in parameter:
                        raise ValueError(
                            f"kinetic parameter {block_name}.{name} lacks {field}"
                        )
                if not isinstance(parameter["allowed_uses"], list):
                    raise ValueError(
                        f"kinetic parameter {block_name}.{name} allowed_uses must be a list"
                    )
            if (
                not parameter["enabled"]
                and parameter["status"] == "unavailable"
                and parameter["value"] == 0
                and "not" not in parameter["limitations"].lower()
            ):
                raise ValueError(f"disabled zero {name!r} lacks a non-evidence warning")

    required_columns = {
        "topic", "symbol", "definition", "value_or_range", "unit", "species",
        "charge_state", "environment", "temperature", "pH", "method", "source",
        "status", "uncertainty", "limitations", "use_class", "code_location",
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
        target = _resolve_json_path(raw, row["code_location"])
        if not isinstance(target, dict):
            raise ValueError(
                f"provenance path for {row['symbol']} must reference an exact JSON object"
            )
        if isinstance(target, dict):
            if (
                "unit" in target
                and row["unit"].lower() not in {"not available", "na"}
                and target["unit"] != row["unit"]
            ):
                raise ValueError(
                    f"provenance unit mismatch for {row['symbol']}: "
                    f"{row['unit']!r} != {target['unit']!r}"
                )
            if "source" in target:
                left, right = _source_token(row["source"]), _source_token(target["source"])
                if left not in {None, "na"} and right not in {None, "na"} and left != right:
                    raise ValueError(f"provenance source mismatch for {row['symbol']}")
            if "value" in target:
                if "status" in target and target["status"] != row["status"]:
                    raise ValueError(f"provenance status mismatch for {row['symbol']}")
                try:
                    csv_value = float(row["value_or_range"])
                except ValueError:
                    csv_value = None
                if csv_value is not None and target["value"] is not None:
                    if not np.isclose(csv_value, float(target["value"]), rtol=1e-12, atol=0):
                        raise ValueError(f"provenance value mismatch for {row['symbol']}")

    g_row = next(row for row in rows if row["symbol"] == "g_SQ")
    configured_g = active["encounter"]["g_radical"]
    if float(g_row["value_or_range"]) != configured_g["value"]:
        raise ValueError("active g_radical disagrees with provenance g_SQ")
    if "10.1016/0003-9861(91)90023-c" not in configured_g["source"].lower():
        raise ValueError("active g_radical source disagrees with provenance")

    for name in ("k_doublet_s", "k_quartet_s", "k_escape_s"):
        if active["encounter"][name]["enabled"]:
            raise ValueError(f"unavailable encounter rate {name!r} must remain disabled")
    records = raw.get("parameter_records", {})
    if not records:
        raise ValueError("parameter_records must not be empty")
    for name, parameter in records.items():
        if parameter.get("parameter_id") != name:
            raise ValueError(f"parameter record key/ID mismatch for {name!r}")
        if not isinstance(parameter.get("allowed_uses"), list):
            raise ValueError(f"parameter record {name!r} lacks allowed_uses")
        if parameter.get("reaction_id") == REACTION_SQ_O2_ET:
            if parameter["unit"] != "M^-1 s^-1":
                raise ValueError("bulk oxygen-reaction constants must remain bimolecular")
            if not CONDITION_FIELDS.issubset(parameter.get("conditions", {})):
                raise ValueError(f"bulk record {name!r} lacks complete conditions")
        if USE_CONDITION_MATCHED in parameter.get("allowed_uses", []):
            unknown = [
                key for key in CONDITION_FIELDS if parameter["conditions"].get(key) is None
            ]
            if unknown:
                raise ValueError(
                    f"condition-matched record {name!r} has unknown conditions: {unknown}"
                )

    for profile_name, profile in raw.get("condition_profiles", {}).items():
        missing = CONDITION_FIELDS - set(profile)
        if missing:
            raise ValueError(
                f"condition profile {profile_name!r} missing fields {sorted(missing)}"
            )

    required_validation = {
        "dataset_id", "unit", "dose", "method", "source", "enabled"
    }
    for name, dataset in raw["validation_datasets"].items():
        missing = required_validation - set(dataset)
        if missing:
            raise ValueError(f"validation dataset {name!r} missing {sorted(missing)}")
        if dataset["enabled"] is not False:
            raise ValueError(f"validation dataset {name!r} must not drive predictions")
        if dataset["dataset_id"] != name:
            raise ValueError(f"validation dataset key/ID mismatch for {name!r}")
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
    json_path = Path(json_path)
    provenance_path = (
        Path(provenance_path)
        if provenance_path is not None
        else json_path.with_name("parameter_provenance.csv")
    )
    raw = validate_authority_bundle(json_path, provenance_path)
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

    for vector_name in ("field_t", "dipolar_axis", "local_field_proxy_rad_s"):
        values[vector_name] = tuple(values[vector_name])
    resolved_parameters = EncounterParameters(**values)
    _validate_encounter_parameters(resolved_parameters)
    if duration is not None:
        duration = _finite_scalar("duration_s", duration)
        if duration <= 0:
            raise ValueError("duration_s must be positive when enabled")
    return ExecutionContext(
        mode=mode,
        parameters=resolved_parameters,
        duration_s=duration,
        authority=raw,
        sensitivity_scenario=sensitivity_scenario,
        non_predictive=mode == SENSITIVITY,
        resolved_inputs={
            name: {
                "value": values[name],
                "unit": active["encounter"].get(name, {}).get("unit"),
                "origin": (
                    "sensitivity_override"
                    if mode == SENSITIVITY
                    and name in raw["sensitivity_scenarios"][sensitivity_scenario].get(
                        "encounter_overrides", {}
                    )
                    else "configured_value"
                    if active["encounter"].get(name, {}).get("enabled")
                    or (
                        mode == SENSITIVITY
                        and name in raw["sensitivity_scenarios"][sensitivity_scenario].get(
                            "use_config_values", []
                        )
                    )
                    else "algebraic_disabled_default"
                ),
            }
            for name in values
        },
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
