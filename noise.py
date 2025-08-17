from qiskit_aer.noise import (NoiseModel, amplitude_damping_error, phase_damping_error)


def noise_mod(gamma_amp: float, phi_frac: float = 0.0) -> NoiseModel:
    """Create an amplitude and phase damping noise model"""
    gamma_amp = gamma_amp * (1.0 - phi_frac)
    gamma_phi = gamma_amp * phi_frac
    noise = NoiseModel()

    error = None
    if gamma_amp > 0:
        error = amplitude_damping_error(gamma_amp)
    if gamma_phi > 0:
        phi_err = phase_damping_error(gamma_phi)
        error = phi_err if error is None else error.compose(phi_err)
    if error:
        noise.add_all_qubit_quantum_error(error, ["id"])
    return noise
