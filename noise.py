from qiskit_aer.noise import (NoiseModel, amplitude_damping_error, 
                                        phase_damping_error)


def make_noise(gamma_amp: float, phi_frac: float = 0.0) -> NoiseModel:
    gamma_amp = gamma_amp * (1.0 - phi_frac)
    gamma_phi = gamma_amp * phi_frac
    noise = NoiseModel()
    if gamma_amp > 0:
        noise.add_all_qubit_quantum_error(amplitude_damping_error(gamma_amp), ["id"])
    if gamma_phi > 0:
        noise.add_all_qubit_quantum_error(phase_damping_error(gamma_phi), ["id"])
    return noise

