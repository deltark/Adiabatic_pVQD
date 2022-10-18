import numpy as np
import scipy as sc

from qiskit import Aer
from pvqd_julien import PVQD
from qiskit.algorithms.optimizers import SPSA, L_BFGS_B, GradientDescent, COBYLA, optimizer
from qiskit.circuit import ParameterVector
from qiskit.circuit.library import EfficientSU2
from qiskit.opflow import X, Z, I, MatrixExpectation, Gradient
from qiskit.quantum_info import Statevector

from qiskit import IBMQ

import matplotlib.pyplot as plt

from ansatze import *

"""Test a simple evolution."""


def exact(final_time, dt, hamiltonian, observable, initial_state):
    """Get the exact values for energy and the observable."""
    energies = []  # list of energies evaluated at timesteps dt
    obs = []  # list of observables
    ts = []  # list of timepoints at which energy/obs are evaluated
    t = 0
    while t <= final_time:
        # get exact state at time t
        exact_state = initial_state.evolve(
            sc.linalg.expm(-1j * t * hamiltonian.to_matrix()))

        # store observables and time
        ts.append(t)
        energies.append(exact_state.expectation_value(
            hamiltonian.to_matrix()))
        obs.append(exact_state.expectation_value(observable.to_matrix()))

        # next timestep
        t += dt

    return ts, energies, obs

time = 3
dt = 0.06

def timedep_ham(t):
    return -t/time * ((Z ^ Z ^ I) + (I ^ Z ^ Z)) - ((I ^ X ^ I) + (X ^ I ^ I) + (I ^ I ^ X))

IBMQ.load_account()
provider = IBMQ.get_provider(
    hub='ibm-q-research-2', group='epfl-2', project='main')
# backend = Aer.get_backend("statevector_simulator")
backend = provider.get_backend('ibmq_qasm_simulator')
# backend = provider.get_backend('ibmq_manila')
expectation = MatrixExpectation()
hamiltonian = timedep_ham
observable = Z ^ Z ^ I

spins = 3
depth = 1

# ansatz = EfficientSU2(2, reps=1)
num_params = depth*spins + depth*(spins-1)
params = ParameterVector('p', num_params)
ansatz = hweff_ansatz_adiab(spins, depth, params)
# optimizer = SPSA(maxiter=300, learning_rate=0.1, perturbation=0.01)
optimizer = {'name': 'QN-SPSA', 'maxiter': 20}
# optimizer = SPSA()
# optimizer = COBYLA()
# optimizer = GradientDescent(learning_rate=0.01)
# optimizer = L_BFGS_B()
initial_parameters = np.zeros(ansatz.num_parameters)
# initial_parameters[-2] = np.pi / 2
# initial_parameters[-4] = np.pi / 2

# run pVQD keeping track of the energy and the magnetization
pvqd = PVQD(ansatz, initial_parameters, optimizer,
            quantum_instance=backend, expectation=expectation, provider=provider)
result = pvqd.evolve(hamiltonian, time, dt, observables=[hamiltonian(time), observable], step_vqe=True)

# get reference results
initial_state = Statevector(ansatz.bind_parameters(initial_parameters))
ref_t, ref_energy, ref_magn = exact(time, dt, hamiltonian(time), observable, initial_state)

_, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.set_title("Energy")
ax1.plot(result.times, result.observables[:, 0], label="pVQD")
ax1.plot(ref_t, ref_energy, label="exact")
ax2.set_title("Magnetization")
ax2.plot(result.times, result.observables[:, 1], label="pVQD")
ax2.plot(ref_t, ref_magn, label="exact")
ax2.set_xlabel(r"time $t$")
plt.tight_layout()
plt.show()
print(result)
