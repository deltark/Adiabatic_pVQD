import numpy as np

from pvqd import PVQD
from ansatze import hweff_ansatz_adiab
from qiskit.circuit import ParameterVector, Parameter

from qiskit.algorithms.state_fidelities import ComputeUncompute
from qiskit.primitives import Estimator, Sampler
from qiskit import BasicAer
# from qiskit.circuit.library import EfficientSU2
from qiskit.quantum_info import Pauli, SparsePauliOp
from qiskit.algorithms.optimizers import GradientDescent

n_spins = 3
depth = 2
n_params = depth*n_spins + depth*(n_spins-1)

ansatz = hweff_ansatz_adiab(n_spins, depth, ParameterVector('p', n_params))
initial_parameters = np.zeros(n_params)

sampler = Sampler()
fidelity = ComputeUncompute(sampler)
estimator = Estimator()

time = 1.0
J = -1
h = -1

hamiltonian = J/time * np.multiply(Parameter('t'), SparsePauliOp(Pauli("ZZI")), Pauli("IZZ")) + h * SparsePauliOp([Pauli("IIX"), Pauli("IXI"), Pauli("XII")])

observable = Pauli("ZZ")
# ansatz = EfficientSU2(2, reps=1)

optimizer = GradientDescent()

# setup the algorithm
pvqd = PVQD(
    fidelity,
    ansatz,
    estimator,
    initial_parameters,
    num_timesteps=2,
    optimizer=optimizer,
)

# specify the evolution problem
problem = EvolutionProblem(
    hamiltonian, time, aux_operators=[hamiltonian, observable]
)

# and evolve!
result = pvqd.evolve(problem)
print(result)