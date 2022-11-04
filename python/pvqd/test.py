import numpy as np

from pvqd import PVQD
from ansatze import hweff_ansatz_adiab
from qiskit.circuit import ParameterVector, Parameter
from qiskit.opflow import PauliSumOp
from qiskit.algorithms.time_evolvers.time_evolution_problem import TimeEvolutionProblem

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
# time_param = Parameter('T')
# coeff_list = 2*[J/time*time_param] + 3*[h]
# hamiltonian = SparsePauliOp(["ZZI", "IZZ", "IIX", "IXI", "XII"], coeffs = np.array(coeff_list))
# hamiltonian = PauliSumOp(J*time_param/time*["ZZI", "IZZ"]) # + h*SparsePauliOp(["IIX", "IXI", "XII"]))
hamiltonian = PauliSumOp(SparsePauliOp(["IIX", "IXI", "XII"]), h)

# print(hamiltonian)

observable = Pauli("ZZI")
# ansatz = EfficientSU2(2, reps=1)

# optimizer = GradientDescent()

def grad_descent(loss_callable, initial_guess, gradient_callable, learn_rate=1.0, max_iter=100, threshold=1e-6):

    count = 0
    res_dict = {}

    loss = loss_callable(initial_guess)
    grad = gradient_callable(initial_guess)
    old_params = initial_guess

    losses = [loss]
    gradients = [grad]
    params = [old_params]

    while count<max_iter and loss>threshold:
        count += 1
        new_params = old_params - learn_rate*grad
        loss = loss_callable(new_params)
        grad = gradient_callable(new_params)

        losses.append(loss)
        gradients.append(grad)
        params.append(new_params)

        old_params = new_params

    res_dict["iter"] = count
    res_dict["fidelities"] = np.clip(1-np.array(losses), 0, 1)
    res_dict["gradients"] = gradients
    res_dict["parameters"] = params

    return res_dict

# setup the algorithm
pvqd = PVQD(
    fidelity,
    ansatz,
    initial_parameters,
    estimator,
    optimizer=grad_descent,
    num_timesteps=2,
)

# specify the evolution problem
problem = TimeEvolutionProblem(
    hamiltonian, time, aux_operators=[hamiltonian]
)

# and evolve!
result = pvqd.evolve(problem)
print(result)