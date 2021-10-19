# This is a test file for the time evolution of the Hamiltonian with parameter assignment method

import numpy as np
import json
import functools
import itertools
import matplotlib.pyplot as plt
from scipy   import  linalg as LA



import qiskit
from qiskit                               import Aer, execute
from qiskit.quantum_info 			      import Pauli
from qiskit.aqua                          import QuantumInstance
from qiskit.aqua.operators 			      import PauliOp, SummedOp, CircuitSampler, StateFn
from qiskit.circuit                       import ParameterVector
from qiskit.aqua.operators.evolutions     import Trotter, PauliTrotterEvolution

from qiskit.aqua.operators.state_fns      import CircuitStateFn
from qiskit.aqua.operators.expectations   import PauliExpectation, AerPauliExpectation, MatrixExpectation
from qiskit.aqua.operators.primitive_ops  import CircuitOp
from qiskit.aqua.operators                import Z, I



from pauli_function import *
from ansatze import *

shots = 1
backend  = Aer.get_backend('statevector_simulator')
instance = QuantumInstance(backend=backend,shots=shots)

nqubits = 3
time_step = 0.05

coup = 0.25
field = 1.0

hzz = generate_ising_Hzz(nqubits, coup)
hx = generate_ising_Hx(nqubits, field)

#initialize in construct_total_circuit

ft = ParameterVector('ft', 2)

vector = np.array([hzz, hx])

step_h = ft*vector*time_step

# print('tut')
# print(step_h)

# step_hzz  = ft[0]*time_step*hzz
# step_hx   = ft[1]*time_step*hx

trotter = PauliTrotterEvolution(reps=1)
# Uzz_dt   = trotter.convert(step_hzz.exp_i()).to_circuit()
# Ux_dt    = trotter.convert(step_hx.exp_i()).to_circuit()

# Uzz_dt   = trotter.convert(step_h[0].exp_i()).to_circuit()
# Ux_dt    = trotter.convert(step_h[1].exp_i()).to_circuit()

U_dt = np.sum([trotter.convert(step_h[j].exp_i()).to_circuit() for j in range(2)])
# print(U_dt_list)
#
# U_dt = np.sum(U_dt_list)

# print(U_dt)

tmax = 2.0
Nt = int(tmax/time_step)
time_func_zz = lambda x : x/tmax
time_func_x = lambda x : 1
# print(time_func_x(1.0))
time_evol_vec = np.array([[time_func_zz(time_step*i) for i in range(Nt+1)], [time_func_x(time_step*i) for i in range(Nt+1)]])

values_dict = [dict(zip(ft[:], time_evol_vec[:,i].tolist())) for i in range(Nt+1)]
# print(values_dict)

ansatz_p = ParameterVector('p',18)
ansatz = hweff_ansatz(ansatz_p)

expectation = PauliExpectation()

circ_wfn = ansatz + U_dt
circ_wfn = expectation.convert(circ_wfn)

print(circ_wfn)

t_wfn = circ_wfn.assign_parameters(values_dict[0])

print(t_wfn)

t_wfn = circ_wfn.assign_parameters(values_dict[1])

print(t_wfn)

# sampler = CircuitSampler(instance)



# U_dt = expectation.convert(U_dt)
# circ_wfn = StateFn(circ_wfn)


# for values in values_dict:
#     sampled_op = sampler.convert(U_dt,params=values)
#
#     print('tut')
#     print(sampled_op)

# tut = np.array([1, 2, 3, 4])
# print(tut[:2])
# print(tut[2:])
