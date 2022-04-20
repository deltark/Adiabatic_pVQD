"""A sample runtime program that submits random circuits for user-specified iterations."""

from qiskit import Aer
from qiskit.providers.ibmq.runtime import UserMessenger
from qiskit.providers.ibmq.runtime.utils import RuntimeEncoder, RuntimeDecoder
from qiskit.providers.ibmq import RunnerResult
import random
import json
import time

# from qiskit import transpile
# from qiskit.circuit.random import random_circuit
# from qiskit.algorithms.optimizers import Optimizer
from qiskit.circuit import ParameterVector, QuantumCircuit
from qiskit.circuit.library.n_local import EfficientSU2

from qiskit import IBMQ
from torch import threshold

import runtime_pVQD

from qiskit.opflow import X, Z, I

from ansatze import *
from pauli_function import *

IBMQ.load_account()
provider = IBMQ.get_provider(
    hub='ibm-q-research-2', group='epfl-2', project='main')

# runtime_backends = provider.backends(input_allowed='runtime')
# print(f"Backends that support Qiskit Runtime: {runtime_backends}")
nqubits = 3
tmax = 3
dt = 0.05
NN = 1
maxiter = 100
shots = 1
# shots = 100000
threshold = 0.99999
# hzz = generate_ising_Hzz(nqubits, -1.0)
# hx = generate_ising_Hx(nqubits, -1.0)
optimizer = 'line search'
ham = ['hzz', 'hx']
anstz = "hweff"
depth = 1

# inital_point = json.load(open(
#     'data/interim_runtime/err_mitig_step10_ibmq_lima_NN1_iter20_shots8000.dat'))
initial_point = None

# h_tfunc = [lambda x: x/tmax]

# callback function for interim results
def interim_result_callback(job_id, interim_result):
    print(interim_result)
    print("callback time: ", time.localtime(time.time()))
    filename = ('data/interim_runtime/err_mitig_step'+str(interim_result["time_slice"][0])+'_'+backend.name() + '_NN'+str(NN) +
                '_iter'+str(maxiter)+'_shots'+str(shots)+'.dat')
    json.dump(interim_result, open(filename, 'w+'))


#pvqd inputs
inputs = {"nqubits": nqubits, "iterations": maxiter, "tmax": tmax, "dt": dt,
          "hamiltonian": ham, "NN": NN, "shots": shots, "ansatz": anstz, "initial_point": initial_point,
          "measurement_error_mitigation": False, "optimizer": optimizer, "threshold": threshold, "depth": depth}

#vqe inputs?
# def dumb_ansatz(n_spins, p):
#     circuit = QuantumCircuit(n_spins)
#     for i in range(n_spins):
#         circuit.rx(p[0], i)
#     return circuit

# params = ParameterVector('p', 12)
# ansatz = EfficientSU2(nqubits, reps=1)
# ham = (Z ^ Z ^ I) + (I ^ Z ^ Z) + (I ^ X ^ I) + (X ^ I ^ I) + (I ^ I ^ X)
# optimizer = {'name': 'QN-SPSA', 'maxiter': maxiter}
# inputs = {"ansatz": ansatz, "operator": ham, "optimizer": optimizer, "measurement_error_mitigation": True}

backend = Aer.get_backend('statevector_simulator')
# backend = Aer.get_backend('qasm_simulator')
user_messenger = UserMessenger()
serialized_inputs = json.dumps(inputs, cls=RuntimeEncoder)
deserialized_inputs = json.loads(serialized_inputs, cls=RuntimeDecoder)

result = runtime_pVQD.main(backend, user_messenger, **deserialized_inputs)
print(result)
print("\n")


# def interim_result_callback(job_id, interim_result):
#     print(f"interim result: {interim_result}")

################################# HARDWARE
# # print("on hardware:")
# # backend = provider.get_backend('ibmq_qasm_simulator')
# backend = provider.get_backend('ibmq_lima')
# # backend = provider.get_backend('ibm_lagos')
# # backend = provider.get_backend('simulator_statevector')
# # Configure backend options
# options = {'backend_name': backend.name()}
# program_id = "p-vqd-xL289veY54"

# # tic = time.time()
# # hour = time.localtime(tic)
# # print("Ding dong, il est ", hour)
# job = provider.runtime.run(program_id=program_id,
#                            options=options,
#                            inputs=inputs,
#                            callback=interim_result_callback)

# # Get runtime job result.
# # begin = time.time()
# result = job.result()
# # end = time.time()-begin
# # tac = time.time()
# # print("total time since submission: ", tac-tic)
# # hour = time.localtime(tac)
# # print("Ding dong, il est ", hour)

# # result['exec_time'] = end
# print(result)
# filename = ('data/test_vqe_lagos_QN-SPSA.dat')
# filename = ('data/VQD/runtime_'+backend.name()+'_NN'+str(NN)+'_iter'+str(maxiter)+'_shots'+str(shots)+'.dat')
# filename = ('data/VQD/runtimepVQD_full_errmitig_'+backend.name()+'.dat')
filename = ('data/VQD/linesearch_optsteps_'+backend.name()+'_T'+str(tmax)+'_dt'+str(dt)+'_shots'+str(shots)+'.dat')
json.dump(result, open(filename, 'w+'))



####TODO
#At every time step, take the best (lowest) configuration parameters instead of the last one
#maybe show different evolution between VQE and pVQD
#For plots: always t/T on abscisse, and 0 at ground energy (E-E_GS on ordonnée)
#Also try different optimisers on first time step and select the best one (see Stefano's message)
#add error mitigation
