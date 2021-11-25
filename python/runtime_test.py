"""A sample runtime program that submits random circuits for user-specified iterations."""

from qiskit import Aer
from qiskit.providers.ibmq.runtime import UserMessenger
from qiskit.providers.ibmq.runtime.utils import RuntimeEncoder, RuntimeDecoder
from qiskit.providers.ibmq import RunnerResult
import random
import json
import time

from qiskit import transpile
from qiskit.circuit.random import random_circuit

from qiskit import IBMQ

import runtime_pVQD

# from pauli_function import *

IBMQ.load_account()
provider = IBMQ.get_provider(
    hub='ibm-q-research-2', group='epfl-2', project='main')

# runtime_backends = provider.backends(input_allowed='runtime')
# print(f"Backends that support Qiskit Runtime: {runtime_backends}")
nqubits = 3
tmax = 3.0
NN = 1
maxiter = 30
shots = 2000
# hzz = generate_ising_Hzz(nqubits, -1.0)
# hx = generate_ising_Hx(nqubits, -1.0)
ham = ['hzz', 'hx']

# h_tfunc = [lambda x: x/tmax]

# callback function for interim results
def interim_result_callback(job_id, interim_result):
    print(interim_result)
    filename = ('data/VQD/interimstep'+str(interim_result["t_step"])+'_runtime_'+backend.name() +
                '_NN'+str(NN)+'_iter'+str(maxiter)+'_shots'+str(shots)+'.dat')
    json.dump(interim_result, open(filename, 'w+'))


inputs = {"nqubits": nqubits, "iterations": maxiter, "tmax": tmax, "dt": tmax,
          "hamiltonian": ham, "NN": NN, "shots": shots}

# backend = Aer.get_backend('qasm_simulator')
user_messenger = UserMessenger()
# serialized_inputs = json.dumps(inputs, cls=RuntimeEncoder)
# deserialized_inputs = json.loads(serialized_inputs, cls=RuntimeDecoder)

# res = runtime_pVQD.main(backend, user_messenger, **deserialized_inputs)
# print(res)
# print("\n")


# def interim_result_callback(job_id, interim_result):
#     print(f"interim result: {interim_result}")


# print("on hardware:")
# backend = provider.get_backend('ibmq_qasm_simulator')
backend = provider.get_backend('ibmq_manila')
# backend = provider.get_backend('ibm_lagos')
# Configure backend options
options = {'backend_name': backend.name()}

# Execute the circuit using the "circuit-runner" program.
job = provider.runtime.run(program_id="p-vqd-xL289veY54",
                           options=options,
                           inputs=inputs,
                           callback=interim_result_callback)

# Get runtime job result.
# begin = time.time()
result = job.result()
# end = time.time()-begin

# result['exec_time'] = end
print(result)
filename = ('data/VQD/runtime_'+backend.name()+'_NN'+str(NN)+'_iter'+str(maxiter)+'_shots'+str(shots)+'.dat')
json.dump(result, open(filename, 'w+'))
