"""A sample runtime program that submits random circuits for user-specified iterations."""

from qiskit import Aer
from qiskit.providers.ibmq.runtime import UserMessenger
from qiskit.providers.ibmq.runtime.utils import RuntimeEncoder, RuntimeDecoder
from qiskit.providers.ibmq import RunnerResult
import random
import json

from qiskit import transpile
from qiskit.circuit.random import random_circuit

from qiskit import IBMQ

import runtime_vqe

IBMQ.load_account()
provider = IBMQ.get_provider(
    hub='ibm-q-research-2', group='epfl-2', project='main')

# runtime_backends = provider.backends(input_allowed='runtime')
# print(f"Backends that support Qiskit Runtime: {runtime_backends}")

inputs = {"iterations": 2}

# backend = Aer.get_backend('qasm_simulator')
# user_messenger = UserMessenger()
# serialized_inputs = json.dumps(inputs, cls=RuntimeEncoder)
# deserialized_inputs = json.loads(serialized_inputs, cls=RuntimeDecoder)

# res = runtime_vqe.main(backend, user_messenger, **deserialized_inputs)
# print(res)
# print("\n")


# def interim_result_callback(job_id, interim_result):
#     print(f"interim result: {interim_result}")


print("on hardware:")
backend = provider.get_backend('ibmq_qasm_simulator')
# Configure backend options
options = {'backend_name': backend.name()}

# Execute the circuit using the "circuit-runner" program.
job = provider.runtime.run(program_id="my-vqe-2",
                           options=options,
                           inputs=inputs)

# Get runtime job result.
result = job.result()
print(result)
