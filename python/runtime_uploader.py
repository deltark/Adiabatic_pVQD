import os
from qiskit import IBMQ

IBMQ.load_account()
# Substitute with your provider.
provider = IBMQ.get_provider(
    hub='ibm-q-research-2', group='epfl-2', project='main')

# provider.runtime.delete_program('my-vqe-2')

sample_program_data = os.path.join(
    os.getcwd(), "runtime_pVQD.py")
sample_program_json = os.path.join(
    os.getcwd(), "runtime_pVQD.json")

# This will fail if a sample-program already exists.
program_id = provider.runtime.upload_program(
    data=sample_program_data,
    metadata=sample_program_json
)
print(program_id)
