import os
from qiskit import IBMQ

IBMQ.load_account()
# Substitute with your provider.
provider = IBMQ.get_provider(
    hub='ibm-q-research-2', group='epfl-2', project='main')

# provider.runtime.delete_program('my-vqe-2')

program_data = os.path.join(
    os.getcwd(), "runtime_pVQD.py")
program_json = os.path.join(
    os.getcwd(), "runtime_pVQD.json")

program_id = "p-vqd-xL289veY54"

provider.runtime.update_program(program_id=program_id,
    data=program_data,
    metadata=program_json
)

# program_id = provider.runtime.upload_program(
#     data=program_data,
#     metadata=program_json
# )
