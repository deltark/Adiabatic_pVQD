import numpy as np
import json

from qiskit import IBMQ, QuantumCircuit
from qiskit.providers.ibmq import RunnerResult

from qiskit import Aer
from qiskit.providers.ibmq.runtime import UserMessenger
from qiskit.providers.ibmq.runtime.utils import RuntimeEncoder, RuntimeDecoder

from qiskit import transpile


from qiskit.circuit import ParameterVector
from qiskit.aqua import QuantumInstance
from qiskit.aqua.operators import CircuitSampler, StateFn
from qiskit.aqua.operators.expectations import PauliExpectation

from qiskit.quantum_info import Pauli
from qiskit.aqua.operators import PauliOp, SummedOp

# from pauli_function import *


def generate_pauli(idx_x, idx_z, n):
	'''
	Args:
		n (integer)
		idx (list)
	Returns:
		tensor product of Pauli operators acting on qubits in idx
	'''

	xmask = [0]*n
	zmask = [0]*n
	for i in idx_x:
		xmask[i] = 1
	for i in idx_z:
		zmask[i] = 1

	a_x = np.asarray(xmask, dtype=np.bool)
	a_z = np.asarray(zmask, dtype=np.bool)

	return Pauli(a_z, a_x)


def generate_ising(n_spins, coup, field):
	'''
	Args:
		n_spins (integer)
		coup    (float)
		field   (float)

	Returns:
		Hamiltonian of Ising model with ZZ interaction a X transverse field
	'''

	int_list = []
	field_list = []

	for i in range(n_spins-1):
		int_list.append(generate_pauli([], [i, i+1], n_spins))

	for i in range(n_spins):
		field_list.append(generate_pauli([i], [], n_spins))

	int_coeff = [coup]*len(int_list)
	field_coeff = [field]*len(field_list)

	H = PauliOp(int_list[0], int_coeff[0])

	for i in range(1, len(int_list)):
		H = H + PauliOp(int_list[i], int_coeff[i])

	for i in range(len(field_list)):
		H = H + PauliOp(field_list[i], field_coeff[i])

	return H


def ei(i, n):
	vi = np.zeros(n)
	vi[i] = 1.0
	return vi[:]

def vqe_ansatz(p):
    n_spins = 3

    circuit = QuantumCircuit(n_spins)
    for i in range(n_spins):
	    circuit.h(i)

    count = 0
    for i in range(n_spins):
        circuit.rx(p[count], i)
        count += 1

    for i in range(n_spins-1):
        circuit.cnot(i, i+1)

    for i in range(n_spins):
        circuit.rx(p[count], i)
        count += 1

    return circuit
    

def main(backend, user_messenger, **kwargs):
    """Main entry point of the program.

    Args:
        backend: Backend to submit the circuits to.
        user_messenger: Used to communicate with the program consumer.
        kwargs: User inputs.
    """
    iterations = kwargs.pop('maxiter', 10)
    n_spins = 3
    shots = 100
    nparameters = 2*n_spins

    instance = QuantumInstance(backend=backend, shots=shots)
    sampler = CircuitSampler(instance)

    H = generate_ising(n_spins, -1, -1)

    empty_params = ParameterVector('ep', nparameters)
    ansatz = vqe_ansatz(empty_params)
    var_params = np.zeros(nparameters)
    
    expectation = PauliExpectation()

    H_prj = StateFn(H, is_measurement=True)
    state_wfn = H_prj @ StateFn(ansatz)
    state_wfn = expectation.convert(state_wfn)

    expectH = []

    for it in range(iterations):
        
        # build dictionary of parameters to values
        # {left[0]: parameters[0], .. ., right[0]: parameters[0] + shift[0], ...}

        # First create the dictionary for expectation value
        values_dict = [dict(zip(empty_params, var_params.tolist()))]

        # Then the values for the gradient
        for i in range(nparameters):
            values_dict.append(dict(zip(empty_params, (var_params + ei(i, nparameters)*np.pi/2.0).tolist())))
            values_dict.append(dict(zip(empty_params, (var_params - ei(i, nparameters)*np.pi/2.0).tolist())))

        results = []

        for values in values_dict:
            sampled_op = sampler.convert(state_wfn, params=values)

            mean = sampled_op.eval()[0].real
            est_err = 0

            if (not instance.is_statevector):
                variance = expectation.compute_variance(sampled_op)[0].real
                est_err = np.sqrt(variance/shots)

            results.append([mean, est_err])

        E = [0, 0]
        g = np.zeros((nparameters, 2))

        E[0], E[1] = results[0]
        expectH.append(E)

        for i in range(nparameters):
            rplus = results[1+2*i]
            rminus = results[2+2*i]
            g[i, :] = (rplus[0]-rminus[0])/2.0, np.sqrt(rplus[1]**2+rminus[1]**2)/2.0

        var_params -= 1e-3*g[:, 0]
    
    output = {"E": expectH[:][0], "E_err": expectH[:][1]}

    return output
