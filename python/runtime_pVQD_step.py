from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
# from logging import log
import numpy as np
import time


from qiskit import Aer, execute
from qiskit.quantum_info import Pauli
from qiskit.utils import QuantumInstance
from qiskit.opflow import PauliOp, SummedOp, CircuitSampler, StateFn
from qiskit.circuit import ParameterVector
from qiskit.opflow.evolutions import PauliTrotterEvolution

from qiskit.opflow.expectations import PauliExpectation, expectation_base

from qiskit.providers.ibmq.runtime import UserMessenger

from qiskit.quantum_info import Pauli
from qiskit.opflow import PauliOp, SummedOp


def ei(i, n):
	vi = np.zeros(n)
	vi[i] = 1.0
	return vi[:]


def compute_gradient(state_wfn, parameters, shift, expectator, sampler, right, left, shots):

	nparameters = len(parameters)

	values_dict = []
	# Then the values for the gradient
	for i in range(nparameters):
		values_dict.append(dict(zip(right[:] + left[:], parameters.tolist() + (
			parameters + shift + ei(i, nparameters)*np.pi/2.0).tolist())))
		values_dict.append(dict(zip(right[:] + left[:], parameters.tolist() + (
			parameters + shift - ei(i, nparameters)*np.pi/2.0).tolist())))

	# Now evaluate the circuits with the parameters assigned

	results = []

	for values in values_dict:
		sampled_op = sampler.convert(state_wfn, params=values)

		mean = sampled_op.eval().real
		est_err = 0
		
		variance = expectator.compute_variance(sampled_op).real
		est_err = np.sqrt(variance/shots)

		results.append([mean, est_err])

	g = np.zeros((nparameters, 2))

	for i in range(nparameters):
		rplus = results[2*i]  # 0+2i
		rminus = results[1+2*i]  # 1+2i
		g[i, :] = (rplus[0]-rminus[0])/2.0, np.sqrt(rplus[1]**2+rminus[1]**2)/2.0

	# self.gradient = g

	return g

def compute_overlap(state_wfn, parameters, shift, sampler, right, left):

	# build dictionary of parameters to values
	# {left[0]: parameters[0], .. ., right[0]: parameters[0] + shift[0], ...}

	# First create the dictionary for overlap
	values_dict = [dict(zip(right[:] + left[:],
						parameters.tolist() + (parameters + shift).tolist()))]

	# Now evaluate the circuits with the parameters assigned

	results = []

	for values in values_dict:
		sampled_op = sampler.convert(state_wfn, params=values)

		mean = sampled_op.eval().real
		est_err = 0

		results.append([mean, est_err])

	E = np.zeros(2)

	E[0], E[1] = results[0]

	# self.overlap = E

	return E


def measure_aux_ops(obs_wfn, obs_params, pauli, parameters, expectator, sampler, instance, shots):

	# This function calculates the expectation value of a given operator

	# Prepare the operator and the parameters
	wfn = StateFn(obs_wfn)
	op = StateFn(pauli, is_measurement=True)
	values_obs = dict(zip(obs_params[:], parameters.tolist()))

	braket = op @ wfn

	grouped = expectator.convert(braket)
	sampled_op = sampler.convert(grouped, params=values_obs)

	mean_value = sampled_op.eval().real
	est_err = 0

	if (not instance.is_statevector):
		variance = expectator.compute_variance(sampled_op).real
		est_err = np.sqrt(variance/shots)

	res = [mean_value, est_err]

	return res


# def run(ths, state_wfn, ansatz, parameters, shift, params_vec, obs_params, hamiltonian, instance, user_messenger: UserMessenger, obs_dict={}, max_iter=100, shots=2000, **kwargs):
def main(backend, user_messenger, **kwargs):

	ths         = kwargs.pop('threshold', 0.99999)
	state_wfn   = kwargs.pop('state_wfn') #circuit
	obs_wfn     = kwargs.pop('obs_wfn')
	# ansatz      = kwargs.pop('ansatz')
	parameters  = kwargs.pop('parameters')
	shift       = kwargs.pop('shift')
	# params_vec  = kwargs.pop('params_vec')
	obs_params  = list(obs_wfn.parameters)
	# hamiltonian = kwargs.pop('hamiltonian')
	obs_dict    = kwargs.pop('obs_dict', {})
	max_iter    = kwargs.pop('max_iter', 50)
	shots       = kwargs.pop('shots', 2000)

	r_circ      = kwargs.pop('r_circ')
	l_circ      = kwargs.pop('l_circ')

	right       = list(r_circ.parameters)
	left        = list(l_circ.parameters)

	# state_wfn.assign_parameters(right + left)

	first_measure = kwargs.pop('first_measure', False)


	instance = QuantumInstance(backend=backend, shots=shots)

	expectation = PauliExpectation()
	sampler = CircuitSampler(instance)

	## Now prepare the state in order to compute the overlap and its gradient
	state_wfn = expectation.convert(state_wfn)
	# right_left = state_wfn._parameters
	# print(right_left)

	## Also the standard state for measuring the observables
	# obs_wfn = ansatz.assign_parameters({params_vec: obs_params})

	#######
	#START
	result = {}

	count = 0
	overlap = [0.01, 0]
	# g_norm = 1

	# and g_norm > len(params)*8e-6:
	if not first_measure:
		while overlap[0] < ths and count < max_iter:

			print("Shift optimizing step:", count+1)
			count = count + 1

			# if opt == 'momentum':
			# 	old_grad = np.asarray(g[:, 0])
			# ## Measure energy and gradient

			overlap_backend = Aer.get_backend('statevector_simulator')
			overlap_instance = QuantumInstance(backend=overlap_backend, shots=1)
			overlap_sampler = CircuitSampler(overlap_instance)
			E = compute_overlap(
				state_wfn, parameters, shift, overlap_sampler, right, left)
			g = compute_gradient(
				state_wfn, parameters, shift, expectation, sampler, right, left, shots)

			overlap = E
			gradient = g

			# tot_steps = tot_steps+1

			# if count == 1:
			# 	initial_fidelities.append(self.overlap[0])
			# 	err_init_fid.append(self.overlap[1])

				# shifts_t.append(list(self.shift))

			print('Overlap', overlap)
			print('Gradient', gradient[:, 0])

			#SGD
			shift = shift + g[:, 0]

			# shifts_t.append(list(self.shift))

			#Norm of the gradient
			# g_vec = np.asarray(g[:, 0])
			# g_norm = np.linalg.norm(g_vec)

			# grad_t.append(list(self.gradient[:, 0]))
			# err_grad_t.append(list(self.gradient[:, 1]))
			# grad_norms_t.append(g_norm)

			# user_messenger.publish(
			# 	{"grad_t": grad_t, "err_grad_t": err_grad_t, "t_step": i+1})

		# Update parameters
		parameters = parameters + shift

		result["new_params"] = parameters
		result["shift"] = shift
		result["overlap"] = overlap
		result["gradient"] = gradient
		result["count"] = count


		print('\n---------------------------------- \n')

		print("Shift after optimizing:",shift)
		print("New parameters:"        ,parameters)

		print("New overlap: "          ,overlap[0])


	#Measure observables
	if len(obs_dict) > 0:
		
		obs_measure = {}
		obs_error = {}

		for (obs_name, obs_pauli) in obs_dict.items():
			run_measure = measure_aux_ops(
				obs_wfn, obs_params, obs_pauli, parameters, expectation, sampler, instance, shots)
			obs_measure[str(obs_name)] = run_measure[0]
			obs_error['err_'+str(obs_name)] = run_measure[1]

	#save and return results

	result["obs_measure"] = obs_measure
	result["obs_error"]   = obs_error

	return result
