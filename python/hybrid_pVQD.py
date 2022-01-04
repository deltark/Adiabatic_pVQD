from logging import log
import numpy as np
import json
import matplotlib.pyplot as plt
from scipy   import  linalg as LA



from qiskit                               import Aer, execute
from qiskit.quantum_info 			      import Pauli
from qiskit.utils                         import QuantumInstance
from qiskit.opflow 			              import PauliOp, SummedOp, CircuitSampler, StateFn
from qiskit.circuit                       import ParameterVector
from qiskit.opflow.evolutions             import Trotter, PauliTrotterEvolution

from qiskit.opflow.state_fns      import CircuitStateFn
from qiskit.opflow.expectations   import PauliExpectation, AerPauliExpectation, MatrixExpectation
from qiskit.opflow.primitive_ops  import CircuitOp
from qiskit.opflow                import Z, I

from qiskit.providers.ibmq.runtime       import UserMessenger
from qiskit.providers.ibmq.runtime.utils import RuntimeEncoder, RuntimeDecoder

from pauli_function import *

import runtime_pVQD_step


# This class aims to simulate the dynamics of a quantum system
# approximating it with a variational ansatz whose parameters
# are varied in order to follow the unitary evolution

# Useful functions

def projector_zero(n_qubits):
	# This function create the global projector |00...0><00...0|
	from qiskit.opflow import Z, I

	prj_list = [0.5*(I+Z) for i in range(n_qubits)]
	prj = prj_list[0]

	for a in range(1, len(prj_list)):
		prj = prj ^ prj_list[a]

	return prj



class pVQD:

	def __init__(self,hamiltonian,ansatz,ansatz_depth,parameters,initial_shift,instance,shots,ham_tfunc=None,ham_integ=None,ham_mag2=None):

		'''
		Args:

		hamiltonian   : [list of operators or operator] list of Hamiltonian parts to be summed (e.g. [Hx, Hzz]),
			or just a single Hamiltonian, time-dependent parts first
		parameters    : [numpy.array] an array containing the parameters of the ansatz
		initial_shift : [numpy.array] an array containing the initial guess of shifts
		ham_tfunc     : [list of lambda functions or lambda function] list of time-dependent functions to be
			multiplied to the Hamiltonian parts in the same order
		ham_integ     : [list of lambda functions or lambda function] list of indefinite integrals of the
			time-dependent Hamiltonian for Magnus expansion
		ham_mag2      : [list of operators or operator] list of terms of the second-order Magnus expansion
		    of the time-dependent Hamiltonian

		'''
		self.hamiltonian     = hamiltonian
		self.ansatz_depth    = ansatz_depth
		self.instance        = instance
		self.parameters      = parameters
		self.num_parameters  = len(parameters)
		self.shift           = initial_shift
		self.shots           = shots
		self.ham_tfunc       = ham_tfunc
		self.ham_integ       = ham_integ
		self.ham_mag2        = ham_mag2

		## Initialize quantities that will be equal all over the calculation
		self.params_vec      = ParameterVector('p',self.num_parameters)
		self.ansatz          = ansatz(self.hamiltonian[0].num_qubits, self.ansatz_depth, self.params_vec)


		# ParameterVector for left and right circuit

		self.left  = ParameterVector('l', self.ansatz.num_parameters)
		self.right = ParameterVector('r', self.ansatz.num_parameters)

		# ParameterVector for measuring abservables
		self.obs_params = ParameterVector('θ',self.ansatz.num_parameters)

		# ParameterVector for time-evolving Trotter circuit (Hamiltonian)
		if self.ham_tfunc is not None:
			if isinstance(self.ham_tfunc, list):
				self.ham_params = ParameterVector('h', len(self.ham_tfunc))
			else:
				self.ham_params = ParameterVector('h', 1)

		self.njobs = 0 #counter for number of jobs submitted to hardware

		# ParameterVector for Magnus expansion
		# if self.ham_integ is not None:
		# 	if isinstance(self.ham_integ, list):
		# 		self.magnus_params = ParameterVector('m', len(self.ham_integ))
		# 	else:
		# 		self.magnus_params = ParameterVector('m', 1)

	def construct_total_circuit(self, time_step):
		## This function creates the circuit that will be used to evaluate overlap and its gradient

		# First, create the Trotter step
		if self.ham_tfunc is not None:

			#former method
			step_h = self.ham_params * \
				np.array(self.hamiltonian[:len(self.ham_tfunc)])*time_step
			step_h = np.append(step_h, np.array(
				self.hamiltonian[len(self.ham_tfunc):])*time_step)

			#Magnus expansion
			# step_h  = self.magnus_params*np.array(self.hamiltonian)
			# step_h  = np.append(step_h, time_step*np.array(self.ham_mag2))

		else:
			step_h = time_step*np.array(self.hamiltonian)
			print(step_h)

		trotter = PauliTrotterEvolution(reps=1)
		# Total Trotter circuit constructed by summing over the Hamiltonian parts
		if len(step_h) > 1:
			U_dt = trotter.convert(step_h[0].exp_i()).to_circuit()
			for j in range(1, len(step_h)):
				U_dt += trotter.convert(step_h[j].exp_i()).to_circuit()
		else:
			U_dt = trotter.convert(step_h[0][0].exp_i()).to_circuit()

		l_circ = self.ansatz.assign_parameters({self.params_vec: self.left})
		r_circ = self.ansatz.assign_parameters({self.params_vec: self.right})

		## Projector
		zero_prj = StateFn(projector_zero(
			self.hamiltonian[0].num_qubits), is_measurement=True)
		state_wfn = zero_prj @ StateFn(r_circ + U_dt + l_circ.inverse())

		return state_wfn

	def construct_hamiltonian(self):
		## This function creates the circuit that will be used to measure the time-dependent energy observable
		ham_circuit = np.sum(self.ham_params*np.array(self.hamiltonian[:len(self.ham_tfunc)])) \
				+ np.sum(np.array(self.hamiltonian[len(self.ham_tfunc):]))

		return ham_circuit



	def run(self,ths,timestep,n_steps,obs_dict={},filename='algo_result.dat',max_iter = 100,initial_point=None):


		## Now prepare the state in order to compute the overlap and its gradient
		state_wfn = self.construct_total_circuit(timestep)

		## Also the standard state for measuring the observables
		obs_wfn   = self.ansatz.assign_parameters({self.params_vec: self.obs_params})

		#######################################################

		times = [i*timestep for i in range(n_steps+1)]

		if self.ham_tfunc is not None:
			## Prepare the time-dependent Hamiltonian parameters
			ham_tfunc_values = np.array([[self.ham_tfunc[i](t) for t in times] for i in range(len(self.ham_tfunc))])
			ham_dict = [dict(zip(self.ham_params[:], ham_tfunc_values[:,i].tolist())) for i in range(n_steps+1)]

			## And the values of the integral parameters for Magnus
			# ham_integ_values = np.array([[self.ham_integ[i](times[j+1]) - self.ham_integ[i](times[j])
			# 				   for j in range(n_steps)] for i in range(len(self.ham_integ))])
			# magnus_dict = [dict(zip(self.magnus_params[:], ham_integ_values[:,i].tolist())) for i in range(n_steps)]

			## And the H(t) operator for observable measurement
			ham_circuit = self.construct_hamiltonian()
			obs_dict['E(t)'] = ham_circuit.assign_parameters(ham_dict[0])

		tot_steps= 0

		if initial_point != None :
			if len(initial_point) != len(self.parameters):
				print("TypeError: Initial parameters are not of the same size of circuit parameters")
				return


			print("\nRestart from: ")
			print(initial_point)
			self.parameters = initial_point

		print("Running the algorithm")

		#It's runtime
		inputs = {}
		inputs["threshold"] = ths
		inputs["state_wfn"] = state_wfn
		inputs["obs_wfn"] = obs_wfn
		inputs["parameters"] = self.parameters
		inputs["shift"] = self.shift
		inputs['obs_params'] = self.obs_params
		inputs['obs_dict'] = obs_dict
		inputs['max_iter'] = max_iter
		inputs['shots'] = self.shots

		inputs['first_measure'] = True

		backend = Aer.get_backend('qasm_simulator')
		user_messenger = UserMessenger()
		serialized_inputs = json.dumps(inputs, cls=RuntimeEncoder)
		deserialized_inputs = json.loads(serialized_inputs, cls=RuntimeDecoder)
		runtime_res = runtime_pVQD_step.main(
                    backend, user_messenger, **deserialized_inputs)

		inputs['first_measure'] = False

		#prepare observables for quench

		if len(obs_dict) > 0:
			obs_measure = {}
			obs_error   = {}

			for (obs_name,obs_pauli) in obs_dict.items():
				obs_measure[str(obs_name)].append(runtime_res["obs_measure"][str(obs_name)])
				obs_error['err_'+str(obs_name)].append(runtime_res["obs_measure"]['err_'+str(obs_name)])


		counter = []
		initial_fidelities = []
		fidelities = []
		err_fin_fid = []
		err_init_fid = []
		params = []

		interm_fidelities = [] #save intermediate fidelities
		shifts = [] #save all parameter shifts
		gradients = [] #save all evaluated gradients
		err_grad = []
		grad_norms = []
		njobs = [] #save the number of jobs submitted to hardware

		params.append(list(self.parameters))


		for i in range(n_steps):

			print('\n================================== \n')
			print("Time slice:",i+1)
			print("Shift before optimizing this step:",self.shift)
			print("Initial parameters:", self.parameters)
			print('\n================================== \n')

			interm_fid_t = [] #will go inside interm_fidelities
			shifts_t = [] #will go inside shifts
			grad_t = []
			err_grad_t = []
			grad_norms_t = []

			self.njobs = 0 #reset njobs

			if self.ham_tfunc is not None:
				#update time dependent Hamiltonian
				state_wfn_Ht = state_wfn.assign_parameters(ham_dict[i+1])
				print(state_wfn_Ht)

				#magnus expansion
				# state_wfn_Ht = state_wfn.assign_parameters(magnus_dict[i])

				#update observable
				obs_dict['E(t)'] = ham_circuit.assign_parameters(ham_dict[i+1])
			else:
				state_wfn_Ht = state_wfn

			self.overlap = [0.01, 0]

			#It's runtime
			inputs["state_wfn"]  = state_wfn_Ht
			inputs["parameters"] = self.parameters
			inputs["shift"]      = self.shift
			inputs['obs_dict']   = obs_dict

			serialized_inputs = json.dumps(inputs, cls=RuntimeEncoder)
			deserialized_inputs = json.loads(serialized_inputs, cls=RuntimeDecoder)
			runtime_res = runtime_pVQD_step.main(backend, user_messenger, **deserialized_inputs)

			self.parameters = runtime_res["new_params"]
			self.shift      = runtime_res["shift"]
			self.overlap    = runtime_res["overlap"]
			self.gradient   = runtime_res["gradient"]
			count           = runtime_res["count"]



			# Measure quantities and save them

			if len(obs_dict) > 0:

				for (obs_name,obs_pauli) in obs_dict.items():
					obs_measure[str(obs_name)].append(runtime_res["obs_measure"][str(obs_name)])
					obs_error['err_'+str(obs_name)].append(runtime_res["obs_measure"]['err_'+str(obs_name)])

			counter.append(count)
			fidelities.append(self.overlap[0])
			err_fin_fid.append(self.overlap[1])

			interm_fidelities.append(list(interm_fid_t))
			shifts.append(list(shifts_t))

			gradients.append(list(grad_t))
			err_grad.append(list(err_grad_t))
			grad_norms.append(list(grad_norms_t))

			njobs.append(self.njobs)

			params.append(list(self.parameters))




		print("Total measurements:",tot_steps)
		print("Measure per step:", tot_steps/n_steps)

		# print("overlap_history ", overlap_history)

		# Save data on file

		log_data = {}
		if len(obs_dict) > 0:
			for (obs_name,obs_pauli) in obs_dict.items():
				log_data[str(obs_name)]        = obs_measure[str(obs_name)]
				log_data['err_'+str(obs_name)] = obs_error['err_'+str(obs_name)]

		log_data['init_F']      = initial_fidelities
		log_data['final_F']     = fidelities
		log_data['err_init_F']  = err_init_fid
		log_data['err_fin_F']   = err_fin_fid
		log_data['iter_number'] = counter
		log_data['times']       = times
		log_data['params']      = list(params)
		log_data['tot_steps']   = [tot_steps]

		log_data['interm_F'] = list(interm_fidelities)
		log_data['shifts'] = list(shifts)
		log_data['gradients'] = list(gradients)
		log_data['err_grad'] = list(err_grad)
		log_data['norm_grad'] = list(grad_norms)
		log_data['njobs'] = list(njobs)

		# log_data['overlap_history'] = overlap_history

		json.dump(log_data, open( filename,'w+'))
