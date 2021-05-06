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
from hamiltonian_op import IsingHamiltonian


# This class aims to simulate the dynamics of a quantum system
# approximating it with a variational ansatz whose parameters
# are varied in order to follow the unitary evolution

# Useful functions

def projector_zero(n_qubits):
	from qiskit.aqua.operators            import Z, I

	prj = (1/np.power(2,n_qubits))*(I+Z)

	for a in range(n_qubits-1):
		prj = prj^(I+Z)

	return prj

def ei(i,n):
	vi = np.zeros(n)
	vi[i] = 1.0
	return vi[:]



class pVQD:

	def __init__(self,hamiltonian,ansatz,parameters,initial_shift,instance,shots,ham_tfunc=None):

		'''
		Args:

		hamiltonian   : [list of operators or operator] list of Hamiltonian parts to be summed (e.g. [Hx, Hzz]),
			or just a single Hamiltonian, time-dependent parts first
		parameters    : [numpy.array] an array containing the parameters of the ansatz
		initial_shift : [numpy.array] an array containing the initial guess of shifts
		ham_tfunc     : [list of lambda functions or lambda function] list of time-dependent functions to be
			multiplied to the Hamiltonian parts in the same order

		'''
		self.hamiltonian     = hamiltonian
		self.instance        = instance
		self.parameters      = parameters
		self.num_parameters  = len(parameters)
		self.shift           = initial_shift
		self.shots           = shots
		self.ham_tfunc       = ham_tfunc

		## Initialize quantities that will be equal all over the calculation
		self.params_vec      = ParameterVector('p',self.num_parameters)
		self.ansatz          = ansatz(self.params_vec)


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



	def construct_total_circuit(self,time_step):
		## This function creates the circuit that will be used to evaluate overlap and its gradient

		# First, create the Trotter step

		# Make the Hamiltonian a list if it's a single one
		if not isinstance(self.hamiltonian, list):
			self.hamiltonian = [self.hamiltonian]

		if self.ham_tfunc is not None:
			step_h  = self.ham_params*np.array(self.hamiltonian[:len(self.ham_tfunc)])*time_step
			step_h  = np.append(step_h, np.array(self.hamiltonian[len(self.ham_tfunc):])*time_step)
		else:
			step_h  = time_step*np.array(self.hamiltonian)

		trotter = PauliTrotterEvolution(reps=1)
		# Total Trotter circuit constructed by summing over the Hamiltonian parts
		U_dt    = np.sum([trotter.convert(step_h[j].exp_i()).to_circuit() for j in range(len(step_h))])


		l_circ  = self.ansatz.assign_parameters({self.params_vec: self.left})
		r_circ  = self.ansatz.assign_parameters({self.params_vec: self.right})

		## Projector
		zero_prj = StateFn(projector_zero(self.hamiltonian[0].num_qubits),is_measurement = True)
		state_wfn = zero_prj @ StateFn(r_circ +U_dt+ l_circ.inverse())


		return state_wfn


	# This function calculate overlap and gradient of the overlap using a global operator on the |0><0| state

	def compute_overlap_and_gradient(self,state_wfn,parameters,shift,expectator,sampler):

		nparameters = len(parameters)
		# build dictionary of parameters to values
		# {left[0]: parameters[0], .. ., right[0]: parameters[0] + shift[0], ...}

		# First create the dictionary for overlap
		values_dict = [dict(zip(self.right[:] + self.left[:], parameters.tolist() + (parameters + shift).tolist()))]


		# Then the values for the gradient
		for i in range(nparameters):
			values_dict.append(dict(zip(self.right[:] + self.left[:], parameters.tolist() + (parameters + shift + ei(i,nparameters)*np.pi/2.0).tolist())))
			values_dict.append(dict(zip(self.right[:] + self.left[:], parameters.tolist() + (parameters + shift - ei(i,nparameters)*np.pi/2.0).tolist())))

		# Now evaluate the circuits with the parameters assigned

		results = []

		for values in values_dict:
			sampled_op = sampler.convert(state_wfn,params=values)

			mean  = sampled_op.eval()[0].real
			#mean  = np.power(np.absolute(mean),2)
			est_err = 0


			if (not self.instance.is_statevector):
				variance = expectator.compute_variance(sampled_op)[0].real
				est_err  = np.sqrt(variance/self.shots)

			results.append([mean,est_err])

		E = np.zeros(2)
		g = np.zeros((nparameters,2))

		E[0],E[1] = results[0]

		for i in range(nparameters):
			rplus  = results[1+2*i]
			rminus = results[2+2*i]
			# G      = (Ep - Em)/2
			# var(G) = var(Ep) * (dG/dEp)**2 + var(Em) * (dG/dEm)**2
			g[i,:] = (rplus[0]-rminus[0])/2.0,np.sqrt(rplus[1]**2+rminus[1]**2)/2.0

		self.overlap  = E
		self.gradient = g

		return E,g

	## this function does the same thing but uses SPSA

	def compute_overlap_and_gradient_spsa(self,state_wfn,parameters,shift,expectator,sampler,count):

		nparameters = len(parameters)
		# build dictionary of parameters to values
		# {left[0]: parameters[0], .. ., right[0]: parameters[0] + shift[0], ...}

		# Define hyperparameters
		c  = 0.1
		a  = 0.16
		A  = 1
		alpha  = 0.602
		gamma  = 0.101

		a_k = a/np.power(A+count,alpha)
		c_k = c/np.power(count,gamma)

		# Determine the random shift

		delta = np.random.binomial(1,0.5,size=nparameters)
		delta = np.where(delta==0, -1, delta)
		delta = c_k*delta

		# First create the dictionary for overlap
		values_dict = [dict(zip(self.right[:] + self.left[:], parameters.tolist() + (parameters + shift).tolist()))]


		# Then the values for the gradient

		values_dict.append(dict(zip(self.right[:] + self.left[:], parameters.tolist() + (parameters + shift + delta).tolist())))
		values_dict.append(dict(zip(self.right[:] + self.left[:], parameters.tolist() + (parameters + shift - delta).tolist())))

		# Now evaluate the circuits with the parameters assigned

		results = []

		for values in values_dict:
			sampled_op = sampler.convert(state_wfn,params=values)

			mean  = sampled_op.eval()[0]
			mean  = np.power(np.absolute(mean),2)
			est_err = 0


			if (not self.instance.is_statevector):
				variance = expectator.compute_variance(sampled_op)[0].real
				est_err  = np.sqrt(variance/self.shots)

			results.append([mean,est_err])

		E = np.zeros(2)
		g = np.zeros((nparameters,2))

		E[0],E[1] = results[0]

		# and the gradient

		rplus  = results[1]
		rminus = results[2]

		for i in range(nparameters):
			# G      = (Ep - Em)/2Δ_i
			# var(G) = var(Ep) * (dG/dEp)**2 + var(Em) * (dG/dEm)**2
			g[i,:] = a_k*(rplus[0]-rminus[0])/(2.0*delta[i]),np.sqrt(rplus[1]**2+rminus[1]**2)/(2.0*delta[i])

		self.overlap  = E
		self.gradient = g

		return E,g


	def measure_aux_ops(self,obs_wfn,pauli,parameters,expectator,sampler):

		# This function calculates the expectation value of a given operator

		# Prepare the operator and the parameters
		wfn = StateFn(obs_wfn)
		op  = StateFn(pauli,is_measurement = True)
		values_obs = dict(zip(self.obs_params[:], parameters.tolist()))

		braket = op @ wfn

		grouped    = expectator.convert(braket)
		sampled_op = sampler.convert(grouped,params = values_obs)

		mean_value = sampled_op.eval()[0].real
		est_err = 0

		if (not self.instance.is_statevector):
			variance = expectator.compute_variance(sampled_op)[0].real
			est_err  = np.sqrt(variance/self.shots)

		res = [mean_value,est_err]

		return res

	def adam_gradient(self,count,m,v,g):
		## This function implements adam optimizer
		beta1 = 0.9
		beta2 = 0.999
		eps   = 1e-8
		alpha = [0.001 for i in range(len(self.parameters))]
		if count == 0:
			count = 1

		new_shift = [0 for i in range(len(self.parameters))]

		for i in range(len(self.parameters)):
			m[i] = beta1 * m[i] + (1 - beta1) * g[i]
			v[i] = beta2 * v[i] + (1 - beta2) * np.power(g[i],2)

			alpha[i] = alpha[i] * np.sqrt(1 - np.power(beta2,count)) / (1 - np.power(beta1,count))

			new_shift[i] = self.shift[i] + alpha[i]*(m[i]/(np.sqrt(v[i])+eps))

		return new_shift


	def run(self,ths,timestep,n_steps,obs_dict={},filename='algo_result.dat',max_iter = 100,opt='sgd',grad='param_shift',initial_point=None):



		## initialize useful quantities once
		if(self.instance.is_statevector):
			expectation = PauliExpectation()
		if(not self.instance.is_statevector):
			expectation = PauliExpectation()

		sampler = CircuitSampler(self.instance)

		## Now prepare the state in order to compute the overlap and its gradient
		state_wfn = self.construct_total_circuit(timestep)
		state_wfn = expectation.convert(state_wfn)


		## Also the standard state for measuring the observables
		obs_wfn   = self.ansatz.assign_parameters({self.params_vec: self.obs_params})

		#######################################################

		times = [i*timestep for i in range(n_steps+1)]

		if self.ham_tfunc is not None:
			## Prepare the time-dependent Hamiltonian parameters
			ham_tfunc_values = np.array([[self.ham_tfunc[i](t) for t in times] for i in range(len(self.ham_tfunc))])
			ham_dict = [dict(zip(self.ham_params[:], ham_tfunc_values[:,i].tolist())) for i in range(n_steps+1)]

		tot_steps= 0

		if initial_point != None :
			if len(initial_point) != len(self.parameters):
				print("TypeError: Initial parameters are not of the same size of circuit parameters")
				return


			print("\nRestart from: ")
			print(initial_point)
			self.parameters = initial_point

		print("Running the algorithm")

		#prepare observables for quench

		if len(obs_dict) > 0:
			obs_measure = {}
			obs_error   = {}

			for (obs_name,obs_pauli) in obs_dict.items():
				first_measure                   = self.measure_aux_ops(obs_wfn,obs_pauli,self.parameters,expectation,sampler)
				obs_measure[str(obs_name)]      = [first_measure[0]]
				obs_error['err_'+str(obs_name)] = [first_measure[1]]


		counter = []
		initial_fidelities = []
		fidelities = []
		err_fin_fid = []
		err_init_fid = []
		params = []

		params.append(list(self.parameters))


		for i in range(n_steps):

			print('\n================================== \n')
			print("Time slice:",i+1)
			print("Shift before optimizing this step:",self.shift)
			print("Initial parameters:", self.parameters)
			print('\n================================== \n')

			if self.ham_tfunc is not None:
				#update time dependent Hamiltonian
				state_wfn_Ht = state_wfn.assign_parameters(ham_dict[i+1])
				print(state_wfn_Ht)
			else:
				state_wfn_Ht = state_wfn

			count = 0
			self.overlap = [0.01,0]
			g_norm = 1

			if opt == 'adam':
				m = np.zeros(len(self.parameters))
				v = np.zeros(len(self.parameters))

			if opt == 'momentum':
				old_grad = np.zeros(len(self.parameters))
				g        = np.zeros((len(self.parameters),2))

			# Save data of overlaps for analysis
			overlap_history = []


			while self.overlap[0] < ths and count < max_iter: # and g_norm > len(params)*8e-6:
				print("Shift optimizing step:",count+1)
				count = count +1

				if opt == 'momentum':
					old_grad = np.asarray(g[:,0])
				## Measure energy and gradient

				if grad == 'param_shift':
					E,g = self.compute_overlap_and_gradient(state_wfn_Ht,self.parameters,self.shift,expectation,sampler)
				if grad == 'spsa':
					E,g = self.compute_overlap_and_gradient_spsa(state_wfn_Ht,self.parameters,self.shift,expectation,sampler,count)

				tot_steps = tot_steps+1

				if count == 1:
					initial_fidelities.append(self.overlap[0])
					err_init_fid.append(self.overlap[1])

				print('Overlap',self.overlap)
				print('Gradient',self.gradient[:,0])

				overlap_history.append(self.overlap[0])

				if opt == 'adam':
					print("\n Adam \n")
					meas_grad = np.asarray(g[:,0])
					self.shift = np.asarray(self.adam_gradient(count,m,v,meas_grad))

				if opt == 'momentum':
					print("Momentum")
					m_grad = np.asarray(g[:,0]) + 0.9*old_grad
					self.shift = self.shift + m_grad

				elif opt== 'sgd':
					self.shift = self.shift + g[:,0]

				#Norm of the gradient
				g_vec = np.asarray(g[:,0])
				g_norm = np.linalg.norm(g_vec)


			# Update parameters


			print('\n---------------------------------- \n')

			print("Shift after optimizing:",self.shift)
			print("New parameters:"        ,self.parameters + self.shift)

			print("New overlap: "          ,self.overlap[0])

			self.parameters = self.parameters + self.shift



			# Measure quantities and save them

			if len(obs_dict) > 0:

				for (obs_name,obs_pauli) in obs_dict.items():
					run_measure   = self.measure_aux_ops(obs_wfn,obs_pauli,self.parameters,expectation,sampler)
					obs_measure[str(obs_name)].append(run_measure[0])
					obs_error['err_'+str(obs_name)].append(run_measure[1])

			counter.append(count)
			fidelities.append(self.overlap[0])
			err_fin_fid.append(self.overlap[1])


			params.append(list(self.parameters))




		print("Total measurements:",tot_steps)
		print("Measure per step:", tot_steps/n_steps)

		print("overlap_history ", overlap_history)

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

		log_data['overlap_history'] = overlap_history

		json.dump(log_data, open( filename,'w+'))
