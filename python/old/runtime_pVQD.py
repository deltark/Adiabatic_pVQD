from logging import raiseExceptions
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
# from logging import log
import numpy as np
import time
# import json
# import functools
# import itertools
import matplotlib.pyplot as plt
# from scipy   import  linalg as LA
import copy



# import qiskit
from qiskit                               import Aer, execute
from qiskit.quantum_info 			      import Pauli
from qiskit.utils                         import QuantumInstance
from qiskit.opflow 			      		  import PauliOp, SummedOp, CircuitSampler, StateFn
from qiskit.circuit                       import ParameterVector
from qiskit.opflow.evolutions     		  import PauliTrotterEvolution

# from qiskit.opflow.state_fns      import CircuitStateFn
from qiskit.opflow.expectations   import PauliExpectation
# from qiskit.opflow.primitive_ops  import CircuitOp
# from qiskit.opflow                import Z, I

from qiskit.providers.ibmq.runtime import UserMessenger

### pauli_function.py
from qiskit.quantum_info import Pauli
from qiskit.opflow import PauliOp, SummedOp

#error mitigation
from qiskit.ignis.mitigation.measurement import CompleteMeasFitter


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


def generate_ising_pbc(n_spins, coup, field):
	'''
	Args:
		n_spins (integer)
		coup    (float)
		field   (float)

	Returns:
		Hamiltonian of Ising model with ZZ interaction a X transverse field, pbc
	'''

	int_list = []
	field_list = []

	int_list.append(generate_pauli([], [0, n_spins-1], n_spins))

	if(n_spins > 2):
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


def generate_ising_Hzz(n_spins, coup):

	int_list = []

	for i in range(n_spins-1):
		int_list.append(generate_pauli([], [i, i+1], n_spins))

	int_coeff = [coup]*len(int_list)

	H = PauliOp(int_list[0], int_coeff[0])

	for i in range(1, len(int_list)):
		H = H + PauliOp(int_list[i], int_coeff[i])

	return H


def generate_ising_Hx(n_spins, field):

	field_list = []

	for i in range(n_spins):
		field_list.append(generate_pauli([i], [], n_spins))

	field_coeff = [field]*len(field_list)

	H = PauliOp(field_list[0], field_coeff[0])

	for i in range(1, len(field_list)):
		H = H + PauliOp(field_list[i], field_coeff[i])

	return H


def generate_magnus_2(n_spins, coup, field):

	listYZ = []
	listZY = []

	for i in range(n_spins-1):
		listYZ.append(generate_pauli([i], [i, i+1], n_spins))
		listZY.append(generate_pauli([i+1], [i, i+1], n_spins))

	coeff = [coup*field]*len(listYZ)

	H2 = PauliOp(listYZ[0], coeff[0])
	H2 += PauliOp(listZY[0], coeff[0])

	for i in range(1, len(listYZ)):
		H2 += PauliOp(listYZ[i], coeff[i])
		H2 += PauliOp(listZY[i], coeff[i])

	return H2
###################################################################

# ansatze.py

## These file will contain all the ansatze used for variational quantum simulation


#=========================================

def hweff_ansatz(p):
	n_spins = 3
	count = 0
	circuit = QuantumCircuit(n_spins)
	depth = 2

	for j in range(depth):

		if(j % 2 == 0):
			# Rx - Rzz block
			for i in range(n_spins):
				circuit.rx(p[count], i)
				count = count + 1

			circuit.barrier()

			for i in range(n_spins-1):
				circuit.rzz(p[count], i, i+1)
				count = count + 1

			circuit.barrier()

		if(j % 2 == 1):
			for i in range(n_spins):
				circuit.ry(p[count], i)
				count = count + 1

			circuit.barrier()

			for i in range(n_spins-1):
				circuit.rzz(p[count], i, i+1)
				count = count + 1

			circuit.barrier()

	# Final block to close the ansatz
	if (depth % 2 == 1):
		for i in range(n_spins):
			circuit.ry(p[count], i)
			count = count + 1
	if (depth % 2 == 0):
		for i in range(n_spins):
			circuit.rx(p[count], i)
			count = count + 1

	return circuit

#==========================================


def hweff_ansatz_adiab(n_spins, depth, p):
	# n_spins = 3
	count = 0
	circuit = QuantumCircuit(n_spins)
	# depth = 3

	for i in range(n_spins):
		circuit.h(i)

	for j in range(depth):

		if(j % 2 == 0):
			# Rzz - Rx block
			for i in range(n_spins-1):
				circuit.rzz(p[count], i, i+1)
				count = count + 1

			circuit.barrier()

			for i in range(n_spins):
				circuit.rx(p[count], i)
				count = count + 1

			circuit.barrier()

		if(j % 2 == 1):
			for i in range(n_spins-1):
				circuit.rzz(p[count], i, i+1)
				count = count + 1

			circuit.barrier()

			for i in range(n_spins):
				circuit.ry(p[count], i)
				count = count + 1

	return circuit

def hweff_ansatz_adiab_firststep(n_spins, depth, p):
	# n_spins = 3
	count = 0
	circuit = QuantumCircuit(n_spins)
	# depth = 3

	for i in range(n_spins):
		circuit.h(i)

	for j in range(depth):

		if(j % 2 == 0):
			# Rzz - Rx block
			for i in range(n_spins):
				circuit.rx(p[count], i)
				count = count + 1

			circuit.barrier()

		if(j % 2 == 1):

			for i in range(n_spins):
				circuit.ry(p[count], i)
				count = count + 1

	return circuit

#==========================================

def general_ansatz(n_spins, depth, p):

	count = 0
	circuit = QuantumCircuit(n_spins)

	for i in range(n_spins):
		circuit.h(i)

	for j in range(depth):

		for i in range(n_spins-1):
			circuit.cnot(i, i+1)

		for i in range(n_spins):
			circuit.rz(p[count], i)
			circuit.rx(p[count+1], i)
			circuit.rz(p[count+2], i)
			count += 3
	
	return circuit

#==========================================

def custom_ansatz(p):

	n_spins = 3
	count = 0
	depth = 2

	t_order = 2

	# def c_gate(qubit0,qubit1,p0,p1,p2,p3,p4):
	# 	circuit.rx(p0,qubit0)
	# 	circuit.ry(p1,qubit0)
	#
	# 	circuit.rx(p2,qubit1)
	# 	circuit.ry(p3,qubit1)
	#
	# 	circuit.rzz(p4,qubit0,qubit1)

	def c_gate(qubit0, qubit1, p0, p1, p2, p3, p4):
		circuit.ry(p0, qubit0)
		circuit.rx(p1, qubit0)

		circuit.ry(p2, qubit1)
		circuit.rx(p3, qubit1)

		circuit.rzz(p4, qubit0, qubit1)

	circuit = QuantumCircuit(n_spins)

	# Prepare initial ground state
	for i in range(n_spins):
		circuit.h(i)

	# Making multiple depths
	for d in range(depth):

		for k in range(1, t_order):
			for r in range(k+1):
				for i in range(n_spins):
					if i % (k+1) == r and i+k < n_spins:
						c_gate(i, i+k, p[count], p[count+1], p[count+2], p[count+3], p[count+4])
						count = count+5
				#circuit.barrier()

		for i in range(n_spins):
			circuit.ry(p[count], i)
			count = count + 1
		#circuit.barrier()

	return circuit

#==========================================


def custom_hweff_ansatz(n_spins, depth, p):

	# n_spins = 3
	count = 0
	# depth = 1
	circuit = QuantumCircuit(n_spins)
	t_order = n_spins

	for i in range(n_spins):
		circuit.h(i)

	for j in range(depth):

		for k in range(1, t_order):
			for r in range(k+1):
				for i in range(n_spins):
					if i % (k+1) == r and i+k < n_spins:
						circuit.rzz(p[count], i, i+k)
						count = count+1

		circuit.barrier()

		if(j % 2 == 0):
			for i in range(n_spins):
				circuit.rx(p[count], i)
				count = count + 1

			circuit.barrier()

		if(j % 2 == 1):

			for i in range(n_spins):
				circuit.ry(p[count], i)
				count = count + 1

			circuit.barrier()

	return circuit

#==========================================

def efficient_SU2(n_spins, depth, p):

	count = 0
	circuit = QuantumCircuit(n_spins)

	for i in range(n_spins):
		circuit.h(i)

	for j in range(depth):

		# Rzz - Rx block
		for i in range(n_spins):
			circuit.ry(p[count], i)
			count = count + 1

		circuit.barrier()

		for i in range(n_spins):
			circuit.rz(p[count], i)
			count = count + 1

		circuit.barrier()

		for i in range(n_spins-1):
			circuit.cnot(i, i+1)
		circuit.cnot(0, n_spins-1)
		
		circuit.barrier()

		
	for i in range(n_spins):
		circuit.ry(p[count], i)
		count = count + 1

	circuit.barrier()

	for i in range(n_spins):
		circuit.rz(p[count], i)
		count = count + 1

	return circuit

#############################################################################################""

# p-VQD

# This class aims to simulate the dynamics of a quantum system
# approximating it with a variational ansatz whose parameters
# are varied in order to follow the unitary evolution

# Useful functions

# def projector_zero(n_qubits):
# 	from qiskit.opflow            import Z, I

# 	prj = (1/np.power(2,n_qubits))*(I+Z)

# 	for a in range(n_qubits-1):
# 		prj = prj^(I+Z)

# 	return prj


def projector_zero(n_qubits):
	# This function create the global projector |00...0><00...0|
	from qiskit.opflow import Z, I

	prj_list = [0.5*(I+Z) for i in range(n_qubits)]
	prj = prj_list[0]

	for a in range(1, len(prj_list)):
		prj = prj ^ prj_list[a]

	return prj

def ei(i,n):
	vi = np.zeros(n)
	vi[i] = 1.0
	return vi[:]



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
		# Make the Hamiltonian a list if it's a single one
		if not isinstance(hamiltonian, list):
			hamiltonian = [hamiltonian]
		self.hamiltonian     = hamiltonian
		self.ansatz_depth 	 = ansatz_depth
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


	def construct_total_circuit(self,time_step):
		## This function creates the circuit that will be used to evaluate overlap and its gradient

		# First, create the Trotter step
		if self.ham_tfunc is not None:

			#former method
			step_h  = self.ham_params*np.array(self.hamiltonian[:len(self.ham_tfunc)])*time_step
			step_h  = np.append(step_h, np.array(self.hamiltonian[len(self.ham_tfunc):])*time_step)

			#Magnus expansion
			# step_h  = self.magnus_params*np.array(self.hamiltonian)
			# step_h  = np.append(step_h, time_step*np.array(self.ham_mag2))

		else:
			step_h  = time_step*np.array(self.hamiltonian)
			print(step_h)

		trotter = PauliTrotterEvolution(reps=1)
		# Total Trotter circuit constructed by summing over the Hamiltonian parts
		if len(step_h)>1:
			U_dt = trotter.convert(step_h[0].exp_i()).to_circuit()
			for j in range(1,len(step_h)):
				U_dt += trotter.convert(step_h[j].exp_i()).to_circuit()
		else:
			U_dt    = trotter.convert(step_h[0][0].exp_i()).to_circuit()


		l_circ  = self.ansatz.assign_parameters({self.params_vec: self.left})
		r_circ  = self.ansatz.assign_parameters({self.params_vec: self.right})

		## Projector
		zero_prj = StateFn(projector_zero(self.hamiltonian[0].num_qubits),is_measurement = True)
		state_wfn = zero_prj @ StateFn(r_circ +U_dt+ l_circ.inverse())


		return state_wfn



	def construct_hamiltonian(self):
		## This function creates the circuit that will be used to measure the time-dependent energy observable
		ham_circuit = np.sum(self.ham_params*np.array(self.hamiltonian[:len(self.ham_tfunc)])) \
					  + np.sum(np.array(self.hamiltonian[len(self.ham_tfunc):]))

		return ham_circuit


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

			mean  = sampled_op.eval().real
			#mean  = np.power(np.absolute(mean),2)
			est_err = 0


			if (not self.instance.is_statevector):
				variance = expectator.compute_variance(sampled_op).real
				est_err  = np.sqrt(variance/self.shots)
				self.njobs += 1

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

	def compute_gradient(self, state_wfn, parameters, shift, expectator, sampler):

		nparameters = len(parameters)
		# build dictionary of parameters to values
		# {left[0]: parameters[0], .. ., right[0]: parameters[0] + shift[0], ...}

		# First create the dictionary for overlap
		# values_dict = [dict(zip(self.right[:] + self.left[:],
        #                   parameters.tolist() + (parameters + shift).tolist()))]

		values_dict = []
		# Then the values for the gradient
		for i in range(nparameters):
			values_dict.append(dict(zip(self.right[:] + self.left[:], parameters.tolist() + (
				parameters + shift + ei(i, nparameters)*np.pi/2.0).tolist())))
			values_dict.append(dict(zip(self.right[:] + self.left[:], parameters.tolist() + (
				parameters + shift - ei(i, nparameters)*np.pi/2.0).tolist())))

		# Now evaluate the circuits with the parameters assigned

		results = []

		for values in values_dict:
			sampled_op = sampler.convert(state_wfn, params=values)

			mean = sampled_op.eval().real
			#mean  = np.power(np.absolute(mean),2)
			est_err = 0

			if (not self.instance.is_statevector):
				variance = expectator.compute_variance(sampled_op).real
				est_err = np.sqrt(variance/self.shots)
				self.njobs += 1

			results.append([mean, est_err])

		# E = np.zeros(2)
		g = np.zeros((nparameters, 2))

		# E[0], E[1] = results[0]

		for i in range(nparameters):
			rplus = results[2*i]  # 0+2i
			rminus = results[1+2*i]  # 1+2i
			# G      = (Ep - Em)/2
			# var(G) = var(Ep) * (dG/dEp)**2 + var(Em) * (dG/dEm)**2
			g[i, :] = (rplus[0]-rminus[0])/2.0, np.sqrt(rplus[1]**2+rminus[1]**2)/2.0

		# self.overlap = E  # evaluate with statevec
		self.gradient = g

		return g

	def compute_overlap(self, state_wfn, parameters, shift, expectator, sampler):

		# nparameters = len(parameters)
		# build dictionary of parameters to values
		# {left[0]: parameters[0], .. ., right[0]: parameters[0] + shift[0], ...}

		# First create the dictionary for overlap
		values_dict = [dict(zip(self.right[:] + self.left[:],
                          parameters.tolist() + (parameters + shift).tolist()))]

		# Then the values for the gradient
		# for i in range(nparameters):
		# 	values_dict.append(dict(zip(self.right[:] + self.left[:], parameters.tolist() + (
		# 		parameters + shift + ei(i, nparameters)*np.pi/2.0).tolist())))
		# 	values_dict.append(dict(zip(self.right[:] + self.left[:], parameters.tolist() + (
		# 		parameters + shift - ei(i, nparameters)*np.pi/2.0).tolist())))

		# Now evaluate the circuits with the parameters assigned

		results = []

		for values in values_dict:
			sampled_op = sampler.convert(state_wfn, params=values)

			mean = sampled_op.eval().real
			#mean  = np.power(np.absolute(mean),2)
			est_err = 0

			if (not self.instance.is_statevector):
				variance = expectator.compute_variance(sampled_op).real
				est_err = np.sqrt(variance/self.shots)

			results.append([mean, est_err])

		E = np.zeros(2)
		# g = np.zeros((nparameters, 2))

		E[0], E[1] = results[0]

		# for i in range(nparameters):
		# 	rplus = results[1+2*i]  # 0+2i
		# 	rminus = results[2+2*i]  # 1+2i
		# 	# G      = (Ep - Em)/2
		# 	# var(G) = var(Ep) * (dG/dEp)**2 + var(Em) * (dG/dEm)**2
		# 	g[i, :] = (rplus[0]-rminus[0])/2.0, np.sqrt(rplus[1]**2+rminus[1]**2)/2.0

		# self.overlap = E
		# self.gradient = g

		return E

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

			mean  = sampled_op.eval()
			mean  = np.power(np.absolute(mean),2)
			est_err = 0


			if (not self.instance.is_statevector):
				variance = expectator.compute_variance(sampled_op).real
				est_err  = np.sqrt(variance/self.shots)
				self.njobs += 1

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

	#def compute_gradient
	#def compute_overlap

	def measure_aux_ops(self,obs_wfn,pauli,parameters,expectator,sampler):

		# This function calculates the expectation value of a given operator

		# Prepare the operator and the parameters
		wfn = StateFn(obs_wfn)
		op  = StateFn(pauli,is_measurement = True)
		values_obs = dict(zip(self.obs_params[:], parameters.tolist()))

		braket = op @ wfn

		grouped    = expectator.convert(braket)
		sampled_op = sampler.convert(grouped,params = values_obs)

		mean_value = sampled_op.eval().real
		est_err = 0

		if (not self.instance.is_statevector):
			variance = expectator.compute_variance(sampled_op).real
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

	def nelder_mead(self, f, simplex,
                alpha=1., gamma=2., rho=-0.5, sigma=0.5):
		'''
			@param f (function): function to optimize, must return a scalar score
				and operate over a numpy array of the same dimensions as x_start
			@param x_start (numpy array): initial position
			@param step (float): look-around radius in initial step
			@no_improv_thr,  no_improv_break (float, int): break after no_improv_break iterations with
				an improvement lower than no_improv_thr
			@max_iter (int): always break after this number of iterations.
				Set it to 0 to loop indefinitely.
			@alpha, gamma, rho, sigma (floats): parameters of the algorithm
				(see Wikipedia page for reference)
			return: tuple (best parameter array, best score)
		'''

		# init
		dim = len(self.shift)
		# prev_best = self.overlap
		# # no_improv = 0
		# res = [[self.shift, prev_best]]

		# for i in range(dim):
		# 	x = copy.copy(self.shift)
		# 	x[i] = x[i] + step
		# 	score = f(x)
		# 	res.append([x, score])

		# # simplex iter
		# iters = 0
		# while 1:
		# order
		# simplex.sort(key=lambda x: -x[1])
		# best = res[0][1]

		# # break after max_iter
		# if max_iter and iters >= max_iter:
		# 	return res[0]
		# iters += 1

		# # break after no_improv_break iterations with no improvement
		# print('...best so far:', best)

		# if best < prev_best - no_improve_thr:
		# 	no_improv = 0
		# 	prev_best = best
		# else:
		# 	no_improv += 1

		# if no_improv >= no_improv_break:
		# 	return res[0]

		# centroid
		x0 = [0.] * dim
		for tup in simplex[:-1]:
			for i, c in enumerate(tup[0]):
				x0[i] += c / (len(simplex)-1)

		# reflection
		xr = x0 + alpha*(x0 - simplex[-1][0])
		rscore = f(xr)
		if simplex[0][1][0] >= rscore[0] > simplex[-2][1][0]:
			del simplex[-1]
			simplex.append([xr, rscore])

		# expansion
		elif rscore[0] > simplex[0][1][0]:
			xe = x0 + gamma*(x0 - simplex[-1][0])
			escore = f(xe)
			if escore[0] > rscore[0]:
				del simplex[-1]
				simplex.append([xe, escore])
			else:
				del simplex[-1]
				simplex.append([xr, rscore])

		else:
			# contraction
			xc = x0 + rho*(x0 - simplex[-1][0])
			cscore = f(xc)
			if cscore[0] > simplex[-1][1][0]:
				del simplex[-1]
				simplex.append([xc, cscore])

			else:
				# reduction
				x1 = simplex[0][0]
				nres = []
				for tup in simplex:
					redx = x1 + sigma*(tup[0] - x1)
					score = f(redx)
					nres.append([redx, score])
				simplex = nres

		simplex.sort(key=lambda x: x[1][0], reverse=True)
		
		return simplex


	def run(self,ths,timestep,n_steps,user_messenger:UserMessenger,obs_dict={},max_iter = 100,opt='sgd',grad='separated_param_shift',initial_point=None):



		## initialize useful quantities once
		expectation = PauliExpectation()
		# if(self.instance.is_statevector):
		# 	expectation = PauliExpectation()
		# if(not self.instance.is_statevector):
		# 	expectation = PauliExpectation()

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

			## And the values of the integral parameters for Magnus
			# ham_integ_values = np.array([[self.ham_integ[i](times[j+1]) - self.ham_integ[i](times[j])
			# 				   for j in range(n_steps)] for i in range(len(self.ham_integ))])
			# magnus_dict = [dict(zip(self.magnus_params[:], ham_integ_values[:,i].tolist())) for i in range(n_steps)]

			## And the H(t) operator for observable measurement
			ham_circuit = self.construct_hamiltonian()
			obs_dict['E(t)'] = ham_circuit.assign_parameters(ham_dict[0])

		tot_steps = 0
		time_slice = 0

		print("Running the algorithm")

		if initial_point != None :
			if len(initial_point["params"][0]) != len(self.parameters):
				print("TypeError: Initial parameters are not of the same size of circuit parameters")
				return

			time_slice = initial_point["time_slice"][0]
			print("\nRestart from: ")
			print("step "+str(time_slice))
			self.parameters = np.array(initial_point["params"][-1])
			self.shift = np.array(initial_point["shifts"][-1][-1])
			# self.gradient = initial_point["gradients"][-1][-1]
			# self.njobs = initial_point["njobs"][-1]

			if len(obs_dict) > 0:
				obs_measure = {}
				obs_error = {}

				for (obs_name, obs_pauli) in obs_dict.items():
					obs_measure[str(obs_name)] = initial_point[str(obs_name)]
					obs_error['err_'+str(obs_name)] = initial_point['err_'+str(obs_name)]

			counter = initial_point["iter_number"]
			initial_fidelities = initial_point["init_F"]
			fidelities = initial_point["final_F"]
			err_fin_fid = initial_point["err_fin_F"]
			err_init_fid = initial_point["err_init_F"]
			params = initial_point["params"]

			interm_fidelities = initial_point["interm_F"]  # save intermediate fidelities
			shifts = initial_point["shifts"]  # save all parameter shifts
			gradients = initial_point["gradients"]  # save all evaluated gradients
			err_grad = initial_point["err_grad"]
			grad_norms = initial_point["norm_grad"]
			njobs = initial_point["njobs"]  # save the number of jobs submitted to hardware

		else:
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

			interm_fidelities = [] #save intermediate fidelities
			shifts = [] #save all parameter shifts
			gradients = [] #save all evaluated gradients
			err_grad = []
			grad_norms = []
			njobs = [] #save the number of jobs submitted to hardware

			params.append(list(self.parameters))


		for i in range(time_slice, n_steps):

			# print('\n================================== \n')
			# print("Time slice:",i+1)
			# print("Shift before optimizing this step:",self.shift)
			# print("Initial parameters:", self.parameters)
			# print('\n================================== \n')

			interm_fid_t = [] #will go inside interm_fidelities
			shifts_t = [list(self.shift)] #will go inside shifts
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

			count = 0
			self.overlap = [0.01,0]
			g_norm = 1

			if opt == 'nelder-mead':
				# init
				dim = len(self.shift)
				# no_improv = 0
				simplex = [[self.shift, self.overlap]]

				for d in range(dim):
					x_nm = copy.copy(self.shift)
					x_nm[d] = x_nm[d] + 0.1
					score = self.compute_overlap(state_wfn_Ht,self.parameters,self.shift,expectation,sampler)
					simplex.append([x_nm, score])

				simplex.sort(key=lambda x: x[1][0], reverse=True)

			if opt == 'adam':
				m = np.zeros(len(self.parameters))
				v = np.zeros(len(self.parameters))

			if opt == 'momentum':
				old_grad = np.zeros(len(self.parameters))
				g        = np.zeros((len(self.parameters),2))

			# Save data of overlaps for analysis
			overlap_history = []

			if opt == 'line search':
				from scipy.optimize import curve_fit
				def f(x, a0, k, x0):
					return a0+(k/2.0)*(x-x0)**2
				# def f(x, a0, A, a, phi, b):
				# 	return a0 + np.exp(-b*x)*A*np.sin(a*x+phi)

				self.overlap = self.compute_overlap(
					state_wfn_Ht, self.parameters, self.shift, expectation, sampler)

				nmesh = 21 #number of points to measure
				optsteps = 8 #number of optimization steps per time step

				for opt_step in range(optsteps):
					print("\n", "opt step: ", opt_step, "\n")
					
					t_mesh = np.linspace(0.1, 5, num=nmesh, endpoint=True)
					E_mesh = []

					g = self.compute_gradient(
						state_wfn_Ht, self.parameters, self.shift, expectation, sampler)
					meas_grad = np.asarray(g[:, 0])

					for t in t_mesh:
						### 1) Construct the total circuit
						### 2) Evaluate only overlap for that circuit (O_t) tuple of value and error
						### 3) E_mesh.append(O_t)
						tempshift = copy.copy(self.shift) + t*meas_grad
						temploss = 1-self.compute_overlap(
							state_wfn_Ht, self.parameters, tempshift, expectation, sampler)
						E_mesh.append(temploss)

					t_opt,E_opt = None,None
					E_average = [E[0] for E in E_mesh]
					E_stdev   = [E[1] for E in E_mesh]
					print(E_average)
					p0 = [min(E_average),1.0,t_mesh[np.argmin(E_average)]]
					# p0 = [1.0, -1.0, t_mesh[np.argmin(E_average)], 0.0, 0.0]
					try:
						p_opt,p_cov = curve_fit(
							f,t_mesh,E_average,p0=p0,sigma=E_stdev,absolute_sigma=True,method='lm',maxfev=1600)
					except RuntimeError:
						print("\n Proceeding to next time step... \n")
						break

					t_opt = [p_opt[2],np.sqrt(p_cov[2,2])]
					# print(t_opt)
					E_opt = [p_opt[0],np.sqrt(p_cov[0,0])]
					# E(t) = E(x0-t*g) ~ a* x^2 + bx + c
					# return t_mesh,E_mesh,t_opt,E_opt

					# plt.figure()
					# plt.plot(t_mesh, E_average, 'o')
					# plt.vlines(t_opt[0], ymin=0, ymax=1, color='black')
					# plt.show()

					self.shift = self.shift + t_opt[0]*meas_grad
					self.overlap = E_opt

			else:
				while self.overlap[0] < ths and count < max_iter: # and g_norm > len(params)*8e-6:
					# print("Shift optimizing step:",count+1)
					count = count +1

					if opt == 'momentum':
						old_grad = np.asarray(g[:,0])
					## Measure energy and gradient
					if opt != 'nelder-mead':

						if grad == 'separated_param_shift':
							overlap_backend = Aer.get_backend('statevector_simulator')
							overlap_instance = QuantumInstance(backend=overlap_backend, shots=1)
							overlap_sampler = CircuitSampler(overlap_instance)
							self.overlap = self.compute_overlap(state_wfn_Ht,self.parameters,self.shift,expectation,overlap_sampler)
							g   = self.compute_gradient(
								state_wfn_Ht, self.parameters, self.shift, expectation, sampler)
						if grad == 'param_shift':
							E,g = self.compute_overlap_and_gradient(state_wfn_Ht,self.parameters,self.shift,expectation,sampler)
						if grad == 'spsa':
							E,g = self.compute_overlap_and_gradient_spsa(state_wfn_Ht,self.parameters,self.shift,expectation,sampler,count)

					tot_steps = tot_steps+1

					if count == 1:
						initial_fidelities.append(self.overlap[0])
						err_init_fid.append(self.overlap[1])

						# shifts_t.append(list(self.shift))

					# print('Overlap',self.overlap)
					# print('Gradient', self.gradient[:, 0])

					overlap_history.append(self.overlap[0])

					interm_fid_t.append(self.overlap[0])

					if opt == 'nelder-mead':
						print("\n Nelder-Mead \n")
						def func_nm(x): return self.compute_overlap(
							state_wfn_Ht, self.parameters, x, expectation, sampler)
						simplex = self.nelder_mead(func_nm, simplex)
						self.shift = np.asarray(simplex[0][0])
						self.overlap = np.asarray(simplex[0][1])
						print("new overlap: ", self.overlap)

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

					shifts_t.append(list(self.shift))

					#Norm of the gradient
					if opt != 'nelder-mead':
						g_vec = np.asarray(g[:,0])
						g_norm = np.linalg.norm(g_vec)

						grad_t.append(list(self.gradient[:, 0]))
						err_grad_t.append(list(self.gradient[:, 1]))
						grad_norms_t.append(g_norm)

					# user_messenger.publish({"grad_t": grad_t, "err_grad_t": err_grad_t, "t_step": i+1})


			# Evaluate final overlap
			if grad == 'separated_param_shift' and opt != 'nelder-mead' and opt != 'line search':
				self.overlap = self.compute_overlap(
					state_wfn_Ht, self.parameters, self.shift, expectation, overlap_sampler)
				interm_fid_t.append(self.overlap[0])

			# Select best configuration
			# best_index = np.argmax(interm_fid_t)
			# self.shift = np.array(shifts_t[best_index])
			# self.overlap = np.array([interm_fid_t[best_index], 0])

			# Update parameters


			# print('\n---------------------------------- \n')

			# print("Shift after optimizing:",self.shift)
			# print("New parameters:"        ,self.parameters + self.shift)

			# print("New overlap: "          ,self.overlap[0])

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

			interm_fidelities.append(list(interm_fid_t))
			shifts.append(list(shifts_t))

			gradients.append(list(grad_t))
			err_grad.append(list(err_grad_t))
			grad_norms.append(list(grad_norms_t))

			njobs.append(self.njobs)

			params.append(list(self.parameters))

			# Save data on file

			log_data = {}
			if len(obs_dict) > 0:
				for (obs_name, obs_pauli) in obs_dict.items():
					log_data[str(obs_name)] = obs_measure[str(obs_name)]
					log_data['err_'+str(obs_name)] = obs_error['err_'+str(obs_name)]

			log_data['init_F'] = initial_fidelities
			log_data['final_F'] = fidelities
			log_data['err_init_F'] = err_init_fid
			log_data['err_fin_F'] = err_fin_fid
			log_data['iter_number'] = counter
			log_data['times'] = times[:i+2]
			log_data['params'] = list(params)
			log_data['tot_steps'] = [tot_steps]

			log_data['interm_F'] = list(interm_fidelities)
			log_data['shifts'] = list(shifts)
			log_data['gradients'] = list(gradients)
			log_data['err_grad'] = list(err_grad)
			log_data['norm_grad'] = list(grad_norms)
			log_data['njobs'] = list(njobs)

			log_data['time_slice'] = [i+1]

			user_messenger.publish(log_data)




		# print("Total measurements:",tot_steps)
		# print("Measure per step:", tot_steps/n_steps)

		# print("overlap_history ", overlap_history)

		# Save data on file

		# log_data = {}
		# if len(obs_dict) > 0:
		# 	for (obs_name,obs_pauli) in obs_dict.items():
		# 		log_data[str(obs_name)]        = obs_measure[str(obs_name)]
		# 		log_data['err_'+str(obs_name)] = obs_error['err_'+str(obs_name)]

		# log_data['init_F']      = initial_fidelities
		# log_data['final_F']     = fidelities
		# log_data['err_init_F']  = err_init_fid
		# log_data['err_fin_F']   = err_fin_fid
		# log_data['iter_number'] = counter
		# log_data['times']       = times
		# log_data['params']      = list(params)
		# log_data['tot_steps']   = [tot_steps]

		# log_data['interm_F'] = list(interm_fidelities)
		# log_data['shifts'] = list(shifts)
		# log_data['gradients'] = list(gradients)
		# log_data['err_grad'] = list(err_grad)
		# log_data['norm_grad'] = list(grad_norms)
		# log_data['njobs'] = list(njobs)

		# log_data['overlap_history'] = overlap_history

		# json.dump(log_data, open( filename,'w+'))
		return log_data

###############################################################################

#MAIN

def main(backend, user_messenger, **kwargs):
	"""Main entry point of the program.

    Args:
        backend: Backend to submit the circuits to.
        user_messenger: Used to communicate with the program consumer.
        kwargs: User inputs.
    """

	#user inputs
	nqubits = kwargs.pop('nqubits', 3)
	V       = kwargs.pop('coupling', -1.0)
	g 		= kwargs.pop('field', -1.0)
	tmax 	= kwargs.pop('tmax', 3.0)
	dt 		= kwargs.pop('dt', 0.05)
	n_steps = int(tmax/dt)
	ths 	= kwargs.pop('threshold', 0.99999)
	depth 	= kwargs.pop('depth', 1)
	anstz   = kwargs.pop('ansatz', 'hweff')
	NN 		= kwargs.pop('NN', 1) #nearest neighbor
	HamDict = kwargs.pop('hamiltonian', ['hzz', 'hx'])
	# ansatz  = kwargs.pop('ansatz', hweff_ansatz_adiab)
	# H_tfunc_bool = kwargs.pop('h_tfunc', [False for i in range(len(HamDict))])
	shots 	= kwargs.pop('shots', 8000)
	opt     = kwargs.pop('optimizer', 'sgd')
	obs     = kwargs.pop('obs', {'E': generate_ising(nqubits, V, g)})
	max_iter= kwargs.pop('iterations', 30)
	initial_point = kwargs.pop('initial_point', None)
	measurement_error_mitigation = kwargs.get(
		"measurement_error_mitigation", False)

	H = []
	for label in HamDict:
		if label == 'hzz':
			H.append(generate_ising_Hzz(nqubits, V))
		elif label == 'hx':
			H.append(generate_ising_Hx(nqubits, g))

	H_tfunc = [lambda x: x/tmax]

	if anstz == 'hweff':
		if NN == 1:
			ansatz = hweff_ansatz_adiab
			ex_params = np.zeros(depth*nqubits + depth*(nqubits-1))
		else:
			ansatz = custom_hweff_ansatz
			ex_params = np.zeros(depth*((nqubits-1) + (nqubits-2) + nqubits))

	elif anstz == 'su2':
		ansatz = efficient_SU2
		ex_params = np.zeros(nqubits*2 + depth*nqubits*2)

	elif anstz == 'general':
		ansatz = general_ansatz
		ex_params = np.zeros(nqubits*3*depth)

	if opt == 'line search':
		np.random.seed(99)
		shift = np.random.normal(0, 0.1, len(ex_params))
	else:
		shift = np.zeros(len(ex_params))
	# shift = np.zeros(len(ex_params))

	# set up quantum instance
	if measurement_error_mitigation:
		instance = QuantumInstance(
			backend,
			shots=shots,
			measurement_error_mitigation_shots=shots,
			measurement_error_mitigation_cls=CompleteMeasFitter,
		)
	else:
		instance = QuantumInstance(backend, shots=shots)

	algo = pVQD(H, ansatz, depth, ex_params,
	            shift, instance, shots, H_tfunc)

	begin = time.time()
	output = algo.run(ths, dt, n_steps, user_messenger, obs_dict=obs, max_iter=max_iter, opt=opt,
					 grad='separated_param_shift', initial_point=initial_point)
	end = time.time() - begin
	output["exec_time"] = end

	return output

