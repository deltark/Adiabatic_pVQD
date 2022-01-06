import numpy as np
import functools
import itertools
import matplotlib.pyplot as plt
from scipy   import  linalg as LA
import json

import time

from qiskit import IBMQ, Aer
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister

from qiskit.aqua                      import QuantumInstance
from qiskit.aqua.operators            import Z, I, X, Y
from qiskit.aqua.operators 			  import PauliOp, SummedOp

from qiskit.providers.aer.noise import NoiseModel

# Build noise model from backend properties
# provider = IBMQ.enable_account("210d86ba5de279939c25b694d52b9125670e5de499f4daaa632e6068936f3f74ece4503dae344774501af1aded7c068bdc948a6062c1983c41bfeb77a1501dcf")
# backend = provider.get_backend('ibmq_santiago')
# noise_model = NoiseModel.from_backend(backend)

if __name__ == "__main__":
	from pauli_function    import *
	from hybrid_pVQD			   import *
	from ansatze           import *

	#tool
	# def second_smallest(numbers):
	#     m1, m2 = float('inf'), float('inf')
	#     for x in numbers:
	#         if x <= m1:
	#             m1, m2 = x, m1
	#         elif x < m2:
	#             m2 = x
	#     return m2

	# Initialize system parameters for Ising

	spins   = 3
	V       = -1.0
	g       = -1.0
	dt      = 0.05
	tmax    = 3.0
	n_steps = int(tmax/dt)
	# n_steps = 1
	# dt = tmax/n_steps

	# Compute the exact ground state of the Hamiltonian
	# Heig = hamiltonian_eig(spins, V, g)
	# exactGS = np.array([np.min(Heig), second_smallest(Heig)])
	# np.save('data/exactGS_J110.npy', exactGS)

	# Algorithm parameters

	# ths = 0.9999999
	ths = 0.99999
	depth = 1


	### Example circ

	# ex_params = np.zeros((depth+1)*spins +depth*(spins-1)) #hweff_ansatz
	ex_params = np.zeros(depth*spins + depth*(spins-1)) #hweff_ansatz_adiab
	# ex_params = np.zeros(depth*(5*(spins-1) + 5*(spins-2) + spins)) #custom_ansatz t_order=3
	# ex_params = np.zeros(depth*(5*(spins-1) + spins)) #custom_ansatz t_order=2
	# ex_params = np.zeros(depth*((spins-1) + (spins-2) + spins)) #custom_hweff_ansatz t_order=3
	# wfn = custom_hweff_ansatz(ex_params)


	### Shift
	# shift  = np.array(len(ex_params)*[0.01])
	# shift = np.random.normal(0.0, 0.01, len(ex_params))
	shift = np.zeros(len(ex_params))

	print("Initial shift:",shift)


	### Generate the Hamiltonian
	Hzz = generate_ising_Hzz(spins, V)
	Hx  = generate_ising_Hx(spins, g)
	H = [Hzz, Hx]
	H_tfunc = [lambda x : x/tmax]

	### indefinite integral for the Magnus expansion
	# H_integral = [lambda x : x**2/(2*tmax), lambda x : x]
	# ### second order Magnus, I have no idea how to make this more general aaaa
	# H_m2 = dt**2/(6*tmax)*np.array([generate_magnus_2(spins, V, g)])
	# print(H_m2)

	# print(wfn)
	# print(H)

	### Backend
	# shots = 400000
	shots = 8000
	# shots = 1

	IBMQ.load_account()
	provider = IBMQ.get_provider(
		hub='ibm-q-research-2', group='epfl-2', project='main')
	# backend = provider.get_backend('ibmq_qasm_simulator')
	backend = provider.get_backend('ibm_lagos')

	# backend = Aer.get_backend('qasm_simulator')
	# backend  = Aer.get_backend('qasm_simulator', noise_model=noise_model)
	# instance = QuantumInstance(backend=backend,shots=shots)

	### Prepare the observables to measure
	obs = {}
	# Magnetization

	# for i in range(spins):
	# 	obs['Sz_'+str(i)]      = PauliOp(generate_pauli([],[i],spins),1.0)
	# 	obs['Sx_'+str(i)]      = PauliOp(generate_pauli([i],[],spins),1.0)
	# 	obs['Sy_'+str(i)]      = PauliOp(generate_pauli([i],[i],spins),1.0)

	#Energy
	obs['E'] = generate_ising(spins,V,g)

	# for (name,pauli) in obs.items():
	# 	print(name)
	# 	print(pauli)


	### Initialize the algorithm

	# Choose a specific set of parameters
	initial_point = None

	# Choose the gradient optimizer: 'sgd', 'adam'
	gradient = 'sgd'


	algo = pVQD(H,hweff_ansatz_adiab,depth,ex_params,shift,provider,backend,shots,H_tfunc)

	begin = time.time()
	algo.run(ths,dt,n_steps, obs_dict = obs,filename= 'data/VQD/hybrid_test_casablanca.dat', max_iter = 30)
	print("Total time : ", time.time()-begin)


#TODO
# étude sur les 3 régimes T grand, ~1/Delta^2, petit, avec abscisse rescalée t/T: exact + trotter + pVQD
# étude ansatz
# comparer avec D-wave ?


#do plots t/T
#correlation exact sim for more qubits
