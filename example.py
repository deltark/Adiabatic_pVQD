import numpy as np
import functools
import itertools
import matplotlib.pyplot as plt
from scipy   import  linalg as LA
import json

from qiskit import IBMQ, Aer
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister

from qiskit.aqua                      import QuantumInstance
from qiskit.aqua.operators            import Z, I, X, Y
from qiskit.aqua.operators 			  import PauliOp, SummedOp

if __name__ == "__main__":
	from pauli_function    import *
	from pVQD			   import *
	from ansatze           import *
	from exact_eigenvalues import hamiltonian_eig

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
	V       = -0.8
	g       = -1.0
	dt      = 0.05
	tmax    = 2.0
	n_steps = int(tmax/dt)

	# Compute the exact ground state of the Hamiltonian
	# Heig = hamiltonian_eig(spins, V, g)
	# exactGS = np.array([np.min(Heig), second_smallest(Heig)])
	# np.save('data/exactGS_J110.npy', exactGS)

	# Algorithm parameters

	ths = 0.99999
	depth = 2


	### Example circ

	# ex_params = np.zeros((depth+1)*spins +depth*(spins-1)) #hweff_ansatz
	# ex_params = np.zeros(depth*spins + depth*(spins-1)) #hweff_ansatz_adiab
	# ex_params = np.zeros(depth*(5*(spins-1) + 5*(spins-2) + spins)) #custom_ansatz t_order=3
	# ex_params = np.zeros(depth*(5*(spins-1) + spins)) #custom_ansatz t_order=2
	ex_params = np.zeros(depth*((spins-1) + (spins-2) + spins)) #custom_hweff_ansatz t_order=3
	wfn = custom_hweff_ansatz(ex_params)


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

	print(wfn)
	print(H)

	### Backend
	# shots = 400000
	shots = 1
	backend  = Aer.get_backend('statevector_simulator')
	instance = QuantumInstance(backend=backend,shots=shots)

	### Prepare the observables to measure
	obs = {}
	# Magnetization

	# for i in range(spins):
	# 	obs['Sz_'+str(i)]      = PauliOp(generate_pauli([],[i],spins),1.0)
	# 	obs['Sx_'+str(i)]      = PauliOp(generate_pauli([i],[],spins),1.0)
	# 	obs['Sy_'+str(i)]      = PauliOp(generate_pauli([i],[i],spins),1.0)

	#Energy
	obs['E'] = generate_ising(spins,V,g)

	for (name,pauli) in obs.items():
		print(name)
		print(pauli)


	### Initialize the algorithm

	# Choose a specific set of parameters
	initial_point = None

	# Choose the gradient optimizer: 'sgd', 'adam'
	gradient = 'sgd'


	algo = pVQD(H,custom_hweff_ansatz,ex_params,shift,instance,shots,H_tfunc)
	algo.run(ths,dt,n_steps, obs_dict = obs,filename= 'data/VQD_J080_dt005_ths5.dat', max_iter = 30, opt = 'sgd')
