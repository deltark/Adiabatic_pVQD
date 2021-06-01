import numpy as np
import functools
import itertools
import matplotlib.pyplot as plt
from scipy   import  linalg as LA
import json

from qiskit import IBMQ, Aer
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit.circuit import ParameterVector

from qiskit.utils                     import QuantumInstance
from qiskit.opflow                    import Z, I, X, Y
from qiskit.opflow         			  import PauliOp, SummedOp

import time

from pauli_function import *
from pVQD			import *
from ansatze        import *

# Initialize system parameters for Ising

spins   = 3
V       = -0.25
g       = -1.0
dt      = 0.05
tmax    = 2.0
n_steps = int(tmax/dt)

# Algorithm parameters

ths = 0.9999999
depth = 2


### Example circ

# ex_params = np.zeros((depth+1)*spins +depth*(spins-1)) #hweff_ansatz
# ex_params = np.zeros(depth*spins + depth*(spins-1)) #hweff_ansatz_adiab
# ex_params = np.zeros(depth*(5*(spins-1) + 5*(spins-2) + spins)) #custom_ansatz t_order=3
# ex_params = np.zeros(depth*(5*(spins-1) + spins)) #custom_ansatz t_order=2
ex_params = np.zeros(depth*((spins-1) + (spins-2) + spins)) #custom_hweff_ansatz t_order=3
wfn = custom_hweff_ansatz(spins,depth,ex_params)


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
shots = 1
# shots = 10000000
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
opt  = 'sgd'
# Choose how to estimate the gradient on hardware: 'param_shift', 'spsa'
grad = 'param_shift'
# Choose which type of cost function use: 'global', 'local'
cost = 'global'


algo = pVQD(hamiltonian   = H,
			ansatz        = custom_hweff_ansatz,
			ansatz_reps   = depth,
			parameters    = ex_params,
			initial_shift = shift,
			instance      = instance,
			shots         = shots,
			ham_tfunc     = H_tfunc)

begin=time.time()
algo.run(ths,dt,n_steps,
	     obs_dict      = obs,
	     filename      = 'data/time_compar_new.dat',
	     max_iter      = 30,
	     opt           = opt,
	     cost_fun      = cost,
	     grad          = grad,
	     initial_point = initial_point)
print("runtime:")
print(time.time()-begin)
