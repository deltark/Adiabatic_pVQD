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

from pauli_function import *
from pVQD			import *
from ansatze        import *
from hamiltonian_op import IsingHamiltonian

# Initialize system parameters for Ising

spins   = 3
V       = 0.25
g       = 1.0
dt      = 0.05
n_steps = 40
tmax = 2.0

# Algorithm parameters

ths = 0.99999
depth = 2


### Example circ

ex_params = np.zeros((depth+1)*spins +depth*(spins-1))
wfn = hweff_ansatz(ex_params)


### Shift
shift  = np.array(len(ex_params)*[0.01])

print("Initial shift:",shift)


### Generate the Hamiltonian
H = IsingHamiltonian(spins, V, g, tmax, timedep=True)

print(wfn)
print(H)

### Backend
shots = 1
backend  = Aer.get_backend('statevector_simulator')
instance = QuantumInstance(backend=backend,shots=shots)

### Prepare the observables to measure
obs = {}
# Magnetization

for i in range(spins):
	obs['Sz_'+str(i)]      = PauliOp(generate_pauli([],[i],spins),1.0)
	obs['Sx_'+str(i)]      = PauliOp(generate_pauli([i],[],spins),1.0)
	obs['Sy_'+str(i)]      = PauliOp(generate_pauli([i],[i],spins),1.0)

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


algo = pVQD(H,hweff_ansatz,ex_params,shift,instance,shots)
algo.run(ths,dt,n_steps, obs_dict = obs,filename= 'data/J_param/025.dat', max_iter = 50)
