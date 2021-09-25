import numpy as np
import json
import functools
import itertools
import matplotlib.pyplot as plt
from scipy   import  linalg as LA



import qiskit
from qiskit                           import QuantumCircuit, ClassicalRegister, QuantumRegister, Aer, execute
from qiskit.aqua                      import QuantumInstance
from qiskit.quantum_info 			  import Pauli
from qiskit.aqua.operators 			  import PauliOp, SummedOp
from qiskit.aqua.operators.evolutions import Trotter, PauliTrotterEvolution

from qiskit.aqua.operators.state_fns     import CircuitStateFn
from qiskit.aqua.operators.primitive_ops import CircuitOp
from qiskit.aqua.operators               import X, Y, Z, I, PauliExpectation, CircuitSampler, StateFn, MatrixExpectation
from qiskit.circuit                      import ParameterVector



from pauli_function        import *


# Create the Hamiltonian of the system

spins   = 3
V       = -1.0
g       = -1.0

# dt = 0.05
tmax = 0.1
# Nt = int(tmax/dt)
Nt = 50
dt = tmax/Nt
times_norm = [i/Nt for i in range(Nt+1)]


### Generate the Hamiltonian
H = generate_ising(spins,V,g)

# Hzz = generate_ising_Hzz(spins, V)
# Hx = generate_ising_Hx(spins, g)
# H = [Hzz, Hx]
# Htfunc = [lambda x : x/tmax]

# Hparams = ParameterVector('H', len(Htfunc))
#
# Hvalues = np.array([[Htfunc[i](t) for t in times_norm] for i in range(len(Htfunc))])
# Hdict = [dict(zip(Hparams[:], Hvalues[:,i].tolist())) for i in range(Nt+1)]

### Prepare the observables to measure

obs = {}
# Magnetization

# for i in range(spins):
# 	obs['Sz_'+str(i)]      = PauliOp(generate_pauli([],[i],spins),1.0)
# 	obs['Sx_'+str(i)]      = PauliOp(generate_pauli([i],[],spins),1.0)
# 	obs['Sy_'+str(i)]      = PauliOp(generate_pauli([i],[i],spins),1.0)

# obs['Sz']      = 0.*I^spins
# obs['Sx']      = 0.*I^spins
# obs['Sy']      = 0.*I^spins
obs['E'] = generate_ising(spins, V, g) #H(T)

# for i in range(spins):
# 	obs['Sz']      += 1/spins*PauliOp(generate_pauli([],[i],spins),1.0)
# 	obs['Sx']      += 1/spins*PauliOp(generate_pauli([i],[],spins),1.0)
# 	obs['Sy']      += 1/spins*PauliOp(generate_pauli([i],[i],spins),1.0)


for (name,pauli) in obs.items():
	print(name)
	print(pauli)

# And the container for the results
obs_measure = {}
obs_error   = {}

## Create a quantum instance

shots     = 1
# n         = 25
backend   = Aer.get_backend('statevector_simulator')
instance  = QuantumInstance(backend=backend,shots=shots)


## Define a function to measure observables

def measure_obs(pauli,wfn,instance,shots):
	'''
	pauli: PauliOp
	wfn:   CircuitStateFn()
	'''

	op = StateFn(pauli,is_measurement = True)

	# Evaluate the aux operator given
	braket = op @ wfn
	grouped = PauliExpectation().convert(braket)
	sampled_op = CircuitSampler(instance).convert(grouped)

	mean_value = sampled_op.eval().real
	est_err = 0

	if (not instance.is_statevector):
		variance = PauliExpectation().compute_variance(sampled_op).real
		est_err  = np.sqrt(variance/shots)

	res = [mean_value,est_err]

	return res


### Create the circuit for different times_norm and measure the observables

# count = 0
for i in range(Nt+1):
	# count +=1

	print("\n-------------------")
	print("Time normalized: "+str(times_norm[i]))
	print("-------------------\n")

	#observable of instantaneous energy H(t)
	obs['E(t)'] = generate_ising(spins, times_norm[i]*V, g)

	#### Initialization
	initialize = QuantumCircuit(spins)
	for j in range(spins):
		initialize.h(j)

	# Now let's create the Trotter operator

	for t in times_norm[:i+1]:
		step_h = dt*generate_ising(spins, t*V, g)
		# trotter = PauliTrotterEvolution(reps=count)
		trotter = PauliTrotterEvolution(reps=1)
		initialize += trotter.convert(step_h.exp_i()).to_circuit()

	#### Total circuit

	wfn = CircuitStateFn(initialize)

	for (obs_name,obs_pauli) in obs.items():

		if i == 0 :
			first_measure                   = measure_obs(obs_pauli, wfn, instance, shots)
			obs_measure[str(obs_name)]      = [first_measure[0]]
			obs_error['err_'+str(obs_name)] = [first_measure[1]]


		else:
			run_measure   = measure_obs(obs_pauli, wfn, instance, shots)
			obs_measure[str(obs_name)].append(run_measure[0])
			obs_error['err_'+str(obs_name)].append(run_measure[1])

## Save the measured observables

log_data = {}
log_data['times']  = times_norm

for (obs_name,obs_pauli) in obs.items():
	log_data[str(obs_name)]        = obs_measure[str(obs_name)]
	log_data['err_'+str(obs_name)] = obs_error['err_'+str(obs_name)]


json.dump(log_data, open('data/trotter/T1_dt1.dat','w+'))
