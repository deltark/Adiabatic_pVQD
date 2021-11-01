## These file will contain all the ansatze used for variational quantum simulation

from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister


#=========================================

def hweff_ansatz(p):
	n_spins = 3
	count = 0
	circuit = QuantumCircuit(n_spins)
	depth = 2

	for j in range(depth):

		if(j%2 == 0):
			# Rx - Rzz block
			for i in range(n_spins):
				circuit.rx(p[count],i)
				count = count +1

			circuit.barrier()

			for i in range(n_spins-1):
				circuit.rzz(p[count],i,i+1)
				count = count +1

			circuit.barrier()

		if(j%2 == 1):
			for i in range(n_spins):
				circuit.ry(p[count],i)
				count = count +1

			circuit.barrier()

			for i in range(n_spins-1):
				circuit.rzz(p[count],i,i+1)
				count = count +1

			circuit.barrier()

	# Final block to close the ansatz
	if (depth%2 == 1):
		for i in range(n_spins):
				circuit.ry(p[count],i)
				count = count +1
	if (depth%2 == 0):
		for i in range(n_spins):
				circuit.rx(p[count],i)
				count = count +1

	return circuit

#==========================================

def hweff_ansatz_adiab(p):
	n_spins = 3
	count = 0
	circuit = QuantumCircuit(n_spins)
	depth = 3

	for i in range(n_spins):
		circuit.h(i)

	for j in range(depth):

		if(j%2 == 0):
			# Rzz - Rx block
			for i in range(n_spins-1):
				circuit.rzz(p[count],i,i+1)
				count = count +1

			circuit.barrier()

			for i in range(n_spins):
				circuit.rx(p[count],i)
				count = count +1

			circuit.barrier()

		if(j%2 == 1):
			for i in range(n_spins-1):
				circuit.rzz(p[count],i,i+1)
				count = count +1

			circuit.barrier()

			for i in range(n_spins):
				circuit.ry(p[count],i)
				count = count +1

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

	def c_gate(qubit0,qubit1,p0,p1,p2,p3,p4):
		circuit.ry(p0,qubit0)
		circuit.rx(p1,qubit0)

		circuit.ry(p2,qubit1)
		circuit.rx(p3,qubit1)

		circuit.rzz(p4,qubit0,qubit1)

	circuit = QuantumCircuit(n_spins)

	# Prepare initial ground state
	for i in range(n_spins):
		circuit.h(i)

	# Making multiple depths
	for d in range(depth):

		for k in range(1,t_order):
			for r in range(k+1):
				for i in range(n_spins):
					if i%(k+1) == r and i+k < n_spins:
						c_gate(i,i+k,p[count],p[count+1],p[count+2],p[count+3],p[count+4])
						count = count+5
				#circuit.barrier()

		for i in range(n_spins):
			circuit.ry(p[count],i)
			count = count +1
		#circuit.barrier()

	return circuit

#==========================================

def custom_hweff_ansatz(p):

	n_spins = 3
	count = 0
	depth = 1
	circuit = QuantumCircuit(n_spins)
	t_order = n_spins

	for i in range(n_spins):
		circuit.h(i)

	for j in range(depth):

		for k in range(1,t_order):
			for r in range(k+1):
				for i in range(n_spins):
					if i%(k+1) == r and i+k < n_spins:
						circuit.rzz(p[count], i, i+k)
						count = count+1

		circuit.barrier()

		if(j%2 == 0):
			for i in range(n_spins):
				circuit.rx(p[count],i)
				count = count +1

			circuit.barrier()

		if(j%2 == 1):

			for i in range(n_spins):
				circuit.ry(p[count],i)
				count = count +1
			
			circuit.barrier()

	return circuit
