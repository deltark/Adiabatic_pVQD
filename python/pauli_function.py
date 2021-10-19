import numpy as np
import matplotlib.pyplot as plt


from qiskit.quantum_info 			  import Pauli
from qiskit.aqua.operators 			  import PauliOp, SummedOp

def generate_pauli(idx_x,idx_z,n):
	'''
	Args:
		n (integer)
		idx (list)
	Returns:
		tensor product of Pauli operators acting on qubits in idx
	'''

	xmask = [0]*n
	zmask = [0]*n
	for i in idx_x : xmask[i] = 1
	for i in idx_z : zmask[i] = 1

	a_x = np.asarray(xmask,dtype =np.bool)
	a_z = np.asarray(zmask,dtype =np.bool)

	return Pauli(a_z,a_x)


def generate_ising_pbc(n_spins,coup,field):
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

	int_list.append(generate_pauli([],[0,n_spins-1],n_spins))

	if(n_spins>2):
		for i in range(n_spins-1):
			int_list.append(generate_pauli([],[i,i+1],n_spins))

	for i in range(n_spins):
		field_list.append(generate_pauli([i],[],n_spins))

	int_coeff = [coup]*len(int_list)
	field_coeff = [field]*len(field_list)

	H = PauliOp(int_list[0],int_coeff[0])

	for i in range(1,len(int_list)):
		H = H + PauliOp(int_list[i],int_coeff[i])

	for i in range(len(field_list)):
		H = H + PauliOp(field_list[i],field_coeff[i])

	return H


def generate_ising(n_spins,coup,field):
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
		int_list.append(generate_pauli([],[i,i+1],n_spins))

	for i in range(n_spins):
		field_list.append(generate_pauli([i],[],n_spins))

	int_coeff = [coup]*len(int_list)
	field_coeff = [field]*len(field_list)

	H = PauliOp(int_list[0],int_coeff[0])

	for i in range(1,len(int_list)):
		H = H + PauliOp(int_list[i],int_coeff[i])

	for i in range(len(field_list)):
		H = H + PauliOp(field_list[i],field_coeff[i])

	return H

def generate_ising_Hzz(n_spins, coup):

	int_list = []

	for i in range(n_spins-1):
		int_list.append(generate_pauli([],[i,i+1],n_spins))

	int_coeff = [coup]*len(int_list)

	H = PauliOp(int_list[0],int_coeff[0])

	for i in range(1,len(int_list)):
		H = H + PauliOp(int_list[i],int_coeff[i])

	return H

def generate_ising_Hx(n_spins, field):

	field_list = []

	for i in range(n_spins):
		field_list.append(generate_pauli([i],[],n_spins))

	field_coeff = [field]*len(field_list)

	H = PauliOp(field_list[0],field_coeff[0])

	for i in range(1, len(field_list)):
		H = H + PauliOp(field_list[i],field_coeff[i])

	return H

def generate_magnus_2(n_spins, coup, field):

	listYZ = []
	listZY = []

	for i in range(n_spins-1):
		listYZ.append(generate_pauli([i],[i,i+1],n_spins))
		listZY.append(generate_pauli([i+1],[i,i+1],n_spins))

	coeff = [coup*field]*len(listYZ)

	H2 = PauliOp(listYZ[0],coeff[0])
	H2 += PauliOp(listZY[0],coeff[0])

	for i in range(1,len(listYZ)):
		H2 += PauliOp(listYZ[i],coeff[i])
		H2 += PauliOp(listZY[i],coeff[i])

	return H2
