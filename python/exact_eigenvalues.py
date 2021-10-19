import numpy as np
import scipy.linalg as LA
import matplotlib.pyplot as plt

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



def second_smallest(numbers):
    m1 = m2 = float('inf')
    for x in numbers:
        if x <= m1:
            m1, m2 = x, m1
        elif x < m2:
            m2 = x
    return m2

def X():
    return np.array([[0, 1],[1, 0]], dtype = np.complex)

def Z():
    return np.array([[1, 0],[0, -1]], dtype = np.complex)

def kronecker(oparray):
    while len(oparray) > 1:
        oparray[0] = np.kron(oparray[0], oparray[1])
        oparray.pop(1)
    return oparray[0]

def ising(nqubits, J, G):
    resZZ = np.zeros([2**nqubits, 2**nqubits], dtype = np.complex)
    resX = np.zeros([2**nqubits, 2**nqubits], dtype = np.complex)

    #interaction part
    for i in range(nqubits-1):
        oplist = [np.eye(2, dtype = np.complex) for j in range(nqubits)]
        oplist[i] = Z()
        oplist[i+1] = Z()
        resZZ += kronecker(oplist)
    oplist = [np.eye(2, dtype = np.complex) for j in range(nqubits)]
    # oplist[0] = Z()
    # oplist[nqubits-1] = Z()
    # resZZ += kronecker(oplist)
    resZZ *= J

    #X part
    for i in range(nqubits):
        oplist = [np.eye(2, dtype = np.complex) for j in range(nqubits)]
        oplist[i] = X()
        resX += kronecker(oplist)
    resX *= G

    return resZZ+resX

def hamiltonian_eig(nqubits, J, G):
    H = ising(nqubits, J, G)
    return LA.eigh(H, eigvals_only=True)


# Create the Hamiltonian of the system

spins   = 3
V       = -1.0
g       = -1.0

# dt = 0.05
tmax = 3.0
# Nt = int(tmax/dt)
Nt = 60
dt = tmax/Nt
times_norm = [i/Nt for i in range(Nt+1)]

eigenvalsGS = np.zeros(Nt+1)
eigenvalsFE = np.zeros(Nt+1)

for i in range(Nt+1):
	# count +=1

	print("\n-------------------")
	print("Time normalized: "+str(times_norm[i]))
	print("-------------------\n")

	E = hamiltonian_eig(spins, times_norm[i]*V, g)
	eigenvalsGS[i] = np.min(E)
	eigenvalsFE[i] = second_smallest(E)

np.save('data/exactGS.npy',eigenvalsGS)
np.save('data/exactFE.npy',eigenvalsFE)
