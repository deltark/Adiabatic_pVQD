import numpy as np
import scipy.linalg as LA

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
