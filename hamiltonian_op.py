from pauli_function import *

class IsingHamiltonian:

    def __init__(self, nqubits, coup, field, tmax, timedep):
        self.num_qubits = nqubits
        self.coup = coup
        self.field = field
        self.time = 0.0
        self.tmax = tmax
        self.timedep = timedep #bool
        if timedep:
            self.operator = generate_ising(nqubits, 0.0, field)
        else:
            self.operator = generate_ising(nqubits, coup, field)

    def update_time(self, newtime):
        if self.timedep:
            self.time = newtime
            ft = self.coup*self.time/self.tmax
            self.operator = generate_ising(self.num_qubits, ft, self.field)
        else:
            print('This Hamiltonian is not time-dependent')
        return self.operator
