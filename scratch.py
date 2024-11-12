from qiskit import QuantumCircuit, transpile, assemble
from qiskit.quantum_info import random_clifford

import time, random

import matplotlib.pyplot as plt
import numpy as np
import ccl


n=4
zero = ccl.StabilizerState.zero_state(n)
circ = ccl.CliffordCircuit(n)
circ.h(0)
clifford = circ.compile()
# clifford = ccl.Clifford.random_clifford(n)
evolved_state = zero.evolve(clifford)
probabilties = [evolved_state.get_probability(ccl.int2basestr(i,2,n)) for i in range(2**n)]
probabilties_exp = evolved_state.to_qiskit().probabilities()

print(np.array_equal([0,0,1],[0,0,1.]))

