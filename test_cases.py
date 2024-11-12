import numpy as np
import random
import copy
import ccl
from collections import Counter


def test_clifford_class():
    """
    Test the ccl.Clifford class methods.
    """
    print("Testing ccl.Clifford class methods.")
    
    # Test 1: Create a Clifford object and test __init__ and __repr__
    print("Test 1: Create a Clifford object and test __repr__")
    n = 3
    clifford = ccl.Clifford.random_clifford(n)
    print("Clifford representation:")
    print(clifford)
    print("Passed")

    print("All tests for ccl.Clifford class passed.\n")

def test_pauli_class():
    """
    Test the ccl.Pauli class methods.
    """
    print("Testing ccl.Pauli class methods.")
    
    # Test 1: Create Pauli objects and test __init__ and __repr__
    print("Test 1: Create Pauli objects and test __repr__")
    pauli_x = ccl.Pauli.from_pystr('X')
    pauli_y = ccl.Pauli.from_pystr('Y')
    pauli_z = ccl.Pauli.from_pystr('Z')
    pauli_i = ccl.Pauli.from_pystr('I')
    print(f"Pauli X: {pauli_x}")
    print(f"Pauli Y: {pauli_y}")
    print(f"Pauli Z: {pauli_z}")
    print(f"Pauli I: {pauli_i}")
    print("Passed")
    
    # Test 2: Test multiply method
    print("Test 2: Test multiply method")
    # X * Y = iZ
    result = pauli_x.multiply(pauli_y)
    assert str(result) == 'iZ', f"Expected iZ, got {result}"
    # Y * Z = iX
    result = pauli_y.multiply(pauli_z)
    assert str(result) == 'iX', f"Expected iX, got {result}"
    # Z * X = iY
    result = pauli_z.multiply(pauli_x)
    assert str(result) == 'iY', f"Expected iY, got {result}"
    print("Passed")
    
    # Test 3: Test conjugate method
    print("Test 3: Test conjugate method")
    # Conjugate Pauli X by Hadamard (H X H = Z)
    circuit = ccl.CliffordCircuit(1)
    circuit.h(0)
    clifford = circuit.compile()
    pauli_x = ccl.Pauli.from_pystr('X')
    conjugated_pauli = pauli_x.conjugate(clifford)
    assert str(conjugated_pauli) == '+Z', f"Expected +Z, got {conjugated_pauli}"
    print("Passed")

    # Test 4: Test commutes method
    print("Test 4: Test commutes method")

    XYYZ = ccl.Pauli.from_pystr("XYYZ")
    XIIX = ccl.Pauli.from_pystr("XIIX")
    ZZZZ = ccl.Pauli.from_pystr("ZZZZ")
    ZIII = ccl.Pauli.from_pystr("ZIII")

    assert not XYYZ.commutes(ZIII), "XYYZ and ZIII should not commute"
    assert ZIII.commutes(ZZZZ), "ZIII and ZZZZ should commute"
    assert ZZZZ.commutes(XIIX), "ZZZZ and XIIX should commute"
    assert not XYYZ.commutes(ZZZZ), "XYYZ and ZZZZ should not commute"
    print("Passed")
    
    
    print("All tests for ccl.Pauli class passed.\n")

def test_stabilizer_state_class():
    """
    Test the ccl.StabilizerState class methods.
    """
    print("Testing ccl.StabilizerState class methods.")
    
    # Test 1: Initialize StabilizerState and test __repr__
    print("Test 1: Initialize StabilizerState and test __repr__")
    n = 3
    stabilizer_state = ccl.StabilizerState.zero_state(n)
    print("StabilizerState representation:")
    print(stabilizer_state)
    print("Passed")
    
    # Test 2: Test get_tableau method
    print("Test 2: Test get_tableau method")
    tableau = stabilizer_state.get_tableau()
    expected_tableau = np.array([[0, 0, 0, 1, 0, 0, 0],
                                 [0, 0, 0, 0, 1, 0, 0],
                                 [0, 0, 0, 0, 0, 1, 0]])
    assert np.array_equal(tableau, expected_tableau), f"Expected tableau {expected_tableau}, got {tableau}"
    print("Passed")
    
    # Test 3: Test evolve method
    print("Test 3: Test evolve method")
    # Evolve the stabilizer state with Hadamard gate
    circuit = ccl.CliffordCircuit(n)
    circuit.h(0)
    circuit.h(1)
    circuit.h(2)
    clifford = circuit.compile()
    evolved_state = ccl.StabilizerState.zero_state(n).evolve(clifford)
    # Now the stabilizer should be X1, X2, X3
    expected_generators = [ccl.Pauli.from_pystr('+XII'), ccl.Pauli.from_pystr('+IXI'),ccl.Pauli.from_pystr('+IIX')]
    for i in range(n):
        assert np.array_equal(evolved_state.generators[i].smolin_vec, expected_generators[i].smolin_vec), "Evolved generator does not match expected"
        assert evolved_state.generators[i].exp == expected_generators[i].exp, "Evolved generator exponent does not match expected"
    print("Passed")


    # Now we test a more complicated clifford circuit
    circuit = ccl.CliffordCircuit(n)
    circuit.h(0)
    circuit.cx(0,1)
    circuit.cx(0,2)
    circuit.s(0)
    circuit.s(1)
    circuit.s(2)
    clifford = circuit.compile()
    evolved_state = ccl.StabilizerState.zero_state(n).evolve(clifford)
    expected_generators = [ccl.Pauli.from_pystr("YYY"),ccl.Pauli.from_pystr("ZZI"),ccl.Pauli.from_pystr("ZIZ")]
    for i in range(n):
        assert np.array_equal(evolved_state.generators[i].smolin_vec, expected_generators[i].smolin_vec), "Evolved generator does not match expected"
        assert evolved_state.generators[i].exp == expected_generators[i].exp, "Evolved generator exponent does not match expected"
    print("Passed")
    
    # Test 4: Test get_probability method
    print("Test 4: Test get_probability method")
    
    for _ in range(10):
        probabilities = [evolved_state.get_probability(ccl.int2basestr(i,2,n)) for i in range(2**n)]
        probabilities_exp = evolved_state.to_qiskit().probabilities()
        assert np.array_equal(probabilities, probabilities_exp), f"Outcome probabilities are incorrect"

        evolved_state = evolved_state.evolve(ccl.Clifford.random_clifford(n))

    print("Passed")
    
    # Test 5: Test measure method
    print("Test 5: Test measure method")

    n = 4  # Number of qubits
    num_measurements = 10000  # Number of measurements to perform

    random_state = ccl.StabilizerState.zero_state(n).evolve(ccl.Clifford.random_clifford(n))

    measurement_results = []

    for _ in range(num_measurements):
        bitstring, _ = random_state.measure()
        measurement_results.append(bitstring)

    outcome_counts = Counter(measurement_results)

    expected_probs = {}
    expected_errors = {}
    flag = True
    for i in range(2**n):
        outcome = "".join(map(str,ccl.int2basestr(i, 2, n)))
        true_prob = random_state.get_probability(outcome)
        err = np.sqrt(true_prob * (1-true_prob)/num_measurements)
        empirical_prob = outcome_counts[outcome]/num_measurements

        assert np.abs(true_prob - empirical_prob) <= 3*err, "Measurement probability error not within 99.7 per cent confidence"

    print("Passed")
    
    print("All tests for ccl.StabilizerState class passed.\n")

def test_clifford_circuit():
    """
    Test various functionalities of the ccl.CliffordCircuit class and related classes.
    """
    # Test 1: Initialize and test the zero state
    print("Test 1: Initialize zero state and check probabilities")
    circuit = ccl.CliffordCircuit(num_qubits=1)
    prob_0 = circuit.stab.get_probability('0')
    prob_1 = circuit.stab.get_probability('1')
    assert abs(prob_0 - 1.0) < 1e-6, f"Initial state probability of '0' is not 1, got {prob_0}"
    assert abs(prob_1 - 0.0) < 1e-6, f"Initial state probability of '1' is not 0, got {prob_1}"
    print("Passed")

    # Test 2: Apply H gate and check probabilities
    print("Test 2: Apply H gate and check probabilities")
    circuit.h(0)
    prob_0 = circuit.stab.get_probability('0')
    prob_1 = circuit.stab.get_probability('1')
    assert abs(prob_0 - 0.5) < 1e-6, f"After H, probability of '0' is not 0.5, got {prob_0}"
    assert abs(prob_1 - 0.5) < 1e-6, f"After H, probability of '1' is not 0.5, got {prob_1}"
    print("Passed")

    # Test 3: Apply H gate twice and check probabilities
    print("Test 3: Apply H gate twice and check probabilities")
    circuit.h(0)
    prob_0 = circuit.stab.get_probability('0')
    prob_1 = circuit.stab.get_probability('1')
    assert abs(prob_0 - 1.0) < 1e-6, f"After H^2, probability of '0' is not 1, got {prob_0}"
    assert abs(prob_1 - 0.0) < 1e-6, f"After H^2, probability of '1' is not 0, got {prob_1}"
    print("Passed")

    # Test 4: Create a Bell state and check probabilities
    print("Test 4: Create a Bell state and check probabilities")
    circuit = ccl.CliffordCircuit(num_qubits=2)
    circuit.h(0)
    circuit.cx(0, 1)
    prob_00 = circuit.stab.get_probability('00')
    prob_11 = circuit.stab.get_probability('11')
    prob_01 = circuit.stab.get_probability('01')
    prob_10 = circuit.stab.get_probability('10')
    assert abs(prob_00 - 0.5) < 1e-6, f"Bell state probability of '00' is not 0.5, got {prob_00}"
    assert abs(prob_11 - 0.5) < 1e-6, f"Bell state probability of '11' is not 0.5, got {prob_11}"
    assert abs(prob_01 - 0.0) < 1e-6, f"Bell state probability of '01' is not 0, got {prob_01}"
    assert abs(prob_10 - 0.0) < 1e-6, f"Bell state probability of '10' is not 0, got {prob_10}"
    print("Passed")

    # Test 5: Measure the Bell state and check outcomes
    print("Test 5: Measure the Bell state and check outcomes")
    counts = {'00': 0, '11': 0}
    num_measurements = 1000
    for _ in range(num_measurements):
        bitstring, _ = circuit.stab.measure()
        if bitstring in counts:
            counts[bitstring] += 1
    assert counts['00'] + counts['11'] == num_measurements, "Measurements should be '00' or '11'"
    ratio_00 = counts['00'] / num_measurements
    ratio_11 = counts['11'] / num_measurements
    assert abs(ratio_00 - 0.5) < 0.05, f"Measurement ratio for '00' is not ~0.5, got {ratio_00}"
    assert abs(ratio_11 - 0.5) < 0.05, f"Measurement ratio for '11' is not ~0.5, got {ratio_11}"
    print("Passed")

    # Test 6: Compile the circuit and compare stabilizer states
    print("Test 6: Compile the circuit and compare stabilizer states")
    compiled_clifford = circuit.compile()
    zero_state = ccl.StabilizerState.zero_state(2)
    evolved_state = zero_state.evolve(compiled_clifford)
    for gen1, gen2 in zip(evolved_state.generators, circuit.stab.generators):
        assert np.array_equal(gen1.smolin_vec, gen2.smolin_vec), "Compiled generators do not match"
        assert gen1.exp == gen2.exp, "Compiled generator exponents do not match"
    print("Passed")

    # Test 7: Test Pauli multiplication
    print("Test 7: Test Pauli multiplication")
    pauli1 = ccl.Pauli.from_pystr('X')
    pauli2 = ccl.Pauli.from_pystr('Z')
    pauli3 = pauli1.multiply(pauli2)
    assert str(pauli3) == '-iY', f"Expected -iY, got {pauli3}"
    print("Passed")

    # Test 8: Test Pauli conjugation
    print("Test 8: Test Pauli conjugation")
    circuit = ccl.CliffordCircuit(1)
    circuit.h(0)
    clifford = circuit.compile()
    pauli = ccl.Pauli.from_pystr('Z')
    pauli_conj = pauli.conjugate(clifford)
    assert str(pauli_conj) == '+X', f"Expected +X, got {pauli_conj}"
    print("Passed")

    # Test 9: Apply X gate and check probabilities
    print("Test 9: Apply X gate and check probabilities")
    circuit = ccl.CliffordCircuit(num_qubits=1)
    circuit.x(0)
    prob_0 = circuit.stab.get_probability('0')
    prob_1 = circuit.stab.get_probability('1')
    assert abs(prob_0 - 0.0) < 1e-6, f"After X, probability of '0' is not 0, got {prob_0}"
    assert abs(prob_1 - 1.0) < 1e-6, f"After X, probability of '1' is not 1, got {prob_1}"
    print("Passed")

    # Test 10: Apply Y gate and check probabilities
    print("Test 10: Apply Y gate and check probabilities")
    circuit = ccl.CliffordCircuit(num_qubits=1)
    circuit.y(0)
    prob_0 = circuit.stab.get_probability('0')
    prob_1 = circuit.stab.get_probability('1')
    assert abs(prob_0 - 0.0) < 1e-6, f"After Y, probability of '0' is not 0, got {prob_0}"
    assert abs(prob_1 - 1.0) < 1e-6, f"After Y, probability of '1' is not 1, got {prob_1}"
    print("Passed")

    # Test 11: Apply Z gate and check probabilities
    print("Test 11: Apply Z gate and check probabilities")
    circuit = ccl.CliffordCircuit(num_qubits=1)
    circuit.z(0)
    prob_0 = circuit.stab.get_probability('0')
    prob_1 = circuit.stab.get_probability('1')
    assert abs(prob_0 - 1.0) < 1e-6, f"After Z, probability of '0' is not 1, got {prob_0}"
    assert abs(prob_1 - 0.0) < 1e-6, f"After Z, probability of '1' is not 0, got {prob_1}"
    print("Passed")

    # Test 12: Apply CZ gate and check entanglement
    print("Test 12: Apply CZ gate and check entanglement")
    circuit = ccl.CliffordCircuit(num_qubits=2)
    circuit.h(0)
    circuit.h(1)
    circuit.cz(0, 1)
    prob_00 = circuit.stab.get_probability('00')
    prob_11 = circuit.stab.get_probability('11')
    prob_01 = circuit.stab.get_probability('01')
    prob_10 = circuit.stab.get_probability('10')
    assert abs(prob_00 - 0.25) < 1e-6, f"After CZ, probability of '00' is not 0.25, got {prob_00}"
    assert abs(prob_11 - 0.25) < 1e-6, f"After CZ, probability of '11' is not 0.25, got {prob_11}"
    assert abs(prob_01 - 0.25) < 1e-6, f"After CZ, probability of '01' is not 0.25, got {prob_01}"
    assert abs(prob_10 - 0.25) < 1e-6, f"After CZ, probability of '10' is not 0.25, got {prob_10}"
    print("Passed")

    print("All tests passed!")

if __name__ == "__main__":
    test_clifford_class()
    test_pauli_class()
    test_stabilizer_state_class()
    test_clifford_circuit()

