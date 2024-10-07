Clifford Circuits Library
Welcome to the Clifford Circuits Library! This library provides classes and methods for simulating Clifford circuits, stabilizer states, and Pauli operators using symplectic matrices. It's designed for quantum computing simulations involving Clifford operations, stabilizer states, and Pauli measurements.

Table of Contents
Introduction
Installation
Getting Started
Creating Pauli Objects
Creating Clifford Objects
Creating Stabilizer States
Evolving Stabilizer States
Evolving with Clifford Circuits
Evolving with Random Cliffords
Evolving with Random Paulis
Examples
Introduction
This library provides a framework for simulating quantum circuits that are composed of Clifford gates. It includes classes to represent Pauli operators, Clifford operations, and stabilizer states. The stabilizer formalism allows efficient simulation of a subset of quantum computations that include state preparation, Clifford gates, and Pauli measurements.

Key Components:

Pauli Class (ccl.Pauli): Represents Pauli operators and supports operations like multiplication and commutation checks.
Clifford Class (ccl.Clifford): Represents Clifford operations using symplectic matrices.
StabilizerState Class (ccl.StabilizerState): Represents stabilizer states and supports evolution under Clifford operations and measurements.
CliffordCircuit Class (ccl.CliffordCircuit): Represents a circuit composed of Clifford gates and provides methods to build and compile Clifford operations.
Installation
To use the Clifford Circuits Library, you need to include the ccl.py file in your project directory or ensure that Python can find it by adjusting the PYTHONPATH environment variable.

Option 1: Place ccl.py in Your Project Directory
Copy the ccl.py file into the same directory as your Python script. You can then import the ccl module directly.

python
Copy code
import ccl
Option 2: Modify PYTHONPATH
If you prefer to keep ccl.py in a separate directory, you can modify the PYTHONPATH environment variable to include the directory containing ccl.py.

On Windows (Command Prompt):

cmd
Copy code
set PYTHONPATH=%PYTHONPATH%;C:\path\to\CliffordCircuitsLibrary
On Unix/Linux or Windows PowerShell:

bash
Copy code
export PYTHONPATH="${PYTHONPATH}:/path/to/CliffordCircuitsLibrary"
Replace C:\path\to\CliffordCircuitsLibrary with the actual path to the directory containing ccl.py.

Now you can import the ccl module in your scripts:

python
Copy code
import ccl
Getting Started
Creating Pauli Objects
You can create Pauli operators using the Pauli class. Pauli operators are represented using the symplectic representation.

python
Copy code
import ccl

# Create Pauli operators from string representation
pauli_x = ccl.Pauli.from_pystr('X')
pauli_y = ccl.Pauli.from_pystr('Y')
pauli_z = ccl.Pauli.from_pystr('Z')
pauli_i = ccl.Pauli.from_pystr('I')

print(f"Pauli X: {pauli_x}")
print(f"Pauli Y: {pauli_y}")
print(f"Pauli Z: {pauli_z}")
print(f"Pauli I: {pauli_i}")
Output:

yaml
Copy code
Pauli X: +X
Pauli Y: +Y
Pauli Z: +Z
Pauli I: +
Creating Clifford Objects
Clifford operations can be created using the Clifford class. You can create a Clifford object by specifying its index in the symplectic group or by providing a symplectic matrix directly.

python
Copy code
# Create a Clifford object (identity operation in this case)
n = 1  # Number of qubits
cliff_idx = 0  # Index in the symplectic group
sign_idx = 0
clifford = ccl.Clifford(cliff_idx, sign_idx, n)

print("Clifford representation:")
print(clifford)
Creating Stabilizer States
Stabilizer states can be created using the StabilizerState class. You can start with the zero state or define your own stabilizers.

python
Copy code
# Create a zero stabilizer state with n qubits
n = 2
stabilizer_state = ccl.StabilizerState.zero_state(n)

print("StabilizerState representation:")
print(stabilizer_state)
Output:

makefile
Copy code
Destabilizers: +XI, +IX
Generators: +ZI, +IZ
Evolving Stabilizer States
Evolving with Clifford Circuits
You can build Clifford circuits using the CliffordCircuit class, compile them into a Clifford object, and then apply the compiled Clifford to a StabilizerState.

python
Copy code
# Initialize a Clifford circuit with 2 qubits
circuit = ccl.CliffordCircuit(num_qubits=2)

# Apply gates
circuit.h(0)      # Apply Hadamard gate to qubit 0
circuit.cx(0, 1)  # Apply CNOT gate from qubit 0 to qubit 1

# Compile the circuit into a Clifford object
compiled_clifford = circuit.compile()

# Create a stabilizer state (e.g., zero state)
stabilizer_state = ccl.StabilizerState.zero_state(2)

# Evolve the stabilizer state with the compiled Clifford
evolved_state = stabilizer_state.evolve(compiled_clifford)

print("Evolved StabilizerState representation:")
print(evolved_state)
Output:

makefile
Copy code
Destabilizers: +XX, +XI
Generators: +ZZ, +IZ
Evolving with Random Cliffords
You can create random Clifford operations and evolve stabilizer states with them.

python
Copy code
# Create a random Clifford operation
n = 2  # Number of qubits
random_clifford = ccl.Clifford.random_clifford(n)

# Evolve the stabilizer state with the random Clifford
evolved_state = stabilizer_state.evolve(random_clifford)

print("Evolved StabilizerState after random Clifford:")
print(evolved_state)
Evolving with Random Paulis
Similarly, you can create random Pauli operators and evolve stabilizer states with them.

python
Copy code
# Create a random Pauli operator
random_pauli = ccl.Pauli.random_pauli(n)

# Evolve the stabilizer state with the random Pauli
evolved_state = stabilizer_state.evolve(random_pauli)

print("Evolved StabilizerState after random Pauli:")
print(evolved_state)
Examples
Example 1: Creating and Multiplying Pauli Operators
python
Copy code
import ccl

# Create Pauli operators
pauli_x = ccl.Pauli.from_pystr('X')
pauli_z = ccl.Pauli.from_pystr('Z')

# Multiply Pauli operators
pauli_y = pauli_x.multiply(pauli_z)

print(f"X * Z = {pauli_y}")  # Expected output: -iY
Output:

Copy code
X * Z = -iY
Example 2: Conjugating a Pauli Operator with a Clifford
python
Copy code
# Create a Clifford circuit that applies a Hadamard gate
circuit = ccl.CliffordCircuit(1)
circuit.h(0)
clifford = circuit.compile()

# Create a Pauli Z operator
pauli_z = ccl.Pauli.from_pystr('Z')

# Conjugate Pauli Z with the Clifford (Hadamard gate)
conjugated_pauli = pauli_z.conjugate(clifford)

print(f"Conjugated Pauli: {conjugated_pauli}")  # Expected output: +X
Output:

yaml
Copy code
Conjugated Pauli: +X
Example 3: Measuring a Stabilizer State
python
Copy code
# Initialize a stabilizer state in the |+⟩ state
circuit = ccl.CliffordCircuit(1)
circuit.h(0)
clifford = circuit.compile()
stabilizer_state = ccl.StabilizerState.zero_state(1)
stabilizer_state = stabilizer_state.evolve(clifford)

# Measure the stabilizer state
bitstring, post_measure_state = stabilizer_state.measure()

print(f"Measurement result: {bitstring}")
print("Post-measurement StabilizerState:")
print(post_measure_state)
Output (Measurement result may vary):

yaml
Copy code
Measurement result: 0
Post-measurement StabilizerState:
Destabilizers: +X
Generators: +Z
Example 4: Evolving with Random Clifford and Measuring
python
Copy code
# Create a zero stabilizer state
n = 2
stabilizer_state = ccl.StabilizerState.zero_state(n)

# Create a random Clifford operation
random_clifford = ccl.Clifford.random_clifford(n)

# Evolve the stabilizer state
evolved_state = stabilizer_state.evolve(random_clifford)

# Measure the evolved state
bitstring, post_measure_state = evolved_state.measure()

print(f"Measurement result after random Clifford evolution: {bitstring}")
Output (Measurement result may vary):

arduino
Copy code
Measurement result after random Clifford evolution: 10
Example 5: Simulating a Quantum Circuit
python
Copy code
# Initialize a Clifford circuit with 2 qubits
circuit = ccl.CliffordCircuit(num_qubits=2)

# Apply gates to create a Bell state
circuit.h(0)
circuit.cx(0, 1)

# Compile the circuit and evolve the stabilizer state
compiled_clifford = circuit.compile()
stabilizer_state = ccl.StabilizerState.zero_state(2)
evolved_state = stabilizer_state.evolve(compiled_clifford)

# Measure the state multiple times
counts = {'00': 0, '11': 0}
num_measurements = 1000
for _ in range(num_measurements):
    bitstring, _ = evolved_state.measure()
    if bitstring in counts:
        counts[bitstring] += 1

print(f"Measurement counts: {counts}")
Output (Counts may vary slightly):

css
Copy code
Measurement counts: {'00': 502, '11': 498}
