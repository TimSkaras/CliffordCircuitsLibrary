# Clifford Circuits Library
Welcome to the Clifford Circuits Library! This library provides classes and methods for simulating Clifford circuits, stabilizer states, and Pauli operators using symplectic matrices. It's designed for quantum computing simulations involving Clifford operations, stabilizer states, and Pauli measurements.

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Getting Started](#getting-started)
- [Examples](#examples)
- [License](#license)

## Introduction

This library provides a framework for simulating quantum circuits that are composed of Clifford gates. It includes classes to represent Pauli operators, Clifford operations, and stabilizer states. The stabilizer formalism allows efficient simulation of a subset of quantum computations that include state preparation, Clifford gates, and Pauli measurements.

### Key Components:

- Pauli Class (`ccl.Pauli`): Represents Pauli operators and supports operations like multiplication and commutation checks.
- Clifford Class (`ccl.Clifford`): Represents Clifford operations using symplectic matrices.
- StabilizerState Class (`ccl.StabilizerState`): Represents stabilizer states and supports evolution under Clifford operations and measurements.
- CliffordCircuit Class (`ccl.CliffordCircuit`): Represents a circuit composed of Clifford gates and provides methods to build and compile Clifford operations.

## Installation

To use the Clifford Circuits Library, you need to include the `ccl.py` file in your project directory or ensure that Python can find it by adjusting the `PYTHONPATH` environment variable.

### Option 1: Place `ccl.py` in Your Project Directory

Copy the `ccl.py` file into the same directory as your Python script. You can then import the ccl module directly.

```python
import ccl
```

### Option 2: Modify `PYTHONPATH`

If you prefer to keep `ccl.py` in a separate directory, you can modify the `PYTHONPATH` environment variable to include the directory containing `ccl.py`.

#### On Windows (Command Prompt):
```cmd
set PYTHONPATH=%PYTHONPATH%;C:\path\to\CliffordCircuitsLibrary
```

#### On Unix/Linux or Windows PowerShell:
```bash
export PYTHONPATH="${PYTHONPATH}:/path/to/CliffordCircuitsLibrary"
```

Then you can import ccl as you would a normal python package

```python
import ccl
```

## Getting Started

### Creating Pauli Objects

You can create Pauli operators using the `Pauli` class. Pauli operators are represented using an array of bits in smolin format. You can create `Pauli` objects by specifying its smolin vector or python string.

#### Smolin Vector
 The most direct method for creating a `Pauli` object is to initialize one with its smolin vector format. 

 | Smolin Vector | Pauli |
|-------------|-------------|
| [0,0]      | I       |
| [0,1]   | Z        |
| [1,0]   | X        |
| [1,1]   | Y        |


```python
import ccl

pauli_I = ccl.Pauli([0,0])
pauli_XZ = ccl.Pauli([1,0,0,1])
```

#### Python String

For something less arcane, you can create a Pauli object by specifying its python string directly.

```python
import ccl

pauli_YIZY = ccl.Pauli.from_pystr("YIZY")
pauli_IIZ = ccl.Pauli.from_pystr("IIZ")
```

### Creating Clifford Objects

Clifford operations can be created using the `Clifford` class. You can create a Clifford object by specifying its index in the symplectic group or by providing a symplectic matrix directly.

```python
import ccl

# Create a Clifford object (identity operation in this case)
n = 1  # Number of qubits
cliff_idx = 3  # Index in the symplectic group
sign_idx = 3
clifford = ccl.Clifford(cliff_idx, sign_idx, n)

print("Clifford representation:")
print(clifford)
```

The clifford index specifies the symplectic matrix of the clifford unitary. This matrix specifies how each Pauli gets mapped when conjugated by the Clifford. It does not specify the sign that each Pauli would pick up, however. The sign is specified in the sign index, which when converted to a bit string of length 2n, tells you whether each Pauli gets a minus sign or not. It is easier to understand with an example.

Suppose we are trying to create a Clifford object to represent the circuit XH, i.e., the circuit applies a hadamard gate and then an X gate. What are the symplectic matrix and sign array for this Clifford? The symplectic matrix is determined by the action of the clifford on X and Z. We have the following relationship

 | P | $UPU^\dagger$ |
|-------------|-------------|
| X      | -Z       |
| Z      |  X     |

which in smolin terms means X (the vector [1,0]) gets mapped to -Z (the vector [0,1] and a sign of -1) and Z (the vector [0,1]) gets mapped to X (the vector [1,0]). In Smolin convention, the first column tells us how $X_1$ transforms, the second column how $Z_1$ transforms, the third column how $X_2$ transforms, and so on. The matrix specifying this linear transformation is just

$$
\begin{bmatrix}
0 & 1 \\
1 & 0
\end{bmatrix}.
$$

You can check this is correct by taking a smolin vector as a column and multiplying this matrix on the left. What about the sign though? The symplectic matrix does not track the signs, so we have to keep track of those in a separate array. In this case, X need a sign flip, so our entire Clifford can be represented

$$
\begin{array}{cc}
0 & 1 \\
1 & 0 \\
\hline
1 & 0
\end{array}
$$


### Creating Stabilizer States

## Examples

### Evolving with Random Cliffords

### Measuring a Stabilizer State

### Getting Probabilities of a Stabilizer State

### Building a Clifford Circuit

## License

This project is licensed under the MIT License.
