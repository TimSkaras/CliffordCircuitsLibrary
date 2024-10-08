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

You can create Pauli operators using the `Pauli` class. Pauli operators are represented using vectors of bits in smolin format. You can create `Pauli` objects by specifying its smolin vector or python string.

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

# Create Pauli operators from string representation
pauli_x = ccl.Pauli.from_pystr('X')
pauli_y = ccl.Pauli.from_pystr('Y')
pauli_z = ccl.Pauli.from_pystr('Z')
pauli_i = ccl.Pauli.from_pystr('I')

print(f"Pauli X: {pauli_x}")
print(f"Pauli Y: {pauli_y}")
print(f"Pauli Z: {pauli_z}")
print(f"Pauli I: {pauli_i}")
```


### Creating Clifford Objects

### Creating Stabilizer States

## Examples

### Evolving with Random Cliffords

### Measuring a Stabilizer State

### Getting Probabilities of a Stabilizer State

### Building a Clifford Circuit

## License

This project is licensed under the MIT License.
