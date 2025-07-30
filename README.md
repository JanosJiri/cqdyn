# Coefficient-based Quantum Dynamics (cQDyn)
`cQDyn` is a Python package for simulating quantum dynamics using state coefficients and Hamiltonian matrices. The code allows inclusion of external field via interaction Hamiltonian. 

## Theory

The wave function $\Psi$ is considered as an expansion of eigenstates $\phi_i$

$$\Psi(t)=\Sigma_i c_i(t) \phi_i$$

where $c_i$ are the time-dependent expansion coefficients. The Hamiltonian is considered in the following form:

$$\hat{H} = \hat{T} + \hat{V} + \hat{V}\_\mathrm{int}\cdot\varepsilon (t) = \hat{H}\_0 + \hat{V}_\mathrm{int}\cdot\varepsilon (t)$$

$\hat{H}\_0$ is the Hamiltonian without interaction, $\hat{V}_\mathrm{int}$ is the interaction strength and $\varepsilon$ is the time-dependent driving field of the interaction, e.g., electric field of a laser pulse. Inserting the ansatz into the time-dependent Schrödinger equations leads to

$$i\hbar \frac{\mathrm{d}}{\mathrm{d}t}\vec{c}(t) = \mathbb{H}\vec{c}$$

with the Hamiltonian matrix $\mathbb{H}\_{ij} = \langle\phi_i|\hat{H}|\phi_j\rangle$ written in components as 

$$\mathbb{H} = \mathbb{H}\_0+\mathbb{V}_\mathrm{int}\cdot\varepsilon (t)$$

Hence, the time-independent Hamiltonian $\mathbb{H}\_{0}$ and the interaction Hamiltonian $\mathbb{V}_\mathrm{int}$ multiplied by a time-dependent field (interaction) $\varepsilon$ are necessary to specify the system.

## Installation
To install the package, clone the repository and run the following command in the terminal in the root directory of the project:

```bash
pip install .
```

After that, you can run the code from the command line as

```bash
cqdyn
```

or import it in your Python scripts.

## Input

The input is given in the JSON format in `input.json` with the following structure
```
{
  "total_time" : 10,
  "dt" : 0.001,
  "print_time" : 0.5,
  "coefficients": ["1+0j", "0+0j"],
  "H_0" : [[2.0000, 0.2000], [0.2000, 1.0000]],
  "V_int" : [[0.0000, 0.2000],[ 0.2000, 0.0000]],
  "field" : "sin(t)"
}
```
where `H_0` is the Hamiltonian without interaction, `V_int` is the interaction Hamiltonian and `field` is the interaction field. Note that `field` is a string with prescription of the time-dependent field where the time is denoted as 't'. `coefficients` hold the state coefficients at the beginning of the dynamics. 

## Output

Python arrays in the binary format are saved in the `cqdyn.npz` file. The file contains the following arrays:
- `time`: time points of the simulation
- `coefficients`: complex coefficients of the wave function at each time point
- `energy`: energy at each time point
Loading the file can be done in Python as
```python
import numpy as np
data = np.load('cqdyn.npz')
time = data['time']
coefficients = data['coefficients']
energy = data['energy']
```

Energies and coefficcients are also save in a human-readable text format in `energy.txt` and `coefficients.txt` files.

The cQDyn output from the terminal is saved to `cqdyn.out`.