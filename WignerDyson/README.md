# WignerDyson Quantum Chaos Simulator

A C++ implementation for simulating quantum chaotic systems using Random Matrix Theory (RMT), focusing on the Wigner-Dyson classification of quantum systems with orthogonal, unitary, and symplectic symmetries.

## Overview

This project analyzes statistical physics behaviors in quantum chaotic systems, with particular focus on:
- Level statistics analysis (Wigner surmise verification)
- Quantum entanglement measures (Concurrence, Bell parameters)
- Transport properties of quantum dots (conductance, shot noise)
- Energy dependence and channel coupling effects

## Prerequisites

- C++11 compatible compiler (g++)
- GNU Make
- Internet connection (for automatic Eigen library download)

Note: The Eigen library is automatically downloaded during the build process, so no manual installation is required.

## Building the Project

To build the project, simply run:

```bash
make
```

This will:
1. Automatically download and extract the Eigen library if not already present
2. Compile all source files

## Simulation Types

### Quantum Chaos Regimes
- `Orthogonal` - Time-reversal symmetric systems (β=1)
- `Unitary` - Systems with broken time-reversal symmetry (β=2)
- `Symplectic` - Systems with spin-orbit coupling (β=4)

### Physical Observables

#### Transport Properties
- `Channel` - Conductance vs channel count
- `Gamma` - Conductance vs gamma parameter
- `Energy` - Conductance vs energy
- `Energy_Gamma` - Conductance vs energy with gamma parameter variation

#### Quantum Entanglement
- `Concurrence` - Entanglement concurrence vs gamma parameter
- `Bell_Parameter_Ress` - Bell parameter vs Ress parameter
- `Bell_Parameter_Gamma` - Bell parameter vs gamma parameter
- `Bell_Parameter_Fixed_Base` - Bell parameter with fixed basis

#### Scattering Properties

### Example Commands

```bash
# Run orthogonal conductance vs channel count simulation
./WignerDyson.exe Orthogonal Channel

# Run unitary Bell parameter vs gamma simulation
./WignerDyson.exe Unitary Bell_Parameter_Gamma

# Run symplectic conductance vs channel count simulation
./WignerDyson.exe Symplectic Channel
```

## Output Files

Simulation results are written to `.dat` files in the current directory:
- Conductance data: `Conductance_Channels.dat`, `Conductance_Gamma.dat`, etc.
- Entanglement data: `Concurrence_Gamma.dat`, `Bell_Parameter_Ress.dat`, etc.

## Data Analysis

The [Data_Analysis](file:///c:/Users/rafam/Documents/Physics_content/Schrodinger_Billiard/WignerDyson/Data_Analysis) directory contains Jupyter Notebooks for analyzing the simulation results:
- `Bell_Parameter.ipynb`
- `Concurrence_Gamma.ipynb`
- `Conductance_Energy.ipynb`
- `Conductance_ShotNoise.ipynb`

## Project Structure

```
.
├── Data_Analysis/          # Jupyter notebooks for data analysis
├── include/                # Header files
```