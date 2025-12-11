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
- wget or curl (for automatic Eigen library download)

Note: The Eigen library is automatically downloaded during the build process, so no manual installation is required.

## Building the Project

### Windows Build

To build the project on Windows, simply run:

```bash
make
```

This will:
1. Automatically download and extract the Eigen library if not already present
2. Compile all source files
3. Generate the executable `WignerDyson.exe`

### Linux Build

To build the project on Linux, use the Linux-specific Makefile:

```bash
make -f Makefile.linux
```

This will:
1. Automatically download and extract the Eigen library if not already present
2. Compile all source files
3. Generate the executable `WignerDyson`

You can also use the provided build script:

```bash
./build.sh
```

To manually download Eigen:

```bash
make -f Makefile.linux eigen-download
```

To clean object files and executables:

```bash
make -f Makefile.linux clean
```

To perform a full clean (including Eigen library):

```bash
make -f Makefile.linux clean-all
```

## Running Simulations

The simulator supports various symmetry classes and simulation types. Run it with the following syntax:

```bash
./WignerDyson.exe <symmetry_class> <simulation_type> [other_parameters]
```

On Linux systems:

```bash
./WignerDyson <symmetry_class> <simulation_type> [other_parameters]
```

### Symmetry Classes

- `Orthogonal` - Gaussian Orthogonal Ensemble (GOE)
- `Unitary` - Gaussian Unitary Ensemble (GUE)
- `Symplectic` - Gaussian Symplectic Ensemble (GSE)

### Simulation Types

#### Conductance Analysis
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
- `Talphabeta` - T_αβ scattering matrix elements

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
- Scattering data: `Talphabeta.dat`

## Data Analysis

The [Data_Analysis](file:///c:/Users/rafam/Documents/Physics_content/Schrodinger_Billiard/WignerDyson/Data_Analysis) directory contains Jupyter Notebooks for analyzing the simulation results:
- `Bell_Parameter.ipynb`
- `Concurrence_Gamma.ipynb`
- `Conductance_Energy.ipynb`
- `Conductance_ShotNoise.ipynb`
- `Talphabeta.ipynb`

## Project Structure

```
.
├── Data_Analysis/          # Jupyter notebooks for data analysis
├── include/                # Header files
│   ├── WignerDyson_AbstractClass_h/
│   │   ├── WignerDyson.h   # Abstract base class
│   │   └── Run_Simulation_*.h # Simulation modules
│   ├── Orthogonal.h        # Orthogonal symmetry implementation
│   ├── Unitary.h           # Unitary symmetry implementation
│   ├── Symplectic.h        # Symplectic symmetry implementation
│   └── Quantum_chaotic_billiard.h
├── src/                    # Source files
│   ├── WignerDyson_AbstractClass_cpp/
│   │   ├── WignerDyson.cpp # Base class implementation
│   │   └── Run_Simulation_*.cpp # Simulation module implementations
│   ├── Orthogonal.cpp      # Orthogonal symmetry implementation
│   ├── Unitary.cpp         # Unitary symmetry implementation
│   ├── Symplectic.cpp      # Symplectic symmetry implementation
│   └── Quantum_chaotic_billiard.cpp
├── Makefile                # Build configuration
├── download_eigen.bat      # Eigen download script for Windows
├── main.cpp                # Main entry point
└── README.md               # This file
```

## Technical Details

- Uses Eigen library for linear algebra operations (automatically downloaded)
- Implements OpenMP parallelization for ensemble averaging
- Supports FMA (Fused Multiply-Add) instruction set optimization
- Cross-platform build system with automatic dependency management

## Troubleshooting

### Eigen download issues
If the automatic Eigen download fails, you can manually download it:
1. Download Eigen 3.4.0 from https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2
2. Extract it in the project root directory
3. Rename the folder to `eigen-3.4.0`

### Compilation errors
Ensure you're using a C++11 compatible compiler with OpenMP support.

## License

This project is intended for scientific research purposes.