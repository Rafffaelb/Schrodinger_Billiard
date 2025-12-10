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

## Quick Start

```bash
# Build the project (automatically downloads required Eigen library)
make

# Run a simulation (example: orthogonal conductance vs channel count)
./WignerDyson.exe Orthogonal Channel

# View results in Data_Analysis Jupyter notebooks
```

## Installing Docker (Windows)

To simplify setup, you can use Docker to run this project in a consistent environment:

Option 1 - Guided Installation:
```cmd
install_docker.bat
```

Option 2 - Automated Installation (requires admin privileges):
```cmd
install_docker_enhanced.bat
```

Follow the instructions in the script to download and install Docker Desktop,
then restart your computer and launch Docker Desktop.

## Using Docker

This project can be built and run using Docker for a consistent environment across different systems.

For detailed Docker usage instructions, see [README_DOCKER.md](README_DOCKER.md).

### Building the Docker Image

Make sure you are in the `WignerDyson` directory before building (this is where the Dockerfile is located):

```bash
# Navigate to the WignerDyson directory where the Dockerfile is located
cd WignerDyson

# Build the Docker image
docker build -t schrodinger-billiard .
```

Alternatively, from the parent directory:
```bash
# Build from the parent directory by specifying the path
docker build -t schrodinger-billiard WignerDyson/
```

Or using the provided helper script (Linux/macOS):
```bash
cd WignerDyson
chmod +x docker_helper.sh
./docker_helper.sh build
```

Or using the provided helper script (Windows):
```cmd
cd WignerDyson
docker_helper.bat build
```

Or using docker-compose:

```bash
cd WignerDyson
docker-compose build
```

### Running Simulations with Docker

```bash
# Run the container interactively
docker run -it schrodinger-billiard

# Run a specific simulation
docker run schrodinger-billiard ./WignerDyson.exe Orthogonal Channel

# Using docker-compose to preserve data
docker-compose run schrodinger-billiard ./WignerDyson.exe Orthogonal Channel

# Using the helper script
./docker_helper.sh run
./docker_helper.sh exec "./WignerDyson.exe Orthogonal Channel"
```

### Accessing Results

When using docker-compose or the volume mount options, the Data_Analysis directory is mounted to your host system, so you can access simulation results directly.

### Running Jupyter Notebooks

```bash
docker-compose run schrodinger-billiard jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser --allow-root
```

Then access the notebook at http://localhost:8888 and use the token displayed in the terminal.

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

## Results Visualization

Results are saved as `.dat` files and can be analyzed using the Jupyter notebooks in `WignerDyson/Data_Analysis/`.

## Project Structure

```
.
├── Data_Analysis/          # Jupyter notebooks for data analysis
│   ├── Bell_Parameter/     # Bell parameter analysis results
│   ├── Channel/            # Conductance vs channel count analysis
│   ├── Concurrence/        # Entanglement concurrence analysis
│   ├── Energy/             # Conductance vs energy analysis
│   ├── Gamma/              # Conductance vs gamma parameter analysis
│   └── Various .ipynb files for data visualization and analysis
├── eigen-3.4.0/            # Eigen library (automatically downloaded, Git-ignored)
├── gpu_experiment/         # Experimental GPU-accelerated implementation (Git-ignored)
├── include/                # Header files for all classes
│   ├── Orthogonal.h        # Orthogonal ensemble (GOE) implementation
│   ├── Unitary.h           # Unitary ensemble (GUE) implementation
│   ├── Symplectic.h        # Symplectic ensemble (GSE) implementation
│   ├── Quantum_chaotic_billiard.h  # Main quantum billiard class
│   └── WignerDyson_AbstractClass_h/  # Abstract base classes
├── src/                    # Source files
│   ├── Orthogonal.cpp      # Orthogonal ensemble implementation
│   ├── Unitary.cpp         # Unitary ensemble implementation
│   ├── Symplectic.cpp      # Symplectic ensemble implementation
│   ├── Quantum_chaotic_billiard.cpp  # Main quantum billiard implementation
│   └── WignerDyson_AbstractClass_cpp/  # Abstract class implementations
├── third_party/            # Third-party installation scripts (Git-ignored)
├── Dockerfile              # Docker configuration for containerization
├── docker-compose.yml      # Docker Compose configuration for multi-container setups
├── .dockerignore           # Files to ignore when building Docker images (Git-ignored)
├── Makefile               # Build configuration
├── main.cpp               # Main program entry point
└── Various .bat files     # Installation scripts for dependencies (Git-ignored)
```

**Note on Git-ignored Directories**: Several directories in this project are not tracked by Git:
- `eigen-3.4.0/`: Third-party Eigen library, automatically downloaded during build
- `gpu_experiment/`: Experimental GPU-accelerated code that is separate from the main project
- `third_party/`: Scripts for installing additional dependencies
- `.dockerignore`: Docker ignore file (not needed in version control)
- `.bat files`: Windows batch scripts for dependency installation

These directories are excluded via `.gitignore` because they contain either third-party code or platform-specific installation scripts. They are essential for building and running the project locally but are not part of the core source code.
