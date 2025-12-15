# Docker Setup for Schrödinger Billiard Simulations

This project includes Docker configuration to easily run the quantum chaotic billiard simulations in a containerized environment.

## Prerequisites

- Docker installed on your system
- Docker Compose (included with Docker Desktop or install separately)

## Building the Docker Image

To build the Docker image, run:

```bash
docker-compose build
```

## Running Simulations

### Using Docker Compose (Recommended)

To run a simulation with specific parameters:

```bash
# Example: Run Orthogonal symmetry class with Channel simulation
docker-compose run --rm schrodinger-billiard ./WignerDyson Orthogonal Channel

# Example: Run Unitary symmetry class with Gamma simulation
docker-compose run --rm schrodinger-billiard ./WignerDyson Unitary Gamma

# Example: Run Symplectic symmetry class with Bell_Parameter_Ress simulation
docker-compose run --rm schrodinger-billiard ./WignerDyson Symplectic Bell_Parameter_Ress
```

The results will be saved in the `Data_Analysis` folder on your host machine.

### Available Parameters

1. **Symmetry Classes:**
   - `Orthogonal` - Gaussian Orthogonal Ensemble (GOE) - Time-reversal symmetric
   - `Unitary` - Gaussian Unitary Ensemble (GUE) - Time-reversal broken
   - `Symplectic` - Gaussian Symplectic Ensemble (GSE) - Time-reversal symmetric with spin-orbit coupling

2. **Simulation Types:**
   - `Channel` - Conductance vs channel count
   - `Gamma` - Conductance vs gamma parameter
   - `Concurrence` - Entanglement concurrence vs gamma
   - `Bell_Parameter_Ress` - Bell parameter vs Ress
   - `Bell_Parameter_Gamma` - Bell parameter vs gamma
   - `Bell_Parameter_Fixed_Base` - Bell parameter with fixed basis
   - `Energy` - Conductance vs energy
   - `Energy_Gamma` - Conductance vs energy with gamma variation

## Running Jupyter Notebooks

To analyze the results using Jupyter notebooks:

```bash
docker-compose up jupyter
```

Then open your browser at http://localhost:8888 to access the Jupyter interface.

## Data Persistence

Simulation results are stored in the `Data_Analysis` directory, which is mounted as a volume in the container. This ensures that your data persists even after the container is stopped or removed.

## Customization

You can customize the UID/GID used in the container by setting environment variables:

```bash
UID=1001 GID=1001 docker-compose up
```

By default, UID and GID are set to 1000.
