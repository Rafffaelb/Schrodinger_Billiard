# Schrödinger's Billiard: Quantum Chaos Simulator

This project simulates quantum billiards - nanoscale devices where electrons behave according to quantum mechanics. Using Random Matrix Theory, we study how quantum chaos affects electron transport properties.

## What are Quantum Billiards?

Quantum billiards (also known as quantum dots) are tiny devices (~nanometers) where electrons are confined in a small region and can move between leads through a chaotic cavity. Think of it like a microscopic pool table where electrons bounce around chaotically.

These systems exhibit quantum behavior:
- Electrons can only occupy discrete energy levels
- Wave-particle duality causes interference effects
- Quantum entanglement can occur between electrons exiting different leads

## Physics Concepts Studied

### Symmetry Classes
We simulate three fundamental symmetry classes that describe different physical conditions:
- **Orthogonal (GOE)**: Time-reversal symmetry preserved, no magnetic field
- **Unitary (GUE)**: Time-reversal symmetry broken, magnetic field applied
- **Symplectic (GSE)**: Time-reversal symmetry preserved but spin-orbit coupling present

### Physical Properties Analyzed
- **Conductance**: How easily electrons flow through the billiard
- **Shot Noise**: Fluctuations in electric current due to discrete electron charges
- **Entanglement**: Quantum correlation between electrons exiting different leads
- **Bell Inequality Violations**: Demonstrating non-classical correlations

## How the Code Works

The simulation is written in C++ and uses Random Matrix Theory to model quantum chaos:

1. **Hamiltonian Generation**: Random matrices representing the quantum system's energy levels
2. **Scattering Matrix Calculation**: Using the Mahaux-Weidenmüller formula to compute how electrons scatter
3. **Physical Property Computation**: Applying formulas like Landauer-Büttiker for conductance
4. **Statistical Analysis**: Computing averages, variances, and distributions over many realizations

## Quick Start

```bash
# Build the project (automatically downloads required Eigen library)
make

# Run a simulation (example: orthogonal conductance vs channel count)
./WignerDyson.exe Orthogonal Channel

# View results in Data_Analysis Jupyter notebooks
```

## Results Visualization

Results are saved as `.dat` files and can be analyzed using the Jupyter notebooks in `WignerDyson/Data_Analysis/`.

See the detailed [technical README](WignerDyson/README.md) for comprehensive build and usage instructions.