#!/bin/bash

# Run script for Schrodinger Billiard project

# Default parameters
SYMMETRY="orthogonal"  # orthogonal, unitary, symplectic
CHANNELS=2             # Number of channels
GAMMA=1.0              # Barrier transparency
SAMPLES=1000           # Number of samples

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--symmetry)
            SYMMETRY="$2"
            shift 2
            ;;
        -c|--channels)
            CHANNELS="$2"
            shift 2
            ;;
        -g|--gamma)
            GAMMA="$2"
            shift 2
            ;;
        -n|--samples)
            SAMPLES="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo "Run the Schrodinger Billiard simulation"
            echo ""
            echo "Options:"
            echo "  -s, --symmetry SYM    Symmetry class: orthogonal, unitary, symplectic (default: orthogonal)"
            echo "  -c, --channels NUM    Number of channels (default: 2)"
            echo "  -g, --gamma VAL       Barrier transparency (default: 1.0)"
            echo "  -n, --samples NUM     Number of samples (default: 1000)"
            echo "  -h, --help            Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0                           # Run with default parameters"
            echo "  $0 -s unitary -c 4 -g 0.5    # Run unitary symmetry with 4 channels and gamma=0.5"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

# Check if the executable exists
if [ ! -f "WignerDyson" ]; then
    echo "Error: WignerDyson executable not found!"
    echo "Please build the project first with:"
    echo "  make -f Makefile.linux"
    exit 1
fi

# Run the simulation
echo "Running Schrodinger Billiard simulation..."
echo "Symmetry: $SYMMETRY"
echo "Channels: $CHANNELS"
echo "Gamma: $GAMMA"
echo "Samples: $SAMPLES"
echo ""

./WignerDyson --symmetry $SYMMETRY --channels $CHANNELS --gamma $GAMMA --samples $SAMPLES