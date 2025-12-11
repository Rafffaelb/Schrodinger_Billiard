#!/bin/bash

# Build script for Schrodinger Billiard project on Linux

echo "Building Schrodinger Billiard project..."

# Compile the project (the Makefile will automatically download Eigen if needed)
echo "Compiling project..."
make -f Makefile.linux

if [ $? -eq 0 ]; then
    echo "Build successful!"
    echo "You can now run the program with:"
    echo "  ./WignerDyson [options]"
    echo ""
    echo "To see available options, run:"
    echo "  ./WignerDyson --help"
else
    echo "Build failed!"
    exit 1
fi