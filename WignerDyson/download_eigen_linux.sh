#!/bin/bash

# Script to download and extract Eigen library on Linux

EIGEN_VERSION="3.4.0"
EIGEN_DIR="eigen-3.4.0"
DOWNLOAD_URL="https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.bz2"

# Check if Eigen already exists
if [ -d "$EIGEN_DIR" ]; then
    echo "Eigen library already exists in $EIGEN_DIR"
    exit 0
fi

echo "Downloading Eigen library version $EIGEN_VERSION..."

# Download Eigen
if command -v wget >/dev/null 2>&1; then
    wget -O eigen.tar.bz2 "$DOWNLOAD_URL"
elif command -v curl >/dev/null 2>&1; then
    curl -L -o eigen.tar.bz2 "$DOWNLOAD_URL"
else
    echo "Error: Neither wget nor curl found. Please install one of them."
    exit 1
fi

# Check if download was successful
if [ ! -f "eigen.tar.bz2" ]; then
    echo "Error: Failed to download Eigen library"
    exit 1
fi

echo "Extracting Eigen library..."

# Extract the archive
if ! tar -xjf eigen.tar.bz2; then
    echo "Error: Failed to extract Eigen library"
    rm -f eigen.tar.bz2
    exit 1
fi

# Clean up
rm -f eigen.tar.bz2

echo "Eigen library successfully downloaded and extracted to $EIGEN_DIR"
echo "You can now build the project with:"
echo "  make -f Makefile.linux"