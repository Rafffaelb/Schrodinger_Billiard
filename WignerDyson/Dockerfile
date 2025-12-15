# Use Ubuntu 22.04 as the base image for better support and security updates
FROM ubuntu:22.04

# Set noninteractive installation to avoid prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages in a single layer for efficiency
RUN apt-get update && apt-get install -y \
    build-essential \
    g++ \
    make \
    wget \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages for data analysis
RUN pip3 install numpy matplotlib jupyter pandas

# Create a non-root user with dynamic UID/GID that we'll map to the host user
ARG USER_ID=1000
ARG GROUP_ID=1000
RUN addgroup --gid ${GROUP_ID} appuser &&\
    adduser --uid ${USER_ID} --gid ${GROUP_ID} --disabled-password --gecos '' appuser

# Set working directory
WORKDIR /app

# Download Eigen first to ensure it's available for compilation
RUN echo "Downloading Eigen library..." && \
    wget -O eigen.tar.bz2 "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2" && \
    tar -xjf eigen.tar.bz2 && \
    rm eigen.tar.bz2 && \
    chown -R ${USER_ID}:${GROUP_ID} eigen-3.4.0

# Copy source files for better layer caching
COPY Makefile.linux ./
COPY src/ src/
COPY include/ include/
COPY main.cpp ./

# Fix the include paths in source files to match the Makefile's -I flags
RUN find src/ -name "*.cpp" -exec sed -i 's|#include <eigen3/Eigen/Dense>|#include <Eigen/Dense>|g' {} \;

# Build the project using the Linux Makefile
RUN make -f Makefile.linux clean && make -f Makefile.linux && \
    chown ${USER_ID}:${GROUP_ID} WignerDyson

# Copy remaining files (data analysis notebooks, etc.)
COPY . .

# Change ownership of the app directory to the appuser
RUN chown -R ${USER_ID}:${GROUP_ID} /app

# Switch to the non-root user
USER ${USER_ID}:${GROUP_ID}

# Set default OpenMP threads to use all available cores
ENV OMP_NUM_THREADS=*

# Expose port for Jupyter notebook (if needed)
EXPOSE 8888

# Default command - show usage
CMD ["./WignerDyson", "--help"]