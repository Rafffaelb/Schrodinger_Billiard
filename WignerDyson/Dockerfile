# Use Ubuntu as the base image
FROM ubuntu:20.04

# Set noninteractive installation to avoid prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages
RUN apt-get update && apt-get install -y \
    build-essential \
    g++ \
    make \
    wget \
    unzip \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages for data analysis
RUN pip3 install numpy matplotlib jupyter pandas

# Set working directory
WORKDIR /app

# Copy the entire project
COPY . .

# Download Eigen if not present (replacing the Windows batch script with direct download)
RUN if [ ! -d "eigen-3.4.0" ]; then \
        echo "Downloading Eigen library..." && \
        wget -O eigen.tar.bz2 "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2" && \
        tar -xjf eigen.tar.bz2 && \
        rm eigen.tar.bz2; \
    fi

# Build the project using the Linux Makefile
RUN make -f Makefile.linux clean && make -f Makefile.linux

# Expose port for Jupyter notebook (if needed)
EXPOSE 8888

# Default command - show usage
CMD ["./WignerDyson.exe", "--help"]