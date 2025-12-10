# Docker Guide for Schrodinger Billiard Project

This guide explains how to use Docker with the Schrodinger Billiard project. Docker provides a consistent environment for building and running the project, regardless of your local system setup.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Building the Docker Image](#building-the-docker-image)
- [Running Simulations](#running-simulations)
- [Accessing Results](#accessing-results)
- [Using Jupyter Notebooks](#using-jupyter-notebooks)
- [Advanced Docker Usage](#advanced-docker-usage)
- [Troubleshooting](#troubleshooting)

## Prerequisites

1. Docker Desktop installed on your system
   - For Windows/Mac: Download from https://www.docker.com/products/docker-desktop
   - For Linux: Follow the official installation guide for your distribution

2. At least 4GB of free disk space

3. Internet connection for initial build (to download base image and dependencies)

## Building the Docker Image

From the `WignerDyson` directory, build the Docker image:

```bash
cd WignerDyson
docker build -t schrodinger-billiard .
```

This will:
- Download the Ubuntu 20.04 base image
- Install all required dependencies (GCC, Python, etc.)
- Copy your project files into the container
- Build the project using the Linux-compatible Makefile

The build process may take several minutes on first run.

## Running Simulations

### Interactive Mode

Run the container interactively to explore available options:

```bash
docker run -it schrodinger-billiard
```

### Direct Simulation Execution

Run a specific simulation directly:

```bash
docker run schrodinger-billiard ./WignerDyson.exe Orthogonal Channel
```

Available symmetry classes:
- `Orthogonal` - Gaussian Orthogonal Ensemble (GOE)
- `Unitary` - Gaussian Unitary Ensemble (GUE)
- `Symplectic` - Gaussian Symplectic Ensemble (GSE)

Available simulation types:
- `Channel` - Conductance vs channel count
- `Gamma` - Conductance vs gamma parameter
- `Concurrence` - Entanglement concurrence vs gamma
- `Bell_Parameter_Ress` - Bell parameter vs Ress
- `Bell_Parameter_Gamma` - Bell parameter vs gamma
- `Bell_Parameter_Fixed_Base` - Bell parameter with fixed basis
- `Energy` - Conductance vs energy
- `Energy_Gamma` - Conductance vs energy with gamma variation

Example simulations:
```bash
# Run unitary conductance vs gamma simulation
docker run schrodinger-billiard ./WignerDyson.exe Unitary Gamma

# Run symplectic Bell parameter vs Ress simulation
docker run schrodinger-billiard ./WignerDyson.exe Symplectic Bell_Parameter_Ress
```

## Accessing Results

To access simulation results on your host machine, use volume mounting:

```bash
docker run -v "$(pwd)/Data_Analysis:/app/Data_Analysis" schrodinger-billiard ./WignerDyson.exe Orthogonal Channel
```

On Windows (PowerShell):
```powershell
docker run -v "${PWD}/Data_Analysis:/app/Data_Analysis" schrodinger-billiard ./WignerDyson.exe Orthogonal Channel
```

On Windows (Command Prompt):
```cmd
docker run -v "%cd%\Data_Analysis:/app/Data_Analysis" schrodinger-billiard ./WignerDyson.exe Orthogonal Channel
```

This mounts the container's `/app/Data_Analysis` directory to your host's `Data_Analysis` directory, allowing you to access the generated `.dat` files.

## Using Jupyter Notebooks

The Docker image includes Jupyter Notebook for analyzing results:

```bash
docker run -p 8888:8888 -v "$(pwd)/Data_Analysis:/app/Data_Analysis" schrodinger-billiard jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser --allow-root
```

Then open your browser and go to http://localhost:8888. You'll need the token displayed in the terminal output.

Available notebooks for analysis:
- `Bell_Parameter.ipynb`
- `Concurrence_Gamma.ipynb`
- `Conductance_Energy.ipynb`
- `Conductance_ShotNoise.ipynb`

## Advanced Docker Usage

### Using Docker Compose

For easier management, use the provided `docker-compose.yml`:

```bash
# Build the image
docker-compose build

# Run a simulation
docker-compose run schrodinger-billiard ./WignerDyson.exe Orthogonal Channel

# Run Jupyter Notebook
docker-compose run -p 8888:8888 schrodinger-billiard jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser --allow-root
```

### Using Helper Scripts

We provide helper scripts for common operations:

Linux/macOS:
```bash
chmod +x docker_helper.sh
./docker_helper.sh build
./docker_helper.sh exec "./WignerDyson.exe Orthogonal Channel"
```

Windows:
```cmd
docker_helper.bat build
docker_helper.bat exec "./WignerDyson.exe Orthogonal Channel"
```

### Running Long Simulations

For long-running simulations, run in detached mode:

```bash
docker run -d --name my-simulation -v "$(pwd)/Data_Analysis:/app/Data_Analysis" schrodinger-billiard ./WignerDyson.exe Unitary Energy_Gamma
```

Check the logs:
```bash
docker logs my-simulation
```

Stop the simulation:
```bash
docker stop my-simulation
```

### Resource Limiting

Limit resources for simulations:

```bash
docker run -d --cpus="2" --memory="4g" --name simulation schrodinger-billiard ./WignerDyson.exe Orthogonal Channel
```

## Troubleshooting

### Common Issues

1. **"Cannot connect to the Docker daemon"**
   - Make sure Docker Desktop is running
   - On Linux, you might need to run with `sudo` or add your user to the docker group

2. **Build fails with file not found errors**
   - Ensure you're running the command from the `WignerDyson` directory
   - Check that all project files are present

3. **Permission denied when accessing results**
   - Files created by Docker may be owned by root
   - On Linux, you might need to adjust permissions: `sudo chown -R $USER:$USER Data_Analysis`

4. **Port already in use when running Jupyter**
   - Change the host port: `-p 8889:8888`

### Getting Help

Show program usage:
```bash
docker run schrodinger-billiard ./WignerDyson.exe --help
```

View container logs:
```bash
docker logs <container_name_or_id>
```

Execute commands in a running container:
```bash
docker exec -it <container_name_or_id> /bin/bash
```

## Cleaning Up

Remove the Docker image:
```bash
docker rmi schrodinger-billiard
```

Remove all stopped containers:
```bash
docker container prune
```

Remove unused volumes:
```bash
docker volume prune
```