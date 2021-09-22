CC = g++

all:
	$(CC) -mfma -fopenmp -O3 main.cpp src/WignerDyson.cpp src/Orthogonal.cpp src/Unitary.cpp src/Symplectic.cpp src/Quantum_chaotic_billiard.cpp -o WignerDyson

clean:
	rm WignerDyson
