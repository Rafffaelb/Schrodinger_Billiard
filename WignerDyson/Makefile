CC = g++

all:
	$(CC) -mfma -fopenmp -O3 main.cpp src/WignerDyson_AbstractClass_cpp/WignerDyson.cpp src/Orthogonal.cpp src/Unitary.cpp src/Symplectic.cpp src/Quantum_chaotic_billiard.cpp src/WignerDyson_AbstractClass_cpp/Run_Simulation_Conductance_Channels.cpp src/WignerDyson_AbstractClass_cpp/Run_Simulation_Conductance_Gamma.cpp src/WignerDyson_AbstractClass_cpp/Run_Simulation_Concurrence_Gamma.cpp src/WignerDyson_AbstractClass_cpp/Run_Simulation_Bell_Parameter_Ress.cpp src/WignerDyson_AbstractClass_cpp/Run_Simulation_Bell_Parameter_Gamma.cpp src/WignerDyson_AbstractClass_cpp/Run_Simulation_Bell_Parameter_Fixed_Base.cpp -o WignerDyson

clean:
	rm WignerDyson
