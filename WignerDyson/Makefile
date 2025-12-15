CC = g++
CXXFLAGS = -mfma -fopenmp -O3 -std=c++11 -I. -Ieigen-3.4.0 -D_USE_MATH_DEFINES
EIGEN_DIR = eigen-3.4.0

# Source files
SRCS = main.cpp \
       src/WignerDyson_AbstractClass_cpp/WignerDyson.cpp \
       src/Orthogonal.cpp \
       src/Unitary.cpp \
       src/Symplectic.cpp \
       src/Quantum_chaotic_billiard.cpp \
       src/WignerDyson_AbstractClass_cpp/Run_Simulation_Conductance_Channels.cpp \
       src/WignerDyson_AbstractClass_cpp/Run_Simulation_Conductance_Gamma.cpp \
       src/WignerDyson_AbstractClass_cpp/Run_Simulation_Concurrence_Gamma.cpp \
       src/WignerDyson_AbstractClass_cpp/Run_Simulation_Bell_Parameter_Ress.cpp \
       src/WignerDyson_AbstractClass_cpp/Run_Simulation_Bell_Parameter_Gamma.cpp \
       src/WignerDyson_AbstractClass_cpp/Run_Simulation_Bell_Parameter_Fixed_Base.cpp \
       src/WignerDyson_AbstractClass_cpp/Run_Simulation_Conductance_Energy.cpp \
       src/WignerDyson_AbstractClass_cpp/Run_Simulation_Conductance_Energy_Gamma.cpp

OBJS = $(SRCS:.cpp=.o)
TARGET = WignerDyson

.PHONY: all clean clean-all eigen-download eigen-check

all: eigen-check $(TARGET)

eigen-check:
	@if [ ! -d "$(EIGEN_DIR)" ]; then \
		echo "Eigen library not found. Downloading..."; \
		make eigen-download; \
	fi

eigen-download:
	@chmod +x ./download_eigen_linux.sh
	@./download_eigen_linux.sh

$(TARGET): $(OBJS)
	$(CC) $(CXXFLAGS) -o $@ $^ -lm

%.o: %.cpp
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o src/*.o src/WignerDyson_AbstractClass_cpp/*.o $(TARGET)

clean-all: clean
	rm -rf $(EIGEN_DIR) $(TARGET)