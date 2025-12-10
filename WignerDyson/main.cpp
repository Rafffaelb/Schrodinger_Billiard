#include <iostream>
#include <cstring>
#include "include/WignerDyson_AbstractClass_h/WignerDyson.h"
#include <cmath>
#include "include/Orthogonal.h"
#include "include/Unitary.h"
#include "include/Symplectic.h"

using namespace std;

// Function to display usage information
void displayUsage(const char* programName) {
    cout << "Usage: " << programName << " <symmetry_class> <simulation_type>\n\n";
    cout << "Symmetry classes:\n";
    cout << "  Orthogonal     Gaussian Orthogonal Ensemble (GOE) - Time-reversal symmetric\n";
    cout << "  Unitary        Gaussian Unitary Ensemble (GUE) - Time-reversal broken\n";
    cout << "  Symplectic     Gaussian Symplectic Ensemble (GSE) - Time-reversal symmetric with spin-orbit coupling\n\n";
    
    cout << "Simulation types:\n";
    cout << "  Channel                    Conductance vs channel count\n";
    cout << "  Gamma                      Conductance vs gamma parameter\n";
    cout << "  Concurrence                Entanglement concurrence vs gamma\n";
    cout << "  Bell_Parameter_Ress        Bell parameter vs Ress\n";
    cout << "  Bell_Parameter_Gamma       Bell parameter vs gamma\n";
    cout << "  Bell_Parameter_Fixed_Base  Bell parameter with fixed basis\n";
    cout << "  Energy                     Conductance vs energy\n";
    cout << "  Energy_Gamma               Conductance vs energy with gamma variation\n";
}

// Function to run orthogonal simulations
void runOrthogonalSimulations(const char* simulationType) {
    double lambda = 0.5;
    int spin_deg = 1;
    
    Orthogonal orthogonal(lambda, spin_deg);
    
    cout << "\n ####### Running Orthogonal Symmetry Class ####### \n" << endl;
    
    if (strcmp(simulationType, "Channel") == 0) {
        cout << "Simulation: Conductance vs Channel Count" << endl;
        orthogonal.Run_Simulation_Conductance_Channels();
    } else if (strcmp(simulationType, "Gamma") == 0) {
        cout << "Simulation: Conductance vs Gamma Parameter" << endl;
        orthogonal.Run_Simulation_Conductance_Gamma();
    } else if (strcmp(simulationType, "Concurrence") == 0) {
        cout << "Simulation: Concurrence vs Gamma Parameter" << endl;
        orthogonal.Run_Simulation_Concurrence_Gamma();
    } else if (strcmp(simulationType, "Bell_Parameter_Ress") == 0) {
        cout << "Simulation: Bell Parameter vs Ress" << endl;
        orthogonal.Run_Simulation_Bell_Parameter_Ress();
    } else if (strcmp(simulationType, "Bell_Parameter_Gamma") == 0) {
        cout << "Simulation: Bell Parameter vs Gamma" << endl;
        orthogonal.Run_Simulation_Bell_Parameter_Gamma();
    } else if (strcmp(simulationType, "Bell_Parameter_Fixed_Base") == 0) {
        cout << "Simulation: Bell Parameter with Fixed Basis" << endl;
        orthogonal.Run_Simulation_Bell_Parameter_Fixed_Base();
    } else if (strcmp(simulationType, "Energy") == 0) {
        cout << "Simulation: Conductance vs Energy" << endl;
        orthogonal.Run_Simulation_Conductance_Energy();
    } else if (strcmp(simulationType, "Energy_Gamma") == 0) {
        cout << "Simulation: Conductance vs Energy with Gamma Variation" << endl;
        orthogonal.Run_Simulation_Conductance_Energy_Gamma();
    } else {
        cout << "Error: Unknown simulation type '" << simulationType << "' for Orthogonal symmetry class.\n" << endl;
        displayUsage("WignerDyson");
    }
}

// Function to run unitary simulations
void runUnitarySimulations(const char* simulationType) {
    double lambda = 0.5;
    int spin_deg = 1;
    
    Unitary unitary(lambda, spin_deg);
    
    cout << "\n ####### Running Unitary Symmetry Class ####### \n" << endl;
    
    if (strcmp(simulationType, "Channel") == 0) {
        cout << "Simulation: Conductance vs Channel Count" << endl;
        unitary.Run_Simulation_Conductance_Channels();
    } else if (strcmp(simulationType, "Gamma") == 0) {
        cout << "Simulation: Conductance vs Gamma Parameter" << endl;
        unitary.Run_Simulation_Conductance_Gamma();
    } else if (strcmp(simulationType, "Concurrence") == 0) {
        cout << "Simulation: Concurrence vs Gamma Parameter" << endl;
        unitary.Run_Simulation_Concurrence_Gamma();
    } else if (strcmp(simulationType, "Bell_Parameter_Ress") == 0) {
        cout << "Simulation: Bell Parameter vs Ress" << endl;
        unitary.Run_Simulation_Bell_Parameter_Ress();
    } else if (strcmp(simulationType, "Bell_Parameter_Gamma") == 0) {
        cout << "Simulation: Bell Parameter vs Gamma" << endl;
        unitary.Run_Simulation_Bell_Parameter_Gamma();
    } else if (strcmp(simulationType, "Bell_Parameter_Fixed_Base") == 0) {
        cout << "Simulation: Bell Parameter with Fixed Basis" << endl;
        unitary.Run_Simulation_Bell_Parameter_Fixed_Base();
    } else if (strcmp(simulationType, "Energy") == 0) {
        cout << "Simulation: Conductance vs Energy" << endl;
        unitary.Run_Simulation_Conductance_Energy();
    } else if (strcmp(simulationType, "Energy_Gamma") == 0) {
        cout << "Simulation: Conductance vs Energy with Gamma Variation" << endl;
        unitary.Run_Simulation_Conductance_Energy_Gamma();
    } else {
        cout << "Error: Unknown simulation type '" << simulationType << "' for Unitary symmetry class.\n" << endl;
        displayUsage("WignerDyson");
    }
}

// Function to run symplectic simulations
void runSymplecticSimulations(const char* simulationType) {
    double lambda = 0.5;
    int spin_deg = 2;
    
    Symplectic symplectic(lambda, spin_deg);
    
    cout << "\n ####### Running Symplectic Symmetry Class ####### \n" << endl;
    
    if (strcmp(simulationType, "Channel") == 0) {
        cout << "Simulation: Conductance vs Channel Count" << endl;
        symplectic.Run_Simulation_Conductance_Channels();
    } else if (strcmp(simulationType, "Gamma") == 0) {
        cout << "Simulation: Conductance vs Gamma Parameter" << endl;
        symplectic.Run_Simulation_Conductance_Gamma();
    } else {
        cout << "Error: Unknown simulation type '" << simulationType << "' for Symplectic symmetry class.\n" << endl;
        cout << "Supported Symplectic simulations: Channel, Gamma\n" << endl;
        displayUsage("WignerDyson");
    }
}

int main(int argc, char **argv) {
    // Check for help flag
    if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)) {
        displayUsage(argv[0]);
        return 0;
    }
    
    // Check minimum arguments
    if (argc != 3) {
        cout << "Error: Incorrect number of arguments.\n" << endl;
        displayUsage(argv[0]);
        return 1;
    }
    
    const char* symmetryClass = argv[1];
    const char* simulationType = argv[2];
    
    // Run the appropriate simulation
    if (strcmp(symmetryClass, "Orthogonal") == 0) {
        runOrthogonalSimulations(simulationType);
    } else if (strcmp(symmetryClass, "Unitary") == 0) {
        runUnitarySimulations(simulationType);
    } else if (strcmp(symmetryClass, "Symplectic") == 0) {
        runSymplecticSimulations(simulationType);
    } else {
        cout << "Error: Unknown symmetry class '" << symmetryClass << "'.\n" << endl;
        displayUsage(argv[0]);
        return 1;
    }
    
    return 0;
}