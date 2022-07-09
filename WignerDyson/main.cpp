#include <iostream>
#include <cstring>
#include "include/WignerDyson_AbstractClass_h/WignerDyson.h"
#include <cmath>
#include "include/Orthogonal.h"
#include "include/Unitary.h"
#include "include/Symplectic.h"

using namespace std;

int main(int argc, char **argv){

	double lambda;
	int num_steps, spin_deg;

	lambda = 0.5;
	num_steps = 2000;

	for (int i = 1; i < argc; i++){
		
		if (strcmp(argv[i],"Orthogonal") == 0){
		
			spin_deg = 1;
			
			Orthogonal orthogonal(lambda, num_steps, spin_deg);

			for (int j = 1; j < argc; j++){
				
				if (strcmp(argv[j],"Channel") == 0){
					
					cout << "\n ####### Running Orthogonal (variable: Channel) ####### \n" << endl;
				
					orthogonal.Run_Simulation_Conductance_Channels();
				}

				if (strcmp(argv[j],"Gamma") == 0){
	
					cout << "\n ###### Running Orthogonal (variable: Gamma) ###### \n" << endl;
				
					orthogonal.Run_Simulation_Conductance_Gamma();
				}

				if (strcmp(argv[j],"Concurrence") == 0){

					cout << "\n ###### Running Orthogonal Concurrence (variable: Gamma) ###### \n" << endl;
				
					orthogonal.Run_Simulation_Concurrence_Gamma();
				}

				if (strcmp(argv[j],"Bell_Parameter_Ress") == 0){
					
					cout << "\n ###### Running Orthogonal Bell Parameter (variable: Ress) ##### \n" << endl;
				
					orthogonal.Run_Simulation_Bell_Parameter_Ress();
				}

				if (strcmp(argv[j],"Bell_Parameter_Gamma") == 0){
		
					cout << "\n ###### Running Orthogonal Bell Parameter (variable: Gamma) #### \n" << endl;
					
					orthogonal.Run_Simulation_Bell_Parameter_Gamma();
				}

				if (strcmp(argv[j],"Bell_Parameter_Fixed_Base") == 0){
				
					cout << "\n ###### Running Orthogonal Bell Parameter Fixed Base ##### \n" << endl;

					orthogonal.Run_Simulation_Bell_Parameter_Fixed_Base();
				}

				if (strcmp(argv[j],"Energy") == 0){
	
					cout << "\n ###### Running Orthogonal (variable: Energy) ###### \n" << endl;
				
					orthogonal.Run_Simulation_Conductance_Energy();
				}
			}
			
			orthogonal.~Orthogonal();
		}
		else{
			if (strcmp(argv[i],"Unitary") == 0){
		
				spin_deg = 1;

				Unitary unitary(lambda, num_steps, spin_deg);

				for (int j = 1; j < argc; j++){
				
					if (strcmp(argv[j],"Channel") == 0){
					
						cout << "\n ####### Running Unitary (variable: Channel) ####### \n" << endl;
					
						unitary.Run_Simulation_Conductance_Channels();
					}
			
					if (strcmp(argv[j],"Gamma") == 0){
	
						cout << "\n ###### Running Unitary (variable: Gamma) ###### \n" << endl;
					
						unitary.Run_Simulation_Conductance_Gamma();
					}
			
					if (strcmp(argv[j],"Concurrence") == 0){

						cout << "\n ###### Running Unitary Concurrence (variable: Gamma) ###### \n" << endl;
						
						unitary.Run_Simulation_Concurrence_Gamma();
					}
					if (strcmp(argv[j],"Bell_Parameter_Ress") == 0){
					
						cout << "\n ###### Running Unitary Bell Parameter (variable: Ress) ##### \n" << endl;
					
						unitary.Run_Simulation_Bell_Parameter_Ress();
					}
					if (strcmp(argv[j],"Bell_Parameter_Gamma") == 0){
						
						cout << "\n ##### Running Unitary Bell Parameter (variable: Gamma) ##### \n" << endl;

						unitary.Run_Simulation_Bell_Parameter_Gamma();
					}

					if (strcmp(argv[j],"Bell_Parameter_Fixed_Base") == 0){

						cout << "\n ##### Running Unitary Bell Parameter Fixed Base ##### \n" << endl;

						unitary.Run_Simulation_Bell_Parameter_Fixed_Base();
					}
				}
				
				unitary.~Unitary();
			}
			else{
				if (strcmp(argv[i],"Symplectic") == 0){
		
					spin_deg = 2;

					Symplectic symplectic(lambda, num_steps, spin_deg);

					for (int j = 1; j < argc; j++){
				
						if (strcmp(argv[j],"Channel") == 0){
					
							cout << "\n ####### Running Symplectic (variable: Channel) ####### \n" << endl;
						
							symplectic.Run_Simulation_Conductance_Channels();
						
						}
						if (strcmp(argv[j],"Gamma") == 0){
						
							cout << "\n ###### Running Symplectic (variable: Gamma) ###### \n" << endl;
						
							symplectic.Run_Simulation_Conductance_Gamma();
						
						}
					}

					symplectic.~Symplectic();
				}
			}
		}
	}
	return 0;
}
