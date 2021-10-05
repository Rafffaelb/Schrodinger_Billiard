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
	num_steps = 100000;

	for (int i = 1; i < argc; i++){
		
		if (strcmp(argv[i],"Orthogonal") == 0){
		
			spin_deg = 1;
			for (int j = 1; j < argc; j++){
				
				if (strcmp(argv[j],"Channel") == 0){
					
					cout << "\n ####### Running Orthogonal (variable: Channel) ####### \n" << endl;
					Orthogonal orthogonal(lambda, num_steps, spin_deg);
					orthogonal.Run_Simulation_Conductance_Channels();
					orthogonal.~Orthogonal();
				}
				if (strcmp(argv[j],"Gamma") == 0){
	
					cout << "\n ###### Running Orthogonal (variable: Gamma) ###### \n" << endl;
					Orthogonal orthogonal(lambda, num_steps, spin_deg);
					orthogonal.Run_Simulation_Conductance_Gamma();
					orthogonal.~Orthogonal();
				}
			}
		}
		else{
			if (strcmp(argv[i],"Unitary") == 0){
		
				spin_deg = 1;
				for (int j = 1; j < argc; j++){
				
					if (strcmp(argv[j],"Channel") == 0){
					
						cout << "\n ####### Running Unitary (variable: Channel) ####### \n" << endl;
						Unitary unitary(lambda, num_steps, spin_deg);
						unitary.Run_Simulation_Conductance_Channels();
						unitary.~Unitary();
					}
					if (strcmp(argv[j],"Gamma") == 0){
	
						cout << "\n ###### Running Unitary (variable: Gamma) ###### \n" << endl;
						Unitary unitary(lambda, num_steps, spin_deg);
						unitary.Run_Simulation_Conductance_Gamma();
						unitary.~Unitary();
					}
				}
			}
			else{
				if (strcmp(argv[i],"Symplectic") == 0){
		
					spin_deg = 2;
					for (int j = 1; j < argc; j++){
				
						if (strcmp(argv[j],"Channel") == 0){
					
							cout << "\n ####### Running Symplectic (variable: Channel) ####### \n" << endl;
							Symplectic symplectic(lambda, num_steps, spin_deg);
							symplectic.Run_Simulation_Conductance_Channels();
							symplectic.~Symplectic();
						}
						if (strcmp(argv[j],"Gamma") == 0){
						
							cout << "\n ###### Running Symplectic (variable: Gamma) ###### \n" << endl;
							Symplectic symplectic(lambda, num_steps, spin_deg);
							symplectic.Run_Simulation_Conductance_Gamma();
							symplectic.~Symplectic();
						}
					}
				}
			}
		}
	}
	return 0;
}
