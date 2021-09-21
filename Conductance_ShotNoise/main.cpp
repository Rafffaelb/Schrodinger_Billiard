#include <iostream>
#include <cstring>
#include "include/WignerDyson.h"
#include <cmath>
#include "include/Orthogonal.h"
#include "include/Unitary.h"
#include "include/Symplectic.h"

using namespace std;

int main(int argc, char **argv){

	double Gamma, lambda, y, V, gama;
	int N1, N2, n, ress, num_steps, spin_deg;

	Gamma = 1;
	ress = 100;
	lambda = 0.5;
	y = sqrt(1.0/Gamma)*(1.0-sqrt(1.0-Gamma));
	V = lambda*lambda/ress;
	num_steps = 100000;

	for (int i = 1; i < argc; i++){
		
		if (strcmp(argv[i],"Orthogonal") == 0){
			
			spin_deg = 1;
			cout << "\n ####### Running Orthogonal ####### \n" << endl;
			Orthogonal orthogonal(Gamma, ress, lambda, num_steps, spin_deg);
			orthogonal.Run_Simulation();
		}
		else{
			if (strcmp(argv[i],"Unitary") == 0){
			
				spin_deg = 1;
				cout << "\n ####### Running Unitary ####### \n" << endl;
				Unitary unitary(Gamma, ress, lambda, num_steps, spin_deg);
				unitary.Run_Simulation();
		}
			else{
				if (strcmp(argv[i],"Symplectic") == 0){
					spin_deg = 2;
					cout << "\n ####### Run Symplectic ###### \n" << endl;
					Symplectic symplectic(Gamma, ress, lambda, num_steps, spin_deg);
					symplectic.Run_Simulation();
				}
			}
		}
	}
	return 0;
}
