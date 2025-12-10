#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <chrono>
#include "../include/WignerDyson_AbstractClass_h/WignerDyson.h"
#include "../include/Quantum_chaotic_billiard.h"

using namespace std;
using namespace Eigen;

void WignerDyson::Run_Simulation_Bell_Parameter_Ress(){

	auto start = chrono::system_clock::now();

	double Gamma, y, V;
       	int ress, ress_idx, N1, N2, n, _num_steps;

	Gamma = 1;
	y = sqrt(double(1.0)/Gamma)*(1.0-sqrt(1.0-Gamma));

	ress = 50;
	_num_steps = 1000000;

	N1 = 2;
	N2 = N1;
	n = N1 + N2;

	// Create_ProjectionMatrices //

	MatrixXcd C1(_spin_deg * 2*N1, _spin_deg * 2*N1); MatrixXcd C2(_spin_deg * 2*N2, _spin_deg * 2*N2);
	C1.setZero(); C2.setZero();
	MatrixXcd *C1_pointer = &C1; MatrixXcd *C2_pointer = &C2;

	Create_ProjectionMatrices(C1_pointer, C2_pointer, N1, N2);
	
	MatrixXd Bell_Parameter_Ress(_num_steps, 11);

	Bell_Parameter_Ress.setZero();

	for (ress_idx = 1; ress_idx < 12; ress_idx++){

		V = _lambda*_lambda/ress;

		// Create W Matrices //

		MatrixXcd W(_spin_deg * ress, _spin_deg * n);
		W.setZero();
		MatrixXcd *W_pointer = &W;

		Create_W(W_pointer, ress, N1, N2, _lambda, y);

		#pragma omp parallel for shared(W, C1, C2)
		for (int step = 1; step < _num_steps + 1; step++){
		
			// Generate Hamiltonian Matrix //

			MatrixXcd H(_spin_deg * ress, _spin_deg * ress);
			H.setZero();
			MatrixXcd *H_pointer = &H;

			Create_H(H_pointer, ress, V);

			// Create billiard setup //

			Quantum_chaotic_billiard billiard_setup(H, W, C1, C2);

			// Scattering Matrix //
		
			billiard_setup.Calculate_Smatrix(0);

			// Bell Parameter //
		
			billiard_setup.Calculate_Bell_Parameter();

			Bell_Parameter_Ress(step-1, ress_idx-1) = billiard_setup.getBell_Parameter();

			if (step % _num_steps == 0){
				std::cout << "\nCurrent number of steps: " << step << "| Current index of Ress: " << ress_idx << std::endl;
			}
		}
		
		//Save Concurrence matrix as txt files //
		Save_txt_files_Bell_Parameter_Ress(Bell_Parameter_Ress, _num_steps);

		ress = ress + 25;
	}

	auto end = chrono::system_clock::now();
	auto elapsed =
		chrono::duration_cast<chrono::minutes>(end-start);
	cout << "\n Simulation time duration: " << elapsed.count() << "\n";
}


