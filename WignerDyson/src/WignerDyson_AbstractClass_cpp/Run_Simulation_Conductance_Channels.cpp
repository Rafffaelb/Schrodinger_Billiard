#include <iostream>
#include <chrono>
#include "../../include/WignerDyson_AbstractClass_h/WignerDyson.h"
#include "../../include/Quantum_chaotic_billiard.h"
#include "../../include/WignerDyson_AbstractClass_h/Run_Simulation_Conductance_Channels.h"
#include <eigen3/Eigen/Dense>
#include <omp.h>

using namespace std;
using namespace Eigen;

void WignerDyson::Run_Simulation_Conductance_Channels(){

	auto start = chrono::system_clock::now();

	double Gamma, y, V;
       	int ress;

	ress = 100;

	Gamma = 1;
	y = sqrt(1.0/Gamma)*(1.0-sqrt(1.0-Gamma));
	V = _lambda*_lambda/ress;

	MatrixXcd G(_num_steps, 10);
	MatrixXcd P(_num_steps, 10);

	G.setZero();
	P.setZero();

	for (int N1 = 1; N1 < 11; N1++){
		
		int N2 = N1;
		int n = N1 + N2;

		// Create W Matrices //

		MatrixXcd W(_spin_deg * ress, _spin_deg * n);
		W.setZero();
		MatrixXcd* W_pointer = &W;

		Create_W(W_pointer, ress, N1, N2, _lambda, y);

		// Create_ProjectionMatrices //

		MatrixXcd C1(_spin_deg * 2*N1, _spin_deg * 2*N1); MatrixXcd C2(_spin_deg * 2*N2, _spin_deg * 2*N2);
		C1.setZero(); C2.setZero();
		MatrixXcd *C1_pointer = &C1; MatrixXcd *C2_pointer = &C2;

		Create_ProjectionMatrices(C1_pointer, C2_pointer, N1, N2);
		
		#pragma omp parallel for shared(W, C1, C2)
		for (int step = 1; step < _num_steps + 1; step++){
		
			// Generate Hamiltonian Matrix //

			MatrixXcd H(_spin_deg * ress, _spin_deg * ress);
			H.setZero();
			MatrixXcd* H_pointer = &H;

			Create_H(H_pointer, ress, V);

			// Create billiard setup //
		
			Quantum_chaotic_billiard billiard_setup(H, W, C1, C2);

			// Scattering Matrix //
		
			billiard_setup.Calculate_Smatrix();

			// Conductance (G) and Power Shot Noise (P) //
		
			billiard_setup.Calculate_G_and_P();

			G(step-1, N1-1) = billiard_setup.getG();
			P(step-1, N1-1) = billiard_setup.getP();
	
			if (step % 50000 == 0){
				std::cout << "\nCurrent number of steps: " << step << "| Current number of open channels (N): " << N1 << std::endl;
			}
		}
	
		//Save G and P matrices as txt files //
	
		Save_txt_files_Channels(G, P, _num_steps);
	}

	auto end = chrono::system_clock::now();
	auto elapsed =
		chrono::duration_cast<chrono::minutes>(end-start);
	cout << "\n Simulation time duration: " << elapsed.count() << "\n";
}

