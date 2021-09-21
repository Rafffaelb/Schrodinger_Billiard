#include <iostream>
#include <cmath>
#include <complex>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <fstream>
#include <ctime>
#include <chrono>
#include "../../include/Ham_O.h"
#include "../../include/W.h"
#include "../../include/ProjectionMatrices.h"
#include "../../include/Quantum_chaotic_billiard.h"
#include "../../include/Save_txt_files.h"
#include "omp.h"

using namespace Eigen;

int main(){

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	std::complex<double> complex_identity(0, 1);
	std::complex<double> number_2(2, 0);

	// Input //
	
	double Gamma, lambda, y, V, gama;
	int N1, N2, n, ress, num_steps;
	
	Gamma = 1;
	ress = 100;
	lambda = 0.5;
	y = sqrt(1.0/Gamma)*(1.0-sqrt(1.0-Gamma));
	V = lambda*lambda/ress;
	num_steps = 100000;
	
	MatrixXcd G(num_steps,10);
	MatrixXcd P(num_steps,10);
	
	G.setZero();
	P.setZero();

	for (int N1 = 1; N1 < 11; N1++ ){
        
		N2 = N1;
    		n = N1+N2;
		
		// Create W Matrices //

		MatrixXcd W(ress,n);
		W.setZero();
		MatrixXcd* W_pointer = &W;
	
		Create_W(W_pointer, ress, N1, N2, lambda, y);

		// Create Projectors //

		MatrixXcd C1(2*N1, 2*N1); MatrixXcd C2(2*N2, 2*N2);
		C1.setZero(); C2.setZero();
		MatrixXcd *C1_pointer = &C1; MatrixXcd *C2_pointer = &C2;

		Create_ProjectionMatrices(C1_pointer, C2_pointer, N1, N2);
		
		#pragma omp parallel for shared(W, C1, C2)	
		for (int step = 1; step < num_steps+1; step++){

			// Generate Hamiltonian Matrix //
			
			MatrixXcd H_O(ress,ress);
			H_O.setZero();
			MatrixXcd* H_pointer = &H_O;

			Create_H(H_pointer, ress, V);

			// Create billiard setup //

			Quantum_chaotic_billiard billiard_setup(H_O, W, C1, C2);

			// Scattering Matrix //
			
			billiard_setup.Calculate_Smatrix();
		
			// Calculate Conductance (G) and Power Shot Noise (P) //

			billiard_setup.Calculate_G_and_P();

			G(step-1, N1-1) = billiard_setup.getG();
			P(step-1, N1-1) = billiard_setup.getP();
			
			if (step % 50000 == 0){
				std::cout << "\nCurrent number of steps: " << step << " | Current number of open channels (N): " << N1 << std::endl;
			}
		}

		// Save G and P matrices as txt files //

		Save_txt_files(G, P, num_steps);

	}

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	std::cout << "\nTime: " << std::chrono::duration_cast<std::chrono::minutes>(end-begin).count() << " [minutes]" << std::endl;

	return 0;
}
