#include <iostream>
#include <cmath>
#include <complex>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <fstream>
#include <ctime>
#include <chrono>
#include "../../include/Ham_O.h"
#include "../../include/W_O.h"
#include "../../include/ProjectionMatrices_O.h"
#include "../../include/Quantum_chaotic_billiard_O.h"
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
	MatrixXcd R(num_steps,10);
	
	G.setZero();
	R.setZero();
	
	for (int N1 = 1; N1 < 11; N1++ ){
        
		N2 = N1;
    		n = N1+N2;

		MatrixXcd identityS = MatrixXcd::Identity(n,n);
		
		// Creating W Matrices //

		MatrixXcd W(ress,n);
		W.setZero();
		MatrixXcd* W_pointer = &W;
	
		Create_W(W_pointer, ress, N1, N2, lambda, y);

		// Creating Projectors //

		MatrixXcd C1(2*N1, 2*N1); MatrixXcd C2(2*N2, 2*N2);
		C1.setZero(); C2.setZero();
		MatrixXcd *C1_pointer = &C1; MatrixXcd *C2_pointer = &C2;

		Create_ProjectionMatrices(C1_pointer, C2_pointer, N1, N2);
		
		#pragma omp parallel for	
		for (int step = 1; step < num_steps+1; step++){

			// Generating Hamiltonian Matrix //
			
			MatrixXcd H(ress,ress);
			H.setZero();
			MatrixXcd* H_pointer = &H;

			Create_H(H_pointer, ress, V);

			Quantum_chaotic_billiard billiard_setup(H, W, C1, C2);
			billiard_setup.Calculate_Smatrix();
			// billiard_setup.calculate_conductance();
			// billiard_setup.calculate_power_shot_noise();

			// MatrixXcd ttdaga = C1*S*C2*(S.adjoint());

			// Conductance and Power Shot Noise //

			// MatrixXcd identityR = MatrixXcd::Identity(ttdaga.rows(),ttdaga.cols());

			// G(step-1, N1-1) = ttdaga.trace();
			// R(step-1, N1-1) = (ttdaga*(identityR-ttdaga)).trace();
			
			if (step % 50000 == 0){
				std::cout << "\nCurrent number of steps: " << step << " | Current number of open channels (N): " << N1 << std::endl;
			}
		}

		std::ofstream output_G("G_O.txt");
		std::ofstream output_R("R_O.txt");
		for (int i = 0; i < num_steps; i++){
			for (int j = 0; j < 10; j++){
				if (j == 9){
					output_G << G(i,j).real() << std::endl;
					output_R << R(i,j).real() << std::endl;
				}
				else{
					output_G << G(i,j).real() << "\t";
					output_R << R(i,j).real() << "\t";
				}
			}
		}
	}

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	std::cout << "\nTime: " << std::chrono::duration_cast<std::chrono::minutes>(end-begin).count() << " [minutes]" << std::endl;

	return 0;
}
