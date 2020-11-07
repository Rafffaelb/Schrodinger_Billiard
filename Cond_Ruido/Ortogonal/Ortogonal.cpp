#include <iostream>
#include <cmath>
#include <complex>
#include <random>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <ctime>
#include <fstream>
#include <chrono>
#include "../../include/Ham_O.h"
#include "../../include/W_O.h"
#include "../../include/Projetores_O.h"

using namespace Eigen;
using namespace std::literals;

int main(){

	std::complex<double> complex_identity(0, 1);
	std::complex<double> number_2(2, 0);

	// Input //

	double Gamma, lambda, y, V, gama;
	int N1, N2, n, ress, num_realization, lambda1;

	Gamma = 1;
	N1 = 1;
        N2 = 1;
    	n = N1+N2;
	ress = 100;
	lambda = 0.5;
	lambda1 = 1;
	y = sqrt(1.0/Gamma)*(1.0-sqrt(1.0-Gamma));
	V = lambda*lambda/ress;
	num_realization = 50000;
	
	MatrixXcd identityS(n,n);
	MatrixXcd G(num_realization,1);
	MatrixXcd R(num_realization,1);

	// Pauli Matrices //

	MatrixXcd matrizpauli1(2,2);
	MatrixXcd matrizpauli2(2,2);
	MatrixXcd matrizpauli3(2,2);

	matrizpauli1.real() << 1, 0, 0, 1;
	matrizpauli1.imag() <<  0, 0, 0,  0;

	matrizpauli2.real() << 0, 0, 0, 0;
	matrizpauli2.imag() <<  0, -1, 1,  0;

	
	matrizpauli3.real() << 1, 0, 0, -1;
	matrizpauli3.imag() <<  0, 0, 0,  0;

	// Creating W Matrices //

	MatrixXcd W(ress,n);
	W.setZero();
	MatrixXcd* W_pointer = &W;
	
	Criando_W(W_pointer, ress, N1, N2, lambda, y);

	// Creating Projectors //

	MatrixXcd C1(2*N1, 2*N1);
	MatrixXcd C2(2*N2, 2*N2);

	C1.setZero();
	C2.setZero();

	MatrixXcd *C1_pointer = &C1;
	MatrixXcd *C2_pointer = &C2;

	Criando_Projetores(C1_pointer, C2_pointer, N1, N2);
	
	for (int realization = 1; realization < num_realization + 1; realization++){

		// Generating Hamiltonian Matrix //

		MatrixXcd H(ress,ress);
		H.setZero();
		MatrixXcd* H_pointer = &H;
		Criando_H(H_pointer, ress, V);
		
		// Inverse Green Function //

		MatrixXcd D = (-H + complex_identity*M_PI*W*(W.adjoint()));
		
		// Scattering Matrix //

		MatrixXcd S = identityS -number_2*complex_identity*M_PI*(W.adjoint())*(D.inverse())*W;

		//std::cout << "\nA matriz S fica:\n" << S << std::endl;
		//std::cout << "\nA matriz S transposta fica:\n" << S.transpose() << std::endl;

		MatrixXcd ttdaga = C1*S*C2*(S.adjoint());

		// Conductance and Power Shot Noise //

		MatrixXcd identityR = MatrixXcd::Identity(ttdaga.rows(),ttdaga.cols());

		G(realization-1, 0) = ttdaga.trace();
		R(realization-1, 0) = (ttdaga*(identityR-ttdaga)).trace();

		if (realization % 10000 == 0){
			std::cout << realization << std::endl;
		}
	}

	std::ofstream output_G("G_O.txt");
	std::ofstream output_R("R_O.txt");
	for (int i = 0; i < num_realization; i++){
		for (int j = 0; j < 1; j++){
			output_G << G(i,j).real() << "\t" << std::endl;
			output_R << R(i,j).real() << "\t" << std::endl;
		}
	}
	return 0;
}
