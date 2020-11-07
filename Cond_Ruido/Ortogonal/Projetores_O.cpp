#include <iostream>
#include <cmath>
#include <complex>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include "../../include/Projetores_O.h"

using namespace Eigen;
using namespace std::literals;

void Criando_Projetores (MatrixXcd *C1_pointer, MatrixXcd *C2_pointer, int N1, int N2){

	MatrixXcd C1tio(2,2);
	MatrixXcd C2tio(2,2);

	MatrixXcd identidade1 = MatrixXcd::Identity(N1,N1);
	MatrixXcd identidade2 = MatrixXcd::Identity(N2,N2);
	
	MatrixXcd C1(2*N1, 2*N1);
	MatrixXcd C2(2*N2, 2*N2);

	for (int i=1; i < 3; i++){
		for (int j=1; j < 3; j++){
			if (i == 1 && j == 1){
				std::complex<double> aux(1,0);
				C1tio(i-1,j-1) = aux;
			}
			else{
				std::complex<double> aux(0,0);
				C1tio(i-1,j-1) = aux;
			}
		}
	}

	for (int i=1; i < 3; i++){
		for (int j=1; j < 3; j++){
			if (i == 2 && j == 2){
				std::complex<double> aux(1,0);
				C2tio(i-1,j-1) = aux;
			}
			else{
				std::complex<double> aux(0,0);
				C2tio(i-1,j-1) = aux;
			}
		}
	}

	for (int i = 1; i < 3; i++){
		for (int j = 1; j < 3; j++){
			C1.block((i-1)*identidade1.rows(), (j-1)*identidade1.cols(), identidade1.rows(), identidade1.cols()) = C1tio(i-1,j-1)*identidade1;
			
		}
	}

	for (int i = 1; i < 3; i++){
		for (int j = 1; j < 3; j++){
			C2.block((i-1)*identidade2.rows(), (j-1)*identidade2.cols(), identidade2.rows(), identidade2.cols()) = C2tio(i-1,j-1)*identidade2;
		}
	}

	*C1_pointer << C1;
	*C2_pointer << C2;

}
