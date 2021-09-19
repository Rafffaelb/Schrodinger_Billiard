#include <iostream>
#include <cmath>
#include <complex>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include "../include/ProjectionMatrices.h"

using namespace Eigen;

void Create_ProjectionMatrices (MatrixXcd *C1_pointer, MatrixXcd *C2_pointer, int N1, int N2){

	MatrixXcd C1tio(2,2);
	MatrixXcd C2tio(2,2);

	MatrixXcd identity1 = MatrixXcd::Identity(N1,N1);
	MatrixXcd identity2 = MatrixXcd::Identity(N2,N2);
	
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
			C1.block((i-1)*identity1.rows(), (j-1)*identity1.cols(), identity1.rows(), identity1.cols()) = C1tio(i-1,j-1)*identity1;
			
		}
	}

	for (int i = 1; i < 3; i++){
		for (int j = 1; j < 3; j++){
			C2.block((i-1)*identity2.rows(), (j-1)*identity2.cols(), identity2.rows(), identity2.cols()) = C2tio(i-1,j-1)*identity2;
		}
	}

	*C1_pointer << C1;
	*C2_pointer << C2;

}

void Create_ProjectionMatrices_Symplectic(MatrixXcd *C1_pointer, MatrixXcd *C2_pointer, int N1, int N2){
	
	MatrixXcd identity2N1 = MatrixXcd::Identity(2*N1,2*N1);
	MatrixXcd identity2N2 = MatrixXcd::Identity(2*N2,2*N2);

	MatrixXcd C1tio(2,2);
	MatrixXcd C2tio(2,2);

	for (int i = 1; i < 3; i++){
		for (int j = 1; j < 3; j++){
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

	for (int i = 1; i < 3; i++){
		for (int j = 1; j < 3; j++){
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

	MatrixXcd C1(identity2N1.rows() * C1tio.rows(), identity2N1.cols() * C1tio.cols());
	MatrixXcd C2(identity2N2.rows() * C2tio.rows(), identity2N2.cols() * C2tio.cols());

	C1.setZero();
	C2.setZero();

	for (int i = 1; i < 3; i++){
		for (int j = 1; j < 3; j++){
			C1.block((i-1)*identity2N1.rows(), (j-1)*identity2N1.cols(), identity2N1.rows(), identity2N1.cols()) = C1tio(i-1,j-1)*identity2N1;
		}
	}

	for (int i = 1; i < 3; i++){
		for (int j = 1; j < 3; j++){
			C2.block((i-1)*identity2N2.rows(), (j-1)*identity2N2.cols(), identity2N2.rows(), identity2N2.cols()) = C2tio(i-1,j-1)*identity2N2;
		}
	}

	*C1_pointer << C1;
	*C2_pointer << C2;

}
