#include <iostream>
#include <cmath>
#include <complex>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include "../include/W.h"

using namespace Eigen;

void Create_W (MatrixXcd *W_pointer, int ress, int N1, int N2, double lambda, double y){

	MatrixXcd W1(ress,N1);
	MatrixXcd W2(ress,N2);
	MatrixXcd W(ress,N1+N2);

	for (int j=1; j < ress+1; j++ ){
		for (int k=1; k < N1+1; k++){
			std::complex<double> aux(y*(sqrt(((2.0*lambda))/(M_PI*(ress+1)))*sin(j*k*M_PI/(ress+1))), 0);
			W1(j-1,k-1) = aux;
		}
	}

	for (int j=1; j < ress+1; j++ ){
		for (int k=1; k < N2+1; k++){
			std::complex<double> aux(y*(sqrt(((2.0*lambda))/(M_PI*(ress+1)))*sin(j*(k+N1)*M_PI/(ress+1))), 0);
			W2(j-1,k-1) = aux;
		}
	}

	W << W1, W2;
	*W_pointer = W;

}

void Create_W_Symplectic (MatrixXcd *W_pointer, int ress, int N1, int N2, double lambda, double y){
	
	MatrixXcd identity2x2(2,2);

	identity2x2 << MatrixXcd::Identity(2,2);

	MatrixXcd W1(ress,N1);
	MatrixXcd W2(ress,N2);
	MatrixXcd W_aux1(W1.rows(), W1.cols() + W2.cols());
	MatrixXcd W_aux2(2*W1.rows(), 2*(W1.cols() + W2.cols()));

	for (int j = 1; j < ress + 1; j++){
		for (int k = 1; k < N1 + 1; k++){
			complex<double> aux(y*(sqrt(((2.0*lambda))/(M_PI*(ress + 1)))*sin(j*k*M_PI/(ress+1))),0);
			W1(j-1,k-1) = aux;
		}
	}

	for (int j = 1; j < ress + 1; j++){
		for (int k = 1; k < N2 + 1; k++){
			complex<double> aux(y*(sqrt(((2.0*lambda))/(M_PI*(ress + 1)))*sin(j*(k+N1)*M_PI/(ress+1))),0);
			W2(j-1,k-1) = aux;
		}
	}

	W_aux1 << W1, W2;

	for (int i = 1; i < W_aux1.rows() + 1; i++){
		for (int j = 1; j < W_aux1.cols() + 1; j++){
			W_aux2.block((i-1)*identity2x2.rows(), (j-1)*identity2x2.cols(), identity2x2.rows(), identity2x2.cols()) = W_aux1(i-1,j-1)*identity2x2;
		}
	}

	MatrixXcd W = W_aux2;
	*W_pointer = W;
}
