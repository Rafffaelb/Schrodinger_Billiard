#include <iostream>
#include <cmath>
#include <complex>
#include <random>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <ctime>
#include <chrono>
#include "../include/Ham_S.h"

using namespace std;
using namespace Eigen;

complex<double> complex_identity1(0,1);

void Create_H (MatrixXcd *H_pointer, int ress, double V){
	
	auto seed = std::chrono::system_clock::now().time_since_epoch().count();	
  	std::normal_distribution<double> distribution(0.0,1.0);
	std::default_random_engine generator(seed);

	MatrixXcd identity2x2_aux = MatrixXcd::Identity(2,2);

	MatrixXcd paulimatrix1(2,2);
	MatrixXcd paulimatrix2(2,2);
	MatrixXcd paulimatrix3(2,2);

	paulimatrix1.real() << 0, 1, 1, 0;
	paulimatrix1.imag() << 0, 0, 0, 0;

	paulimatrix2.real() << 0, 0, 0, 0;
	paulimatrix2.imag() << 0, -1, 1, 0;

	paulimatrix3.real() << 1, 0, 0, -1;
	paulimatrix3.imag() << 0, 0, 0, 0;

	MatrixXcd A(ress,ress); MatrixXcd C(ress,ress);
	A.setZero(); C.setZero();

	for (int i = 1; i < ress + 1; i++){
		for (int j = 1; j < ress + 1; j++){
			complex<double> aux_1 = distribution(generator);
			A(i-1,j-1) = aux_1;
		}
	}

	for (int i = 1; i < ress + 1; i++){
		for (int j = 1; j < ress + 1; j++){
			complex<double> aux_1 = distribution(generator);
			C(i-1,j-1) = aux_1;
		}
	}

	MatrixXcd H0(ress,ress); MatrixXcd H1(ress,ress);
	MatrixXcd H2(ress, ress); MatrixXcd H3(ress,ress);

	H0.setZero(); H1.setZero();
	H2.setZero(); H3.setZero();

	for (int i = 1; i < ress + 1; i++){
		H0(i-1,i-1) = A(i-1,i-1)*sqrt(V/(2.0));
		H1(i-1,i-1) = 0.0;
		H2(i-1,i-1) = 0.0;
		H3(i-1,i-1) = 0.0;
		for (int j = i + 1; j < ress + 1; j++){
			H0(i-1,j-1) = A(i-1,j-1)*sqrt(V/(4.0));
			H1(i-1,j-1) = A(j-1,i-1)*sqrt(V/(4.0));
			H2(i-1,j-1) = C(j-1,i-1)*sqrt(V/(4.0));
			H3(i-1,j-1) = C(i-1,j-1)*sqrt(V/(4.0));
		}
	}

	
	MatrixXcd Symmetric(ress,ress); MatrixXcd Antisymmetric1(ress,ress);
	MatrixXcd Antisymmetric2(ress,ress); MatrixXcd Antisymmetric3(ress,ress);

	Symmetric.setZero(); Antisymmetric1.setZero();
	Antisymmetric2.setZero(); Antisymmetric3.setZero();

	Symmetric << H0 + H0.adjoint(); Antisymmetric1 << H1 - H1.adjoint();
	Antisymmetric2 << H2 - H2.adjoint(); Antisymmetric3 << H3 - H3.adjoint();

	MatrixXcd Q0(2*ress,2*ress); MatrixXcd Q1(2*ress,2*ress);
	MatrixXcd Q2(2*ress,2*ress); MatrixXcd Q3(2*ress,2*ress);

	for (int i = 1; i < Symmetric.rows() + 1; i++){
		for (int j = 1; j < Symmetric.cols() + 1; j++){
			Q0.block((i-1)*identity2x2_aux.rows(), (j-1)*identity2x2_aux.cols(), identity2x2_aux.rows(), identity2x2_aux.cols()) = Symmetric(i-1,j-1)*identity2x2_aux;
			Q1.block((i-1)*paulimatrix1.rows(), (j-1)*paulimatrix1.cols(), paulimatrix1.rows(), paulimatrix1.cols()) = Antisymmetric1(i-1,j-1)*paulimatrix1;
			Q2.block((i-1)*paulimatrix2.rows(), (j-1)*paulimatrix2.cols(), paulimatrix2.rows(), paulimatrix2.cols()) = Antisymmetric2(i-1,j-1)*paulimatrix2;
			Q3.block((i-1)*paulimatrix3.rows(), (j-1)*paulimatrix3.cols(), paulimatrix3.rows(), paulimatrix3.cols()) = Antisymmetric3(i-1,j-1)*paulimatrix3;
		}
	}

	MatrixXcd H(2*ress,2*ress);
	H << Q0 + complex_identity1*(Q1 + Q2 + Q3);
	*H_pointer = H;
}

