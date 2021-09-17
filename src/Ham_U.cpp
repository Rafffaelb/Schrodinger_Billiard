#include <iostream>
#include <cmath>
#include <complex>
#include <random>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <ctime>
#include <chrono>
#include "../include/Ham_U.h"

using namespace std;
using namespace Eigen;

void Create_H (MatrixXcd *H_pointer, int ress, double V){

	complex<double> complex_identity(0,1); 	
		
	auto seed = std::chrono::system_clock::now().time_since_epoch().count();	
  	std::normal_distribution<double> distribution(0.0,1.0);
	std::default_random_engine generator(seed);
	
	MatrixXcd A(ress,ress); 
	MatrixXcd H1(ress,ress); MatrixXcd H2(ress,ress);
	MatrixXcd Symmetric(ress,ress); MatrixXcd Antisymmetric(ress,ress);

	A.setZero();
	H1.setZero(); H2.setZero();
	Symmetric.setZero(); Antisymmetric.setZero();
	
	for (int i = 1; i < ress + 1; i++){
		for (int j = 1; j < ress + 1; j++){
			double aux = distribution(generator);
			A(i-1,j-1) = aux;
		}
	}

	for (int i = 1; i < ress + 1; i++){
		H1(i-1,i-1) = A(i-1,i-1)*sqrt(V/(2.0));
		H2(i-1,i-1) = 0.0;
		for (int j = i + 1; j < ress + 1; j++){
			H1(i-1,j-1) = A(i-1,j-1)*sqrt(V/(2.0));
			H2(i-1,j-1) = A(j-1,i-1)*sqrt(V/(2.0));
		}
	}

	Symmetric << H1 + H1.adjoint();
	Antisymmetric << H2 - H2.adjoint();

	MatrixXcd H = Symmetric + complex_identity*Antisymmetric;
	*H_pointer = H;
	
}

