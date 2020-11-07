#include <iostream>
#include <cmath>
#include <complex>
#include <random>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <ctime>
#include <chrono>
#include "../include/Ham_O.h"

using namespace std;
using namespace Eigen;
using namespace std::literals;

void Criando_H (MatrixXcd *H_pointer, int ress, double V){
	
	auto seed = std::chrono::system_clock::now().time_since_epoch().count();	
  	std::normal_distribution<double> distribution(0.0,1.0);
	std::default_random_engine generator(seed);
	
	MatrixXcd A(ress,ress); 
	MatrixXcd H1(ress,ress);
	MatrixXcd Simetrica(ress,ress);

	A.setZero();
	H1.setZero();
	Simetrica.setZero();
	
	for (int i = 1; i < ress + 1; i++){
		for (int j = 1; j < ress + 1; j++){
			double aux = distribution(generator);
			A(i-1,j-1) = aux;
		}
	}

	for (int i = 1; i < ress + 1; i++){
		H1(i-1,i-1) = A(i-1,i-1)*sqrt(V/(2.0));
		for (int j = i + 1; j < ress + 1; j++){
			H1(i-1,j-1) = A(i-1,j-1)*sqrt(V/(1.0));
		}
	}

	Simetrica << H1 + H1.adjoint();

	MatrixXcd H = Simetrica;
	*H_pointer = H;
	
}

