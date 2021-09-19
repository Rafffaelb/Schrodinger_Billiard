#include <iostream>
#include <cmath>
#include <complex>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include "../include/Quantum_chaotic_billiard.h"

using namespace std;

Quantum_chaotic_billiard::Quantum_chaotic_billiard(MatrixXcd H, MatrixXcd W, MatrixXcd C1, MatrixXcd C2){
	Set_Setup(H, W, C1, C2);
}

void Quantum_chaotic_billiard::Set_Setup(MatrixXcd H, MatrixXcd W, MatrixXcd C1, MatrixXcd C2)
{
	_H = H;
	_W = W;
	_C1 = C1;
	_C2 = C2;
}

void Quantum_chaotic_billiard::Calculate_Smatrix(){
	
	complex<double> number_2(2,0);
	complex<double> complex_identity(0,1);

	int ress = _H.rows();
	int N1 = (_C1.rows())/2;
	int N2 = (_C2.rows())/2;
	int n = N1+N2;
	
	MatrixXcd identityS = MatrixXcd::Identity(n,n);

	MatrixXcd D(_H.rows(), _H.cols());

	D << (-_H + complex_identity*M_PI*_W*(_W.adjoint()));
	PartialPivLU<MatrixXcd> lu(D);
	MatrixXcd D_inv_W = lu.inverse()*_W;

	// Scattering Matrix //

	MatrixXcd S(n,n);

	S << identityS - number_2*complex_identity*M_PI*(_W.adjoint())*D_inv_W;
	
	this -> _S = S;
	
}

void Quantum_chaotic_billiard::Calculate_G_and_P(){

	MatrixXcd ttdaga = _C1*_S*_C2*(_S.adjoint());

	MatrixXcd identityP = MatrixXcd::Identity(ttdaga.rows(),ttdaga.cols());

	_G = ttdaga.trace();
	_P = (ttdaga*(identityP-ttdaga)).trace();

}

complex<double> Quantum_chaotic_billiard::getG(){
	return this -> _G;
}

complex<double> Quantum_chaotic_billiard::getP(){
	return this -> _P;
}
