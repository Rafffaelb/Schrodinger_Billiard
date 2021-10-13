#include <iostream>
#include <cmath>
#include <ctime>
#include <chrono>
#include <random>
#include <complex>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/QR>
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

void Quantum_chaotic_billiard::Calculate_Concurrence(){

	const int N1 = (_C1.rows())/2;
	const int N2 = (_C2.rows())/2;

	MatrixXcd t = _S.block(N1,0,N2,N1);

	MatrixXcd ttdaga = t*t.adjoint();

	VectorXcd eigenvalues_ttdaga = ttdaga.eigenvalues();

	double tau_1 = eigenvalues_ttdaga(0).real();
	double tau_2 = eigenvalues_ttdaga(1).real();

	_Concurrence = 2*(sqrt(tau_1*(1-tau_1)*tau_2*(1-tau_2))/(tau_1+tau_2-2*tau_1*tau_2));

	_Entanglement = -((1+sqrt(1-pow(_Concurrence,2)))/2)*log2((1+sqrt(1-pow(_Concurrence,2)))/2) - (1-(1+sqrt(1-pow(_Concurrence,2)))/2)*log2(1-(1+sqrt(1-pow(_Concurrence,2)))/2);
}

double Quantum_chaotic_billiard::getConcurrence(){
	return this -> _Concurrence;
}

double Quantum_chaotic_billiard::getEntanglement(){
	return this -> _Entanglement;
}

void Haar_measure();

void Quantum_chaotic_billiard::Calculate_Bell_Parameter_Ress(){

	const int N1 = (_C1.rows())/2;
	const int N2 = (_C2.rows())/2;

	MatrixXcd t = _S.block(N1,0,N2,N1);
	MatrixXcd r = _S.block(0,0,N1,N1);

	Haar_measure();

}

double Quantum_chaotic_billiard::getBell_Parameter_Ress(){
	return this -> _Concurrence;
}

void Haar_measure(){

	complex<double> complex_identity(0,1);

	MatrixXcd Z(2,2), A(2,2), B(2,2), Q(2,2), R(2,2);
	MatrixXcd Diag_R, Delta;

	ColPivHouseholderQR<MatrixXcd> qr(Z.rows(), Z.cols());

	auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::normal_distribution<double> distribution(0.0,1.0);
	std::default_random_engine generator(seed);
	
	for (int i = 1; i < 3; i++){
		for (int j = 1; j < 3; j++){
			double aux = distribution(generator);
			A(i-1,j-1) = aux;
		}
	}

	for (int i = 1; i < 3; i++){
			for (int j = 1; j < 3; j++){
			double aux = distribution(generator);
			B(i-1,j-1) = aux;
		}
	}

	Z = (1/sqrt(2))*(A+complex_identity*B);

	qr.compute(Z);

	Q = qr.householderQ().setLength(qr.nonzeroPivots());
	R = qr.matrixR().template triangularView<Upper>();

	Diag_R = R.diagonal().matrix().asDiagonal();

	Delta = Diag_R*Diag_R.cwiseAbs().inverse();

	Q = Q*Delta; // Unitary random matrix distributed with Haar measure //


}


