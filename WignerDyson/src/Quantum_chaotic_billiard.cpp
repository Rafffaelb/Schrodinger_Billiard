#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <chrono>
#include <random>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include "../include/Quantum_chaotic_billiard.h"

// Define M_PI if not available
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

MatrixXcd Kronecker_Product(MatrixXcd A, MatrixXcd B){

	MatrixXcd C(A.rows() * B.rows(), A.cols() * B.rows());

	for (int i = 0; i < A.rows(); i++){
		for (int j = 0; j < A.cols(); j++){
			C.block(i*B.rows(), j*B.cols(), B.rows(), B.cols()) = A(i,j)*B;
		}
	}

	return C;
}

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

void Quantum_chaotic_billiard::Calculate_Smatrix(double Energy){
	
	complex<double> number_2(2,0);
	complex<double> complex_identity(0,1);

	int ress = _H.rows();
	int N1 = (_C1.rows())/2;
	int N2 = (_C2.rows())/2;
	int n = N1+N2;
	
	MatrixXcd identityS = MatrixXcd::Identity(n,n);

	MatrixXcd D(ress, ress);

	D << (Energy*MatrixXcd::Identity(ress,ress) -_H + complex_identity*M_PI*_W*(_W.adjoint()));
	PartialPivLU<MatrixXcd> lu(D);
	MatrixXcd D_inv_W = lu.inverse()*_W;

	// Test //

	// MatrixXcd x = D.partialPivLu().solve(_W);
	// double relative_error = (D*x - _W).norm() / _W.norm(); // norm() is L2 norm
	// std::cout << "\n The relative error is:\n" << relative_error << std::endl;
	

	// Scattering Matrix //

	MatrixXcd S(n,n);

	S << identityS - number_2*complex_identity*M_PI*(_W.adjoint())*D_inv_W;
	
	this -> _S = S;
	//cout << "\nH = \n" << _H << endl;
	//cout << "\nS = \n" << S << endl;

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

MatrixXcd Create_Unitary_Random_Matrix();
complex<double> Calculate_Noise(MatrixXcd r, MatrixXcd t, MatrixXcd U_L, MatrixXcd U_R);
complex<double> Calculate_Noise_Fixed_Base(MatrixXcd r, MatrixXcd t, MatrixXcd U_L, MatrixXcd U_R);

void Quantum_chaotic_billiard::Calculate_Bell_Parameter(){

	MatrixXcd r, t, U_L(2,2), U_R(2,2), U_Lprime(2,2), U_Rprime(2,2), paulimatrix_z(2,2);
	complex<double> C_a_b, C_a_bprime, C_aprime_b, C_aprime_bprime;

	paulimatrix_z << 1, 0,
			 0, -1;

	const int N1 = (_C1.rows())/2;
	const int N2 = (_C2.rows())/2;

	t = _S.block(N1,0,N2,N1);
	r = _S.block(0,0,N1,N1);

	U_L = Create_Unitary_Random_Matrix();
	U_R = Create_Unitary_Random_Matrix();
	U_Lprime = Create_Unitary_Random_Matrix();
	U_Rprime = Create_Unitary_Random_Matrix();

	C_a_b = Calculate_Noise(r, t, U_L, U_R);
	C_a_bprime = Calculate_Noise(r, t, U_L, U_Rprime);
	C_aprime_b = Calculate_Noise(r, t, U_Lprime, U_R);
	C_aprime_bprime = Calculate_Noise(r, t, U_Lprime, U_Rprime);

	_Bell_Parameter = abs(((C_a_b + C_aprime_b + C_a_bprime - C_aprime_bprime).real()));
	_Bell_Parameter_Dephase = 2*abs((paulimatrix_z*r*t.adjoint()*paulimatrix_z*t*r.adjoint()).trace())/((r.adjoint()*r*t.adjoint()*t).trace()).real();
}

void Quantum_chaotic_billiard::Calculate_Bell_Parameter_Fixed_Base(){

	MatrixXcd r, t, U_L(2,2), U_R(2,2), U_Lprime(2,2), U_Rprime(2,2), paulimatrix_x(2,2), paulimatrix_z(2,2);
	complex<double> C_a_b, C_a_bprime, C_aprime_b, C_aprime_bprime;

	paulimatrix_x << 0, 1,
		      	 1, 0;

	paulimatrix_z << 1, 0,
			 0, -1;

	const int N1 = (_C1.rows())/2;
	const int N2 = (_C2.rows())/2;

	t = _S.block(N1,0,N2,N1);
	r = _S.block(0,0,N1,N1);

	U_L = paulimatrix_z;
	U_R = -(1/sqrt(2))*(paulimatrix_x + paulimatrix_z);
	U_Lprime = paulimatrix_x;
	U_Rprime = (1/sqrt(2))*(paulimatrix_z - paulimatrix_x);

	C_a_b = Calculate_Noise_Fixed_Base(r, t, U_L, U_R);
	C_a_bprime = Calculate_Noise_Fixed_Base(r, t, U_L, U_Rprime);
	C_aprime_b = Calculate_Noise_Fixed_Base(r, t, U_Lprime, U_R);
	C_aprime_bprime = Calculate_Noise_Fixed_Base(r, t, U_Lprime, U_Rprime);

	_Bell_Parameter = abs(((C_a_b + C_aprime_b + C_a_bprime - C_aprime_bprime).real()));
}

double Quantum_chaotic_billiard::getBell_Parameter(){
	return this -> _Bell_Parameter;
}

double Quantum_chaotic_billiard::getBell_Parameter_Dephase(){
	return this -> _Bell_Parameter_Dephase;
}

MatrixXcd Create_Unitary_Random_Matrix(){

	// Function to Create Unitary Random Matrix distributed with Haar Measure //

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

	Q = Q*Delta;
       
	return Q; // Unitary random matrix distributed with Haar measure //
}

complex<double> Calculate_Noise(MatrixXcd r, MatrixXcd t, MatrixXcd U_L, MatrixXcd U_R){
	
	complex<double> C, Const_Norm;

	MatrixXcd paulimatrix_z(2,2), a_dot_sigma(2,2), b_dot_sigma(2,2);

	paulimatrix_z << 1, 0,
			 0, -1;

	a_dot_sigma = U_L.adjoint()*paulimatrix_z*U_L;
	b_dot_sigma = U_R.adjoint()*paulimatrix_z*U_R;

	Const_Norm = ((r*r.adjoint()).trace())*((t*t.adjoint()).trace())-(r*t.adjoint()*t*r.adjoint()).trace();

	C = ((((a_dot_sigma)*r*r.adjoint()).trace())*(((b_dot_sigma)*t*t.adjoint()).trace()) - (((a_dot_sigma)*r*t.adjoint()*(b_dot_sigma)*t*r.adjoint()).trace()))/Const_Norm;

	return C;
}

complex<double> Calculate_Noise_Fixed_Base(MatrixXcd r, MatrixXcd t, MatrixXcd U_L, MatrixXcd U_R){
	
	complex<double> C, Const_Norm;

	MatrixXcd paulimatrix_z(2,2), a_dot_sigma(2,2), b_dot_sigma(2,2);

	paulimatrix_z << 1, 0,
			 0, -1;

	a_dot_sigma = U_L;
	b_dot_sigma = U_R;

	Const_Norm = ((r*r.adjoint()).trace())*((t*t.adjoint()).trace())-(r*t.adjoint()*t*r.adjoint()).trace();

	C = ((((a_dot_sigma)*r*r.adjoint()).trace())*(((b_dot_sigma)*t*t.adjoint()).trace()) - (((a_dot_sigma)*r*t.adjoint()*(b_dot_sigma)*t*r.adjoint()).trace()))/Const_Norm;

	return C;
}