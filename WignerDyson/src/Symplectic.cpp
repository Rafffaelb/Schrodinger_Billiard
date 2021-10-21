#include <iostream>
#include "../include/Symplectic.h"
#include <cmath>
#include <complex>
#include <ctime>
#include <chrono>
#include <random>
#include <fstream>
#include <string>

using namespace std;

Symplectic::Symplectic(double lambda, int num_steps, int spin_deg){

	this -> _lambda = lambda;
	this -> _num_steps = num_steps;
	this -> _spin_deg = spin_deg;
}

Symplectic::~Symplectic() {}

void Symplectic::Create_W(MatrixXcd* W_pointer, int ress, int N1, int N2, double lambda, double y){

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


void Symplectic::Create_ProjectionMatrices(MatrixXcd* C1_pointer, MatrixXcd* C2_pointer, int N1, int N2){
	

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

void Symplectic::Create_H(MatrixXcd* H_pointer, int ress, double V){

	auto seed = std::chrono::system_clock::now().time_since_epoch().count();	
  	std::normal_distribution<double> distribution(0.0,1.0);
	std::default_random_engine generator(seed);

	complex<double> complex_identity(0,1);

	MatrixXcd identity2x2 = MatrixXcd::Identity(2,2);

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
			Q0.block((i-1)*identity2x2.rows(), (j-1)*identity2x2.cols(), identity2x2.rows(), identity2x2.cols()) = Symmetric(i-1,j-1)*identity2x2;
			Q1.block((i-1)*paulimatrix1.rows(), (j-1)*paulimatrix1.cols(), paulimatrix1.rows(), paulimatrix1.cols()) = Antisymmetric1(i-1,j-1)*paulimatrix1;
			Q2.block((i-1)*paulimatrix2.rows(), (j-1)*paulimatrix2.cols(), paulimatrix2.rows(), paulimatrix2.cols()) = Antisymmetric2(i-1,j-1)*paulimatrix2;
			Q3.block((i-1)*paulimatrix3.rows(), (j-1)*paulimatrix3.cols(), paulimatrix3.rows(), paulimatrix3.cols()) = Antisymmetric3(i-1,j-1)*paulimatrix3;
		}
	}

	MatrixXcd H(2*ress,2*ress);
	H << Q0 + complex_identity*(Q1 + Q2 + Q3);
	*H_pointer = H;		
}

void Symplectic::Save_txt_files_Channels(MatrixXcd G, MatrixXcd P, int num_steps){
	std::ofstream output_G("Data_Analysis/Channel/G_S_Channel.txt");
	std::ofstream output_P("Data_Analysis/Channel/P_S_Channel.txt");

	for(int i = 0; i < num_steps; i++){
		for (int j = 0; j < 10; j++){
			if (j == 9){
				output_G << G(i,j).real() << std::endl;
				output_P << P(i,j).real() << std::endl;
			}
			else{
				output_G << G(i,j).real() << "\t";
				output_P << P(i,j).real() << "\t";
			}
		}
	}	
}

void Symplectic::Save_txt_files_Gamma(MatrixXcd G, MatrixXcd P, int num_steps, int N1){
	std::ofstream output_G("Data_Analysis/Gamma/G_S_Gamma_N"+to_string(N1)+".txt");
	std::ofstream output_P("Data_Analysis/Gamma/P_S_Gamma_N"+to_string(N1)+".txt");

	for(int i = 0; i < num_steps; i++){
		for (int j = 0; j < 21; j++){
			if (j == 20){
				output_G << G(i,j).real() << std::endl;
				output_P << P(i,j).real() << std::endl;
			}
			else{
				output_G << G(i,j).real() << "\t";
				output_P << P(i,j).real() << "\t";
			}
		}
	}	
}

void Symplectic::Save_txt_files_Concurrence_Gamma(MatrixXd Concurrence, MatrixXd Entanglement, int num_steps){

}

void Symplectic::Save_txt_files_Bell_Parameter_Ress(MatrixXd Bell_Parameter_Ress, int num_steps){

}

void Symplectic::Save_txt_files_Bell_Parameter_Gamma(MatrixXd Bell_Parameter_Gamma, MatrixXd Bell_Parameter_Dephase_Gamma, int num_steps){

}
