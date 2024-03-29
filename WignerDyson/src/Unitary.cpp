#include <iostream>
#include "../include/Unitary.h"
#include <cmath>
#include <complex>
#include <random>
#include <ctime>
#include <chrono>
#include <fstream>
#include <string>

using namespace std;

Unitary::Unitary(double lambda, int spin_deg){

	this -> _lambda = lambda;
	this -> _spin_deg = spin_deg;
}

Unitary::~Unitary() {}

void Unitary::Create_W(MatrixXcd* W_pointer, int ress, int N1, int N2, double lambda, double y){

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

void Unitary::Create_ProjectionMatrices(MatrixXcd* C1_pointer, MatrixXcd* C2_pointer, int N1, int N2){

	MatrixXcd identity1 = MatrixXcd::Identity(N1,N1);
	MatrixXcd identity2 = MatrixXcd::Identity(N2,N2);
	
	MatrixXcd C1((N1+N2), (N1+N2));
	MatrixXcd C2((N1+N2), (N1+N2));

	C1.block(0, 0, N1, N1) << identity1; C1.block(0, N1, N1, N2) << MatrixXcd::Zero(N1, N2);
	C1.block(N1, 0, N2, N1) << MatrixXcd::Zero(N2, N1); C1.block(N1, N1, N2, N2) << MatrixXcd::Zero(N2, N2);

	C2.block(0, 0, N1, N1) << MatrixXcd::Zero(N1, N1); C2.block(0, N1, N1, N2) << MatrixXcd::Zero(N1, N2);
	C2.block(N1, 0, N2, N1) << MatrixXcd::Zero(N2, N1); C2.block(N1, N1, N2, N2) << identity2;

	*C1_pointer << C1;
	*C2_pointer << C2;
}

void Unitary::Create_H(MatrixXcd* H_pointer, int ress, double V){

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

void Unitary::Save_txt_files_Channels(MatrixXcd G, MatrixXcd P, int num_steps){
	
	std::ofstream output_G("Data_Analysis/Channel/G_U_Channel.txt");
	std::ofstream output_P("Data_Analysis/Channel/P_U_Channel.txt");

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

void Unitary::Save_txt_files_Gamma(MatrixXcd G, MatrixXcd P, int num_steps, int N1){
	std::ofstream output_G("Data_Analysis/Gamma/G_U_Gamma_N"+to_string(N1)+".txt");
	std::ofstream output_P("Data_Analysis/Gamma/P_U_Gamma_N"+to_string(N1)+".txt");

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

void Unitary::Save_txt_files_Concurrence_Gamma(MatrixXd Concurrence, MatrixXd Entanglement, int num_steps){

	std::ofstream output_Concurrence("Data_Analysis/Concurrence/Concurrence_U_Gamma.txt");
	std::ofstream output_Entanglement("Data_Analysis/Concurrence/Entanglement_U_Gamma.txt");

	for(int i = 0; i < num_steps; i++){
		for (int j = 0; j < 21; j++){
			if (j == 20){
				output_Concurrence << Concurrence(i,j) << std::endl;
				output_Entanglement << Entanglement(i,j) << std::endl;
			}
			else{
				output_Concurrence << Concurrence(i,j) << "\t";
				output_Entanglement << Entanglement(i,j) << "\t";
			}
		}
	}	
}

void Unitary::Save_txt_files_Bell_Parameter_Ress(MatrixXd Bell_Parameter_Ress, int num_steps){

	std::ofstream output_Bell_Parameter_Ress("Data_Analysis/Bell_Parameter/Bell_Ress/Bell_Parameter_U_Ress.txt");

	for(int i = 0; i < num_steps; i++){
		for (int j = 0; j < 11; j++){
			if (j == 10){
				output_Bell_Parameter_Ress << Bell_Parameter_Ress(i,j) << std::endl;
			}
			else{
				output_Bell_Parameter_Ress << Bell_Parameter_Ress(i,j) << "\t";
			}
		}
	}
}

void Unitary::Save_txt_files_Bell_Parameter_Gamma(MatrixXd Bell_Parameter_Gamma, MatrixXd Bell_Parameter_Dephase_Gamma, int num_steps){

	std::ofstream output_Bell_Parameter_Gamma("Data_Analysis/Bell_Parameter/Bell_Gamma/Bell_Parameter_U_Gamma.txt");
	std::ofstream output_Bell_Parameter_Dephase_Gamma("Data_Analysis/Bell_Parameter/Bell_Gamma/Bell_Parameter_Dephase_U_Gamma.txt");

	for(int i = 0; i < num_steps; i++){
		for (int j = 0; j < 21; j++){
			if (j == 20){
				output_Bell_Parameter_Gamma << Bell_Parameter_Gamma(i,j) << std::endl;
				output_Bell_Parameter_Dephase_Gamma << Bell_Parameter_Dephase_Gamma(i,j) << std::endl;
			}
			else{
				output_Bell_Parameter_Gamma << Bell_Parameter_Gamma(i,j) << "\t";
				output_Bell_Parameter_Dephase_Gamma << Bell_Parameter_Dephase_Gamma(i,j) << "\t";
			}
		}
	}
}

void Unitary::Save_txt_files_Bell_Parameter_Fixed_Base(MatrixXd Bell_Parameter_Fixed_Base, int num_steps){

	std::ofstream output_Bell_Parameter_Fixed_Base("Data_Analysis/Bell_Parameter/Bell_Fixed_Base/Bell_Parameter_U_Fixed_Base.txt");

	for(int i = 0; i < num_steps; i++){
		for (int j = 0; j < 2; j++){
			if (j == 1){
				output_Bell_Parameter_Fixed_Base << Bell_Parameter_Fixed_Base(i,j) << std::endl;
			}
			else{
				output_Bell_Parameter_Fixed_Base << Bell_Parameter_Fixed_Base(i,j) << "\t";
			}
		}
	}
}

void Unitary::Save_txt_files_Energy(MatrixXcd G, int num_steps, int N1){

	std::ofstream output_G("Data_Analysis/Energy/Energy_Channel/G_U_Gamma_N"+to_string(N1)+".txt");

	for(int i = 0; i < num_steps; i++){
		for (int j = 0; j < 61; j++){
			if (j == 60){
				output_G << G(i,j).real() << std::endl;
			}
			else{
				output_G << G(i,j).real() << "\t";
			}
		}
	}	
}

void Unitary::Save_txt_files_Energy_Gamma(MatrixXcd G, int num_steps, int N1, int gamma_idx){

	std::ofstream output_G("Data_Analysis/Energy/Energy_Gamma/G_U_Gamma_N"+to_string(N1)+"_"+to_string(gamma_idx)+".txt");

	for(int i = 0; i < num_steps; i++){
		for (int j = 0; j < 61; j++){
			if (j == 60){
				output_G << G(i,j).real() << std::endl;
			}
			else{
				output_G << G(i,j).real() << "\t";
			}
		}
	}	
}

