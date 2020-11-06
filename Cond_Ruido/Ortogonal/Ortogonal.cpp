#include <iostream>
#include <cmath>
#include <complex>
#include <random>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <ctime>
#include <fstream>
#include <chrono>

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

void Criando_W (MatrixXcd *W_pointer, int ress, int N1, int N2, double lambda, double y){

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

void Criando_Projetores (MatrixXcd *C1_pointer, MatrixXcd *C2_pointer, int N1, int N2){

	MatrixXcd C1tio(2,2);
	MatrixXcd C2tio(2,2);

	MatrixXcd identidade1 = MatrixXcd::Identity(N1,N1);
	MatrixXcd identidade2 = MatrixXcd::Identity(N2,N2);
	
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
			C1.block((i-1)*identidade1.rows(), (j-1)*identidade1.cols(), identidade1.rows(), identidade1.cols()) = C1tio(i-1,j-1)*identidade1;
			
		}
	}

	for (int i = 1; i < 3; i++){
		for (int j = 1; j < 3; j++){
			C2.block((i-1)*identidade2.rows(), (j-1)*identidade2.cols(), identidade2.rows(), identidade2.cols()) = C2tio(i-1,j-1)*identidade2;
		}
	}

	*C1_pointer << C1;
	*C2_pointer << C2;

}

int main(){

	std::complex<double> complex_identity(0, 1);
	std::complex<double> number_2(2, 0);

	// Input //

	double Gamma, lambda, y, V, gama;
	int N1, N2, n, ress, num_realization, lambda1;

	Gamma = 1;
	N1 = 1;
        N2 = 1;
    	n = N1+N2;
	ress = 100;
	lambda = 0.5;
	lambda1 = 1;
	y = sqrt(1.0/Gamma)*(1.0-sqrt(1.0-Gamma));
	V = lambda*lambda/ress;
	num_realization = 50000;
	
	MatrixXcd identityS(n,n);
	MatrixXcd G(num_realization,1);
	MatrixXcd R(num_realization,1);
	// Pauli Matrices //

	MatrixXcd matrizpauli1(2,2);
	MatrixXcd matrizpauli2(2,2);
	MatrixXcd matrizpauli3(2,2);

	matrizpauli1.real() << 1, 0, 0, 1;
	matrizpauli1.imag() <<  0, 0, 0,  0;

	matrizpauli2.real() << 0, 0, 0, 0;
	matrizpauli2.imag() <<  0, -1, 1,  0;

	
	matrizpauli3.real() << 1, 0, 0, -1;
	matrizpauli3.imag() <<  0, 0, 0,  0;

	// Creating W Matrices //

	MatrixXcd W(ress,n);
	W.setZero();
	MatrixXcd* W_pointer = &W;
	
	Criando_W(W_pointer, ress, N1, N2, lambda, y);

	// Creating Projectors //

	MatrixXcd C1(2*N1, 2*N1);
	MatrixXcd C2(2*N2, 2*N2);

	C1.setZero();
	C2.setZero();

	MatrixXcd *C1_pointer = &C1;
	MatrixXcd *C2_pointer = &C2;

	Criando_Projetores(C1_pointer, C2_pointer, N1, N2);
	
	for (int realization = 1; realization < num_realization + 1; realization++){

		// Generating Hamiltonian Matrix //

		MatrixXcd H(ress,ress);
		H.setZero();
		MatrixXcd* H_pointer = &H;
		Criando_H(H_pointer, ress, V);

		// Inverse Green Function //

		MatrixXcd D = (-H + complex_identity*M_PI*W*(W.adjoint()));
		
		// Scattering Matrix //

		MatrixXcd S = identityS -number_2*complex_identity*M_PI*(W.adjoint())*(D.inverse())*W;

		MatrixXcd ttdaga = C1*S*C2*(S.adjoint());

		// Conductance and Power Shot Noise //

		MatrixXcd identityR = MatrixXcd::Identity(ttdaga.rows(),ttdaga.cols());

		G(realization-1, 0) = ttdaga.trace();
		R(realization-1, 0) = (ttdaga*(identityR-ttdaga)).trace();

		if (realization % 10000 == 0){
			std::cout << realization << std::endl;
		}
	}

	std::ofstream output_G("G_O.txt");
	std::ofstream output_R("R_O.txt");
	for (int i = 0; i < num_realization; i++){
		for (int j = 0; j < 1; j++){
			output_G << G(i,j).real() << " " << std::endl;
			output_R << R(i,j).real() << " " << std::endl;
		}
	}
	return 0;
}
