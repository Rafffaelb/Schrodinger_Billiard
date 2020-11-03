#include <iostream>
#include <cmath>
#include <complex>
#include <random>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <ctime>
#include <fstream>


using namespace Eigen;
using namespace std::literals;

int main(){

	std::complex<double> complex_identity(0, 1);
	std::complex<double> number_2(2, 0);
	std::default_random_engine generator1;
  	std::normal_distribution<double> distribution1(0.0,1.0);

	double Gamma, lambda, y, V, gama;
	int N1, N2, n, ress, num_realization, lambda1;

	Gamma = 1;
	N1 = 1;
    	N2 = 1;
    	n = N1+N2;
	ress = 50;
	lambda = 0.5;
	lambda1 = 1;
	y = sqrt(1.0/Gamma)*(1.0-sqrt(1.0-Gamma));
	V = lambda*lambda/ress;
	num_realization = 100000;

	MatrixXcd G(num_realization,1);
	MatrixXcd R(num_realization,1);
	MatrixXcd Teste(num_realization,1);

	// Pauli Matrices //

    	MatrixXcd identidade2x2 = MatrixXcd::Identity(2,2);
    	MatrixXcd identidade2N1 = MatrixXcd::Identity(2*N1,2*N1);
	MatrixXcd identidade2N2 = MatrixXcd::Identity(2*N2,2*N2);
   	MatrixXcd identidade_canal = MatrixXcd::Identity(n,n);

	MatrixXcd matrizpauli1(2,2);
	MatrixXcd matrizpauli2(2,2);
	MatrixXcd matrizpauli3(2,2);
    	MatrixXcd Sigma_y(2*n,2*n);

	matrizpauli1.real() << 0, 1, 1, 0;
	matrizpauli1.imag() <<  0, 0, 0,  0;

	matrizpauli2.real() << 0, 0, 0, 0;
	matrizpauli2.imag() <<  0, -1, 1,  0;

	
	matrizpauli3.real() << 1, 0, 0, -1;
	matrizpauli3.imag() <<  0, 0, 0,  0;

    	for (int i = 1; i < identidade_canal.rows() + 1; i++){
        	for (int j = 1; j < identidade_canal.cols() + 1; j++){
            		Sigma_y.block((i-1)*matrizpauli2.rows(),(j-1)*matrizpauli2.cols(), matrizpauli2.rows(), matrizpauli2.cols()) = identidade_canal(i-1,j-1)*matrizpauli2;
        	}
    	}
    	// Creating Projectors //

	Eigen::MatrixXcd C1tio(2,2);
	Eigen::MatrixXcd C2tio(2,2);

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

	MatrixXcd C1(identidade2N1.rows() * C1tio.rows(),identidade2N1.cols() * C1tio.cols() );
	MatrixXcd C2(identidade2N2.rows() * C2tio.rows(),identidade2N2.cols() * C2tio.cols() );

	C1.setZero();
	C2.setZero();

	for (int i = 1; i < 3; i++){
		for (int j = 1; j < 3; j++){
			C1.block((i-1)*identidade2N1.rows(), (j-1)*identidade2N1.cols(), identidade2N1.rows(), identidade2N1.cols()) = C1tio(i-1,j-1)*identidade2N1;
		}
	}

	for (int i = 1; i < 3; i++){
		for (int j = 1; j < 3; j++){
			C2.block((i-1)*identidade2N2.rows(), (j-1)*identidade2N2.cols(), identidade2N2.rows(), identidade2N2.cols()) = C2tio(i-1,j-1)*identidade2N2;
		}
	}

	// Creating W Matrices //

	Eigen::MatrixXcd W1(ress,N1);
	Eigen::MatrixXcd W2(ress,N2);
	Eigen::MatrixXcd W(W1.rows(), W1.cols() + W2.cols());
	Eigen::MatrixXcd W_simpl(2*W1.rows(),2*(W1.cols() + W2.cols()));

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

    	for (int i = 1; i < W.rows()+1; i++){
        	for (int j = 1; j < W.cols()+1; j++){
            	W_simpl.block((i-1)*identidade2x2.rows(), (j-1)*identidade2x2.cols(), identidade2x2.rows(), identidade2x2.cols()) = W(i-1,j-1)*identidade2x2;
        	}
    	}

	MatrixXcd A_real(ress,ress); MatrixXcd A_imag(ress,ress);
	MatrixXcd C_real(ress,ress); MatrixXcd C_imag(ress,ress);
    	MatrixXcd H(ress,ress); MatrixXcd H1(ress,ress);
    	MatrixXcd H2(ress,ress); MatrixXcd H3(ress,ress);
    	MatrixXcd Q0(2*ress,2*ress); MatrixXcd Q1(2*ress,2*ress);
    	MatrixXcd Q2(2*ress,2*ress); MatrixXcd Q3(2*ress,2*ress);
	MatrixXcd Simetrica(ress,ress); MatrixXcd Antisimetrica1(ress,ress);
	MatrixXcd Antisimetrica2(ress,ress); MatrixXcd Antisimetrica3(ress,ress);
    	MatrixXcd identityS(W_simpl.cols(),W_simpl.cols());
	MatrixXcd D(2*ress,2*ress); MatrixXcd Ham(2*ress,2*ress);
	MatrixXcd S(W_simpl.cols(),W_simpl.cols());

	for (int realization = 1; realization < num_realization + 1; realization++){

		// Generating Hamiltonian Matrix //
		A_real.setZero();
		A_imag.setZero();
		C_real.setZero();
		C_imag.setZero();

		for (int i = 1; i < ress + 1; i++){
			for (int j = 1; j < ress + 1; j++){
				std::complex<double> aux_1 = distribution1(generator1);
				A_real(i-1,j-1) = aux_1;
			}
		}

		for (int i = 1; i < ress + 1; i++){
			for (int j = 1; j < ress + 1; j++){
				std::complex<double> aux_1 = distribution1(generator1);
				C_real(i-1,j-1) = aux_1;
			}
		}

		MatrixXcd A = A_real;
		MatrixXcd C = C_real;

		for (int i = 1; i < ress + 1; i++){
			H(i-1,i-1) = A(i-1,i-1)*sqrt(V/(2.0));
			H1(i-1,i-1) = 0.0;
            		H2(i-1,i-1) = 0.0;
            		H3(i-1,i-1) = 0.0;
			for (int j = i + 1; j < ress + 1; j++){
				H(i-1,j-1) = A(i-1,j-1)*sqrt(V/(4.0));
				H1(i-1,j-1) = A(j-1,i-1)*sqrt(V/(4.0));
                		H2(i-1,j-1) = C(j-1,i-1)*sqrt(V/(4.0));
                		H3(i-1,j-1) = C(i-1,j-1)*sqrt(V/(4.0));
			}
		}

		Simetrica << H + H.adjoint(); Antisimetrica1 << H1 - H1.adjoint();
		Antisimetrica2 << H2 - H2.adjoint(); Antisimetrica3 << H3 - H3.adjoint();

		for (int i = 1; i < Simetrica.rows()+1; i++){
            		for (int j = 1; j < Simetrica.cols()+1; j++){
                		Q0.block((i-1)*identidade2x2.rows(), (j-1)*identidade2x2.cols(), identidade2x2.rows(), identidade2x2.cols()) = Simetrica(i-1,j-1)*identidade2x2;
                		Q1.block((i-1)*matrizpauli1.rows(), (j-1)*matrizpauli1.cols(), matrizpauli1.rows(), matrizpauli1.cols()) = Antisimetrica1(i-1,j-1)*matrizpauli1;
                		Q2.block((i-1)*matrizpauli2.rows(), (j-1)*matrizpauli2.cols(), matrizpauli2.rows(), matrizpauli2.cols()) = Antisimetrica2(i-1,j-1)*matrizpauli2;
                		Q3.block((i-1)*matrizpauli3.rows(), (j-1)*matrizpauli3.cols(), matrizpauli3.rows(), matrizpauli3.cols()) = Antisimetrica3(i-1,j-1)*matrizpauli3;
            		}
        	}

        	Ham << Q0 + complex_identity*(Q1+Q2+Q3);
		
		D << (-Ham + complex_identity*M_PI*W_simpl*(W_simpl.adjoint()));
		// Scattering Matrix //

		identityS << MatrixXcd::Identity(W_simpl.cols(),W_simpl.cols());
		S << identityS - number_2*complex_identity*M_PI*(W_simpl.adjoint())*(D.inverse())*W_simpl;

		MatrixXcd Teste = Sigma_y*(S.transpose())*Sigma_y;
        	
		MatrixXcd ttdaga = C1*S*C2*(S.adjoint());

		// Conductance and Power Shot Noise //

		MatrixXcd identityR = MatrixXcd::Identity(ttdaga.rows(),ttdaga.cols());

		G(realization-1, 0) = ttdaga.trace();
		R(realization-1, 0) = (ttdaga*(identityR-ttdaga)).trace();

		if (realization % 10000 == 0){
			std::cout << realization << std::endl;
		}
	}

	std::ofstream output_G("G_S.txt");
	std::ofstream output_R("R_S.txt");
	for (int i = 0; i < num_realization; i++){
		for (int j = 0; j < 1; j++){
			output_G << G(i,j).real() << " " << std::endl;
			output_R << R(i,j).real() << " " << std::endl;
		}
	}
	return 0;
}
