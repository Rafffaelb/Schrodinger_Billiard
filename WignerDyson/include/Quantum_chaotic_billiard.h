#include <iostream>
#ifndef QUANTUM_CHAOTIC_BILLIARD_H
#define QUANTUM_CHAOTIC_BILLIARD_H

using namespace std;
using namespace Eigen;

class Quantum_chaotic_billiard
{
	private:
		complex<double> _G;
		complex<double> _P;
		double _Concurrence;
		double _Entanglement;
		double _Bell_Parameter;
		double _Bell_Parameter_Dephase;

		MatrixXcd _H;
		MatrixXcd _W;
		MatrixXcd _C1;
		MatrixXcd _C2;
		MatrixXcd _S;
	
	public:
		Quantum_chaotic_billiard(MatrixXcd H, MatrixXcd W, MatrixXcd C1, MatrixXcd C2);

		void Set_Setup(MatrixXcd H, MatrixXcd W, MatrixXcd C1, MatrixXcd C2);

		void Calculate_Smatrix();

		void Calculate_G_and_P();

		complex<double> getG();

		complex<double> getP();

		void Calculate_Concurrence();

		double getConcurrence();
		
		double getEntanglement();

		void Calculate_Bell_Parameter();

		void Calculate_Bell_Parameter_Fixed_Base();

		double getBell_Parameter();

		double getBell_Parameter_Dephase();
};

#endif
