#include <iostream>
#ifndef QUANTUM_CHAOTIC_BILLIARD_O_H
#define QUANTUM_CHAOTIC_BILLIARD_O_H

using namespace std;
using namespace Eigen;

class Quantum_chaotic_billiard
{
	private:
		MatrixXcd _H;
		MatrixXcd _W;
		MatrixXcd _C1;
		MatrixXcd _C2;
	
	public:
		Quantum_chaotic_billiard(MatrixXcd H, MatrixXcd W, MatrixXcd C1, MatrixXcd C2);

		void Set_Setup(MatrixXcd H, MatrixXcd W, MatrixXcd C1, MatrixXcd C2);

		void Calculate_Smatrix();
};

#endif
