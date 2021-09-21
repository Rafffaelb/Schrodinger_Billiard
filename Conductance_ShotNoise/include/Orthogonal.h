#ifndef ORTHOGONAL_H
#define ORTHOGONAL_H

#include "WignerDyson.h"

class Orthogonal: public WignerDyson{

	public:
		Orthogonal(double Gamma, int ress, double lambda, int num_steps, int spin_deg);

		~Orthogonal();

		void Create_W(MatrixXcd* W_pointer, int _ress, int N1, int N2, double _lambda, double _y);

		void Create_ProjectionMatrices(MatrixXcd *C1_pointer, MatrixXcd* C2_pointer, int N1, int N2);

		void Create_H(MatrixXcd* H_pointer, int ress, double V);
		void Save_txt_files(MatrixXcd G, MatrixXcd P, int num_steps);
};

#endif