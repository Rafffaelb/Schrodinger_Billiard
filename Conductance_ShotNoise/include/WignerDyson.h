#include <eigen3/Eigen/Dense>
#ifndef WIGNERDYSON_H
#define WIGNERDYSON_H

using namespace Eigen;

class WignerDyson{

	protected:

		int _ress;
		int _num_steps;
		int _spin_deg;
		double _Gamma;
		double _lambda;


	public:
		WignerDyson();

		virtual ~WignerDyson();

		virtual void Create_W(MatrixXcd* W_pointer, int _ress, int N1, int N2, double _lambda, double y) = 0;
		
		void Run_Simulation();

		virtual void Create_ProjectionMatrices(MatrixXcd *C1_pointer, MatrixXcd* C2_pointer, int N1, int N2) = 0;

		virtual void Create_H(MatrixXcd* H_pointer, int _ress, double V) = 0;
		virtual void Save_txt_files(MatrixXcd G, MatrixXcd P, int num_setps) = 0;
};

#endif