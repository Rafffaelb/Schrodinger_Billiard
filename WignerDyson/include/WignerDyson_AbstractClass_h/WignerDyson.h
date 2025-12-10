#include <Eigen/Dense>
#ifndef WIGNERDYSON_H
#define WIGNERDYSON_H

using namespace Eigen;

class WignerDyson{

	protected:

		int _spin_deg;
		double _lambda;

	public:

		WignerDyson();

		virtual ~WignerDyson();

		virtual void Create_W(MatrixXcd* W_pointer, int _ress, int N1, int N2, double _lambda, double y) = 0;
		
		void Run_Simulation_Conductance_Channels();
		void Run_Simulation_Conductance_Gamma();
		void Run_Simulation_Concurrence_Gamma();
		void Run_Simulation_Bell_Parameter_Ress();
		void Run_Simulation_Bell_Parameter_Gamma();
		void Run_Simulation_Bell_Parameter_Fixed_Base();
		void Run_Simulation_Conductance_Energy();
		void Run_Simulation_Conductance_Energy_Gamma();

		virtual void Create_ProjectionMatrices(MatrixXcd *C1_pointer, MatrixXcd* C2_pointer, int N1, int N2) = 0;

		virtual void Create_H(MatrixXcd* H_pointer, int _ress, double V) = 0;
	
		virtual void Save_txt_files_Channels(MatrixXcd G, MatrixXcd P, int num_steps) = 0;
		virtual void Save_txt_files_Gamma(MatrixXcd G, MatrixXcd P, int num_steps, int N1) = 0;
		virtual void Save_txt_files_Concurrence_Gamma(MatrixXd Concurrence, MatrixXd Entanglement, int num_steps) = 0;
		virtual void Save_txt_files_Bell_Parameter_Ress(MatrixXd Bell_Parameter_Ress, int num_steps) = 0;
		virtual void Save_txt_files_Bell_Parameter_Gamma(MatrixXd Bell_Parameter_Gamma, MatrixXd Bell_Parameter_Dephase_Gamma, int num_steps) = 0;
		virtual void Save_txt_files_Bell_Parameter_Fixed_Base(MatrixXd Bell_Parameter_Fixed_Base, int num_steps) = 0;
		virtual void Save_txt_files_Energy(MatrixXcd G, int num_steps, int N1) = 0;
		virtual void Save_txt_files_Energy_Gamma(MatrixXcd G, int num_steps, int N1, int gamma_idx) = 0;

};

#endif