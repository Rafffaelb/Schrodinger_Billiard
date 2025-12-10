#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <chrono>
#include <iomanip>
#include "../include/WignerDyson_AbstractClass_h/WignerDyson.h"
#include "../include/Quantum_chaotic_billiard.h"

using namespace std;
using namespace Eigen;

void WignerDyson::Run_Simulation_Bell_Parameter_Ress(){
	
	auto start = chrono::system_clock::now();

	int ress, _num_steps, N1, N2, n;

	ress = 50;
	_num_steps = 100000;

	N1 = 1;
	N2 = N1;
	n = N1 + N2;

	double y = 0;
	double V = _lambda*_lambda/ress;

	MatrixXd Bell_Parameter_Ress(_num_steps, 5);

	Bell_Parameter_Ress.setZero();

	for (int ress_idx = 1; ress_idx < 6; ress_idx++){
		
		cout << "\n===== BELL PARAMETER (RESS) SIMULATION =====" << endl;
		cout << "Processing Ress index " << ress_idx << " of 5" << endl;
		cout << "Ress value: " << ress << endl;
		cout << "Number of steps per Ress: " << _num_steps << endl;
		cout << "==========================================" << endl;
		
		auto ress_start = chrono::system_clock::now();

		// Create W Matrices //

		MatrixXcd W(_spin_deg * ress, _spin_deg * n);
		W.setZero();
		MatrixXcd* W_pointer = &W;

		Create_W(W_pointer, ress, N1, N2, _lambda, y);

		// Create_ProjectionMatrices //

		MatrixXcd C1(_spin_deg * 2*N1, _spin_deg * 2*N1); MatrixXcd C2(_spin_deg * 2*N2, _spin_deg * 2*N2);
		C1.setZero(); C2.setZero();
		MatrixXcd *C1_pointer = &C1; MatrixXcd *C2_pointer = &C2;

		Create_ProjectionMatrices(C1_pointer, C2_pointer, N1, N2);

		#pragma omp parallel for shared(W, C1, C2)
		for (int step = 1; step < _num_steps + 1; step++){

			// Generate Hamiltonian Matrix //

			MatrixXcd H(_spin_deg * ress, _spin_deg * ress);
			H.setZero();
			MatrixXcd* H_pointer = &H;

			Create_H(H_pointer, ress, V);

			// Create billiard setup //

			Quantum_chaotic_billiard billiard_setup(H, W, C1, C2);

			// Scattering Matrix //
		
			billiard_setup.Calculate_Smatrix(0);

			// Bell Parameter //
		
			billiard_setup.Calculate_Bell_Parameter();

			Bell_Parameter_Ress(step-1, ress_idx-1) = billiard_setup.getBell_Parameter();

			// Improved progress reporting
			if (step % 10000 == 0 || step == _num_steps){
				auto current_time = chrono::system_clock::now();
				auto elapsed = chrono::duration_cast<chrono::milliseconds>(current_time - ress_start);
				double progress = (double)step / _num_steps * 100;
				long long eta_ms = (long long)(elapsed.count() * (_num_steps - step) / (double)step);
				auto eta = chrono::milliseconds(eta_ms);
				
				int hours = chrono::duration_cast<chrono::hours>(eta).count();
				int minutes = chrono::duration_cast<chrono::minutes>(eta % chrono::hours(1)).count();
				int seconds = chrono::duration_cast<chrono::seconds>(eta % chrono::minutes(1)).count();
				
				#pragma omp critical
				{
					cout << fixed << setprecision(1);
					cout << "\rRess " << ress_idx << "/5 | ";
					cout << "Step: " << step << "/" << _num_steps << " (" << progress << "%) | ";
					cout << "Elapsed: " << elapsed.count()/1000.0 << "s | ";
					cout << "ETA: " << setw(2) << setfill('0') << hours << ":"
					     << setw(2) << setfill('0') << minutes << ":"
					     << setw(2) << setfill('0') << seconds;
					if (step == _num_steps) cout << " [COMPLETE]";
					cout.flush();
				}
			}
		}
		
		// Print completion for this ress value
		auto ress_end = chrono::system_clock::now();
		auto ress_elapsed = chrono::duration_cast<chrono::milliseconds>(ress_end - ress_start);
		cout << "\nRess " << ress << " (index " << ress_idx << ") completed in " << ress_elapsed.count()/1000.0 << " seconds" << endl;
		
		//Save Concurrence matrix as txt files //
		Save_txt_files_Bell_Parameter_Ress(Bell_Parameter_Ress, _num_steps);

		ress = ress + 25;
	}

	auto end = chrono::system_clock::now();
	auto elapsed =
		chrono::duration_cast<chrono::minutes>(end-start);
	cout << "\n===== SIMULATION COMPLETE =====" << endl;
	cout << "Total simulation time: " << elapsed.count() << " minutes" << endl;
	cout << "===============================" << endl;
}