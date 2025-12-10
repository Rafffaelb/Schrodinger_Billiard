#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <chrono>
#include <iomanip>  // 用于格式化输出
#include "../include/WignerDyson_AbstractClass_h/WignerDyson.h"
#include "../include/Quantum_chaotic_billiard.h"

using namespace std;
using namespace Eigen;

void WignerDyson::Run_Simulation_Conductance_Gamma(){

	auto start = chrono::system_clock::now();

	double Gamma, y, V;
       	int ress, N1, N2, n, _num_steps;

	ress = 100;
	_num_steps = 1000000;
	V = _lambda*_lambda/ress;

	for (int i = 1; i < 5; i++){
		if (i == 4){
			N1 = 10;
		}
		else{
			N1 = i;
		}

		N2 = N1;
		n = N1 + N2;

		// Create_ProjectionMatrices //
		MatrixXcd C1(_spin_deg * 2*N1, _spin_deg * 2*N1); MatrixXcd C2(_spin_deg * 2*N2, _spin_deg * 2*N2);
		C1.setZero(); C2.setZero();
		MatrixXcd *C1_pointer = &C1; MatrixXcd *C2_pointer = &C2;
		Create_ProjectionMatrices(C1_pointer, C2_pointer, N1, N2);

		// Gamma Vector //
		for (int gamma_idx = 1; gamma_idx < 22; gamma_idx++){
			if (gamma_idx == 1) {
				Gamma = 0.0001;
			} else {
				Gamma = double(gamma_idx - 1) / double(20);
			}
			
			cout << "\n===== GAMMA SIMULATION =====" << endl;
			cout << "Processing Gamma index " << gamma_idx << " of 21" << endl;
			cout << "Gamma value: " << Gamma << endl;
			cout << "Number of steps: " << _num_steps << endl;
			cout << "Channel count N: " << N1 << endl;
			cout << "===========================" << endl;
			
			auto gamma_start = chrono::system_clock::now();

			y = sqrt(double(1.0)/Gamma)*(1.0-sqrt(1.0-Gamma));

			// Create W Matrices //
			MatrixXcd W(_spin_deg * ress, _spin_deg * n);
			W.setZero();
			MatrixXcd *W_pointer = &W;
			Create_W(W_pointer, ress, N1, N2, _lambda, y);

			MatrixXcd G(_num_steps, 21);
			MatrixXcd P(_num_steps, 21);

			G.setZero();
			P.setZero();

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

				// Conductance (G) and Power Shot Noise (P) //
				billiard_setup.Calculate_G_and_P();

				G(step-1, gamma_idx-1) = billiard_setup.getG();
				P(step-1, gamma_idx-1) = billiard_setup.getP();
				
				// Improved progress reporting
				if (step % 500000 == 0 || step == _num_steps){
					auto current_time = chrono::system_clock::now();
					auto elapsed = chrono::duration_cast<chrono::milliseconds>(current_time - gamma_start);
					double progress = (double)step / _num_steps * 100;
					long long eta_ms = (long long)(elapsed.count() * (_num_steps - step) / (double)step);
					auto eta = chrono::milliseconds(eta_ms);
					
					int hours = chrono::duration_cast<chrono::hours>(eta).count();
					int minutes = chrono::duration_cast<chrono::minutes>(eta % chrono::hours(1)).count();
					int seconds = chrono::duration_cast<chrono::seconds>(eta % chrono::minutes(1)).count();
					
					#pragma omp critical
					{
						cout << fixed << setprecision(1);
						cout << "\rN=" << N1 << " | Gamma " << gamma_idx << "/21 | ";
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
			
			// Print completion for this gamma value
			auto gamma_end = chrono::system_clock::now();
			auto gamma_elapsed = chrono::duration_cast<chrono::milliseconds>(gamma_end - gamma_start);
			cout << "\nGamma " << Gamma << " (index " << gamma_idx << ") with N=" << N1 << " completed in " << gamma_elapsed.count()/1000.0 << " seconds" << endl;
	
			//Save G and P matrices as txt files //
			Save_txt_files_Gamma(G, P, _num_steps, N1);
		}
	}

	auto end = chrono::system_clock::now();
	auto elapsed =
		chrono::duration_cast<chrono::minutes>(end-start);
	cout << "\n===== SIMULATION COMPLETE =====" << endl;
	cout << "Total simulation time: " << elapsed.count() << " minutes" << endl;
	cout << "===============================" << endl;
}
