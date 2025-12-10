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

void WignerDyson::Run_Simulation_Conductance_Energy_Gamma(){
	
	auto start = chrono::system_clock::now();

	int ress, _num_steps, N1, N2, n;

	ress = 100;
	_num_steps = 10000;

	N1 = 2;
	N2 = N1;
	n = N1 + N2;

	// Gamma Vector //

	VectorXd Gamma(10);
	Gamma << 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0;

	// Loop over Gamma values
	for (int gamma_idx = 1; gamma_idx < 11; gamma_idx++){
		
		cout << "\n===== ENERGY-GAMMA SIMULATION =====" << endl;
		cout << "Processing Gamma index " << gamma_idx << " of 10" << endl;
		cout << "Gamma value: " << Gamma(gamma_idx-1) << endl;
		cout << "Number of steps per Gamma: " << _num_steps << endl;
		cout << "Energy points per step: 61" << endl;
		cout << "===============================" << endl;
		
		auto gamma_start = chrono::system_clock::now();

		double Delta = 2*M_PI/(sqrt(double(ress))*_lambda);
		double small_gamma;

		// Create W Matrices //

		MatrixXcd W(_spin_deg * ress, _spin_deg * n);
		W.setZero();
		MatrixXcd* W_pointer = &W;

		double y = sqrt(1.0/Gamma(gamma_idx-1))*(1.0-sqrt(1.0-Gamma(gamma_idx-1)));
		double V = _lambda*_lambda/ress;

		Create_W(W_pointer, ress, N1, N2, _lambda, y);

		// Create_ProjectionMatrices //

		MatrixXcd C1(_spin_deg * 2*N1, _spin_deg * 2*N1); MatrixXcd C2(_spin_deg * 2*N2, _spin_deg * 2*N2);
		C1.setZero(); C2.setZero();
		MatrixXcd *C1_pointer = &C1; MatrixXcd *C2_pointer = &C2;

		Create_ProjectionMatrices(C1_pointer, C2_pointer, N1, N2);

		MatrixXcd G(_num_steps, 61);
		G.setZero();

		double Energy;

		#pragma omp parallel for private(Energy) shared(W, C1, C2, G, small_gamma, Delta)
		for (int step = 1; step < _num_steps + 1; step++){

			// Generate Hamiltonian Matrix //

			MatrixXcd H(_spin_deg * ress, _spin_deg * ress);
			H.setZero();
			MatrixXcd* H_pointer = &H;

			Create_H(H_pointer, ress, V);

			// Create billiard setup //

			Quantum_chaotic_billiard billiard_setup(H, W, C1, C2);

			small_gamma = (gamma_idx < 5) ? (N1*Delta/(8*M_PI)) : (N1*Gamma(gamma_idx-1)*Delta/(8*M_PI));

			#pragma omp parallel for shared(W, C1, C2, H) firstprivate(billiard_setup)
			for (int energy_idx = 1; energy_idx < 62; energy_idx++){
				
				Energy = small_gamma*((double)energy_idx-31);
				
				// Scattering Matrix //
			
				billiard_setup.Calculate_Smatrix(Energy);

				// Conductance (G) and Power Shot Noise (P) //
	
				billiard_setup.Calculate_G_and_P();

				G(step-1, energy_idx-1) = billiard_setup.getG();
			}

			// Improved progress reporting
			if (step % 1000 == 0 || step == _num_steps){
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
					cout << "\rGamma " << gamma_idx << "/10 | ";
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
		cout << "\nGamma " << Gamma(gamma_idx-1) << " (index " << gamma_idx << ") completed in " << gamma_elapsed.count()/1000.0 << " seconds" << endl;
	
		//Save G matrix as txt files //
		Save_txt_files_Energy_Gamma(G, _num_steps, N1, gamma_idx);
	}

	auto end = chrono::system_clock::now();
	auto elapsed =
		chrono::duration_cast<chrono::minutes>(end-start);
	cout << "\n===== SIMULATION COMPLETE =====" << endl;
	cout << "Total simulation time: " << elapsed.count() << " minutes" << endl;
	cout << "===============================" << endl;
}