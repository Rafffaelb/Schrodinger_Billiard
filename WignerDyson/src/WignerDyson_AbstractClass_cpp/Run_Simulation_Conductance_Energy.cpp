#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <chrono>
#include <iomanip>  // For formatted output
#include "../include/WignerDyson_AbstractClass_h/WignerDyson.h"
#include "../include/Quantum_chaotic_billiard.h"

using namespace std;
using namespace Eigen;

void WignerDyson::Run_Simulation_Conductance_Energy(){

	auto start = chrono::system_clock::now();

	double Gamma, Delta, y, V, Energy;
       	int N1, N2, n, _num_steps;

	Delta = 0.01;
	Gamma = 1;
	y = sqrt(double(1.0)/Gamma)*(1.0-sqrt(1.0-Gamma));
	int ress = 100;
	_num_steps = 5000; 
	V = _lambda*_lambda/ress;

	for (int i = 1; i < 11; i++){

		N1 = i;

		N2 = N1;
		n = N1 + N2;

		// Create_ProjectionMatrices //

		MatrixXcd C1(_spin_deg * 2*N1, _spin_deg * 2*N1); MatrixXcd C2(_spin_deg * 2*N2, _spin_deg * 2*N2);
		C1.setZero(); C2.setZero();
		MatrixXcd *C1_pointer = &C1; MatrixXcd *C2_pointer = &C2;

		Create_ProjectionMatrices(C1_pointer, C2_pointer, N1, N2);
	
		MatrixXcd G(_num_steps, 61);

		G.setZero();
		
		// Create W Matrices //

		MatrixXcd W(_spin_deg * ress, _spin_deg * n);
		W.setZero();
		MatrixXcd *W_pointer = &W;

		Create_W(W_pointer, ress, N1, N2, _lambda, y);
		
		auto energy_start = chrono::system_clock::now();
		
		// Report energy simulation start
		double energy_min = ((N1*Gamma*Delta/M_PI)*(((double)1-31)/2))/8;
		double energy_max = ((N1*Gamma*Delta/M_PI)*(((double)61-31)/2))/8;
		
		cout << "\n===== ENERGY SIMULATION =====" << endl;
		cout << "Processing open channels N = " << N1 << endl;
		cout << "Energy range: " << energy_min << " to " << energy_max << endl;
		cout << "Number of steps per energy: " << _num_steps << endl;
		cout << "=============================" << endl;

		for (int step = 1; step < _num_steps + 1; step++){
		
			// Generate Hamiltonian Matrix //

			MatrixXcd H(_spin_deg * ress, _spin_deg * ress);
			H.setZero();
			MatrixXcd* H_pointer = &H;

			Create_H(H_pointer, ress, V);

			// Create billiard setup //

			Quantum_chaotic_billiard billiard_setup(H, W, C1, C2);
			
			#pragma omp parallel for shared(H, W, C1, C2) firstprivate(billiard_setup)
			for (int energy_idx = 1; energy_idx < 62; energy_idx++){

				Energy = ((N1*Gamma*Delta/M_PI)*(((double)energy_idx-31)/2))/8;

				// Scattering Matrix //
				
				billiard_setup.Calculate_Smatrix(Energy);

				// Conductance (G) and Power Shot Noise (P) //
		
				billiard_setup.Calculate_G_and_P();

				G(step-1, energy_idx-1) = billiard_setup.getG();
			}

			// Progress reporting with percentage and ETA
			if (step % 100 == 0 || step == _num_steps){
				auto current_time = chrono::system_clock::now();
				auto elapsed = chrono::duration_cast<chrono::milliseconds>(current_time - energy_start);
				double progress = (double)step / _num_steps * 100;
				long long eta_ms = (long long)(elapsed.count() * (_num_steps - step) / (double)step);
				auto eta = chrono::milliseconds(eta_ms);
				
				int hours = chrono::duration_cast<chrono::hours>(eta).count();
				int minutes = chrono::duration_cast<chrono::minutes>(eta % chrono::hours(1)).count();
				int seconds = chrono::duration_cast<chrono::seconds>(eta % chrono::minutes(1)).count();
				
				#pragma omp critical
				{
					cout << fixed << setprecision(1);
					cout << "\rEnergy N=" << N1 << "/10 | ";
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
		
		// Print completion for this energy value
		auto energy_end = chrono::system_clock::now();
		auto energy_elapsed = chrono::duration_cast<chrono::milliseconds>(energy_end - energy_start);
		cout << "\n\nEnergy N=" << N1 << " completed in " << energy_elapsed.count()/1000.0 << " seconds" << endl;

		//Save G matrix as txt files //
		Save_txt_files_Energy(G, _num_steps, N1);

	}

	auto end = chrono::system_clock::now();
	auto elapsed =
		chrono::duration_cast<chrono::minutes>(end-start);
	cout << "\n Simulation time duration: " << elapsed.count() << "\n";
}
