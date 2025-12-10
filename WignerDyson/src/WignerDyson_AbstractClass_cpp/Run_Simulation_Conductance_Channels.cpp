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

void WignerDyson::Run_Simulation_Conductance_Channels(){

	auto start = chrono::system_clock::now();

	double Gamma, y, V;
       	int ress, _num_steps;

	ress = 100;
	_num_steps = 100000;

	Gamma = 1;
	y = sqrt(1.0/Gamma)*(1.0-sqrt(1.0-Gamma));
	V = _lambda*_lambda/ress;

	MatrixXcd G(_num_steps, 10);
	MatrixXcd P(_num_steps, 10);

	G.setZero();
	P.setZero();

	for (int N1 = 1; N1 < 11; N1++){
		
		int N2 = N1;
		int n = N1 + N2;

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
		
		// Print channel information
		cout << "\n===== CHANNEL SIMULATION =====" << endl;
		cout << "Processing channel " << N1 << " of 10" << endl;
		cout << "Number of steps per channel: " << _num_steps << endl;
		cout << "=============================" << endl;
		
		auto channel_start = chrono::system_clock::now();
		
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

			G(step-1, N1-1) = billiard_setup.getG();
			P(step-1, N1-1) = billiard_setup.getP();
	
			// Improved progress reporting
			if (step % 10000 == 0 || step == _num_steps){
				auto current_time = chrono::system_clock::now();
				auto elapsed = chrono::duration_cast<chrono::milliseconds>(current_time - channel_start);
				double progress = (double)step / _num_steps * 100;
				long long eta_ms = (long long)(elapsed.count() * (_num_steps - step) / (double)step);
				auto eta = chrono::milliseconds(eta_ms);
				
				int hours = chrono::duration_cast<chrono::hours>(eta).count();
				int minutes = chrono::duration_cast<chrono::minutes>(eta % chrono::hours(1)).count();
				int seconds = chrono::duration_cast<chrono::seconds>(eta % chrono::minutes(1)).count();
				
				#pragma omp critical
				{
					cout << fixed << setprecision(1);
					cout << "\rChannel " << N1 << "/10 | ";
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
		
		// Print completion for this channel
		auto channel_end = chrono::system_clock::now();
		auto channel_elapsed = chrono::duration_cast<chrono::milliseconds>(channel_end - channel_start);
		cout << "\nChannel " << N1 << " completed in " << channel_elapsed.count()/1000.0 << " seconds" << endl;
	
		//Save G and P matrices as txt files //
	
		Save_txt_files_Channels(G, P, _num_steps);
	}

	auto end = chrono::system_clock::now();
	auto elapsed =
		chrono::duration_cast<chrono::minutes>(end-start);
	cout << "\n===== SIMULATION COMPLETE =====" << endl;
	cout << "Total simulation time: " << elapsed.count() << " minutes" << endl;
	cout << "===============================" << endl;
}