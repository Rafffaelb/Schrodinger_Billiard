#include <iostream>
#include <chrono>
#include "../../include/WignerDyson_AbstractClass_h/WignerDyson.h"
#include "../../include/Quantum_chaotic_billiard.h"
#include <eigen3/Eigen/Dense>
#include <omp.h>

using namespace std;
using namespace Eigen;

WignerDyson::WignerDyson() {};

WignerDyson::~WignerDyson() {};

void WignerDyson::Create_W(MatrixXcd* W_pointer, int _ress, int N1, int N2, double _lambda, double y) {};

void WignerDyson::Create_ProjectionMatrices(MatrixXcd *C1_pointer, MatrixXcd* C2_pointer, int N1, int N2) {};

void WignerDyson::Create_H(MatrixXcd* H_pointer, int _ress, double V) {};

void WignerDyson::Save_txt_files_Channels(MatrixXcd G, MatrixXcd P, int num_steps) {};

void WignerDyson::Save_txt_files_Gamma(MatrixXcd G, MatrixXcd P, int num_steps, int N1) {};

void WignerDyson::Save_txt_files_Concurrence_Gamma(MatrixXd Concurrence, MatrixXd Entanglement, int num_steps) {};

void WignerDyson::Save_txt_files_Bell_Parameter_Ress(MatrixXd Bell_Parameter_Ress, int num_steps) {};

void WignerDyson::Save_txt_files_Bell_Parameter_Gamma(MatrixXd Bell_Parameter_Gamma, MatrixXd Bell_Parameter_Dephase_Gamma, int num_steps) {};

void WignerDyson::Save_txt_files_Bell_Parameter_Fixed_Base(MatrixXd Bell_Parameter_Fixed_Base, int num_steps) {};

void WignerDyson::Save_txt_files_Energy(MatrixXcd G, int num_steps, int N1) {};

void WignerDyson::Save_txt_files_Energy_Gamma(MatrixXcd G, int num_steps, int N1, int gamma_idx) {};

