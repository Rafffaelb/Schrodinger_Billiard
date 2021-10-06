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

void WignerDyson::Save_txt_files_Concurrence_Gamma(MatrixXd, int num_steps) {};
