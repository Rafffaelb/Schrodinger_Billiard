#include <iostream>
#ifndef PROJECTIONMATRICES_H
#define PROJECTIONMATRICES_H

using namespace std;
using namespace Eigen;

void Create_ProjectionMatrices (MatrixXcd *C1_pointer, MatrixXcd *C2_pointer, int N1, int N2);

void Create_ProjectionMatrices_S (MatrixXcd *C1_pointer, MatrixXcd *C2_pointer, int N1, int N2);

#endif
