#ifndef W_H
#define W_H

using namespace std;
using namespace Eigen;

void Create_W(MatrixXcd *W_pointer, int ress, int N1, int N2, double lambda, double y);

void Create_W_S(MatrixXcd *W_pointer, int ress, int N1, int N2, double lambda, double y);

#endif
