#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <fstream>

using namespace std;
using namespace Eigen;

void Save_txt_files(MatrixXcd G, MatrixXcd P, int num_steps){

	std::ofstream output_G("G.txt");
	std::ofstream output_P("P.txt");

	for(int i = 0; i < num_steps; i++){
		for (int j = 0; j < 10; j++){
			if (j == 9){
				output_G << G(i,j).real() << std::endl;
				output_P << P(i,j).real() << std::endl;
			}
			else{
				output_G << G(i,j).real() << "\t";
				output_P << P(i,j).real() << "\t";
			}
		}
	}
}
