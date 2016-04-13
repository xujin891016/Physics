/* Jin Xu, a graduate student */
/* Department of Physics, University of Alberta, Edmonton, Alberta, T6G 2E1, Canada */
/* PHYS 580, the term project */

/* This library is for the exact diagonalization */

#include <Eigen/Dense>
using namespace Eigen;

int index_min(VectorXd A);

VectorXi sort_vector(VectorXd A);

double groundstate(MatrixXd H, VectorXd& evec);

void diagonalization(MatrixXd D, VectorXd& eval, MatrixXd& evec);