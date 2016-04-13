/* Jin Xu, a graduate student */
/* Department of Physics, University of Alberta, Edmonton, Alberta, T6G 2E1, Canada */
/* PHYS 580, the term project */

/* This library is for matrix write and matrix read */

#include <Eigen/Dense>
using namespace Eigen;

#include <vector>
using namespace std;

void save(const char *filename, const MatrixXd& m);

void load(const char *filename, MatrixXd& m);

void infinite_write(vector<MatrixXd> correlation, int system_size);

void finite_write(vector<MatrixXd> correlation, int system_size, int iter);

void finite_read(vector<MatrixXd>& correlation, int system_size, int iter);

void printf_bonds(vector<MatrixXd> nearest_correlation, VectorXd evec, int iter, int system_size);