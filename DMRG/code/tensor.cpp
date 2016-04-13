/* Jin Xu, a graduate student */
/* Department of Physics, University of Alberta, Edmonton, Alberta, T6G 2E1, Canada */
/* PHYS 580, the term project */

/* This library is for tensor operations */

#include <Eigen/Dense>
using namespace Eigen;

int power(int x, int y)
{
	int i;
	int value=1;
	for(i=0;i<y;i++)
		value=value*x;
	return value;
} // Both x and y are positive integers

MatrixXd tensor_product(MatrixXd A, MatrixXd B)
{
	MatrixXd C(A.rows()*B.rows(), A.cols()*B.cols());
	int i, j;
	for(i=0; i<A.rows(); i++)
		for(j=0; j<A.cols(); j++)
			C.block(i*B.rows(), j*B.cols(), B.rows(), B.cols())=A(i,j)*B;
	return C;
} // The tensor product of two matrices
