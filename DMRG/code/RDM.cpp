/* Jin Xu, a graduate student */
/* Department of Physics, University of Alberta, Edmonton, Alberta, T6G 2E1, Canada */
/* PHYS 580, the term project */

/* This library is for the reduced density matrix */

#include <Eigen/Dense>
using namespace Eigen;

#include <cstdlib>
using namespace std;

#include <iostream>
using namespace std;

#include "tensor.hpp"
#include "ED.hpp"

MatrixXd reduced_density_matrix(VectorXd evec, int n)
{
	// The size of the superblock is 2n
	// The vector evec represents the ground state of the superblock
	// The size of the vector evec has to be power(2,2n)
	// The function will return a reduced density matrix whose size is power(2,n)
	
	int L=power(2,n);
	if(evec.size()!=L*L)
	{
		cout << "The size of the vector representing the ground state of the superblock should be power(2,2n)! " << endl;
		exit(0);
	}
	
	MatrixXd density(L, L);
	int i, j, k;
	for(i=0;i<L;i++)
		for(j=0;j<L;j++)
		{
			density(i,j)=0;
			for(k=0;k<L;k++)
				density(i,j)=density(i,j)+evec(i*L+k)*evec(j*L+k);
		}
	
	return density;
}

MatrixXd transformation_matrix(MatrixXd density, int m)
{
	// The dimension of the matrix density has to be power(2,m+1)xpower(2,m+1)
	// The function will return a power(2,m+1)xpower(2,m) transformation matrix T
	
	int L=power(2,m);
	if(density.rows()!=2*L)
	{
		cout << "density.rows() should be equal to power(2,m+1)!" << endl;
		exit(0);
	}
	if(density.cols()!=2*L)
	{
		cout << "density.cols() should be equal to power(2,m+1)!" << endl;
		exit(0);
	}

	VectorXd eval(2*L);
	MatrixXd evec(2*L,2*L);
	diagonalization(density, eval, evec);
	
	MatrixXd T(2*L,L);
	T=evec.block(0, L, 2*L, L);
	
	return T;
}
	
	
	
	
	
	
	

