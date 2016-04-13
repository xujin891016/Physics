/* Jin Xu, a graduate student */
/* Department of Physics, University of Alberta, Edmonton, Alberta, T6G 2E1, Canada */
/* PHYS 580, the term project */

/* This library is for the exact diagonalization */

#include <Eigen/Dense>
using namespace Eigen;

int index_min(VectorXd A)
{
	// Find the index for the minimum value in the vector
	int index;
	int i;
	double mini;
	
	index=0;
	mini=A(0);
	for(i=1; i<A.size(); i++)
	{
		if(A(i)<mini)
		{
			index=i;
			mini=A(i);
		}
	}
	
	return index;
}

VectorXi sort_vector(VectorXd A)
{
	// Sort a vector by returning the order of indices
	double replacement=A.maxCoeff()+0.5;
	int i;
	VectorXi index(A.size());
	
	for(i=0;i<A.size();i++)
	{
		index(i)=index_min(A);
		A(index(i))=replacement;
	}
	
	return index;
}

double groundstate(MatrixXd H, VectorXd& evec)
{
	// Diagonalize the Hamiltonian matrix and get the ground state
	SelfAdjointEigenSolver<MatrixXd> eigensolver(H);
	if(eigensolver.info()!=Success)
		abort();
	MatrixXd states=eigensolver.eigenvectors();
	VectorXd values=eigensolver.eigenvalues();
	int index=index_min(values);
	double groundstate_energy=values(index);
	evec=states.col(index);
	return groundstate_energy;
}

void diagonalization(MatrixXd D, VectorXd& eval, MatrixXd& evec)
{
	// Diagonalize the density matrix and sort them
	SelfAdjointEigenSolver<MatrixXd> eigensolver(D);
	if(eigensolver.info()!=Success)
		abort();
	MatrixXd states=eigensolver.eigenvectors();
	VectorXd values=eigensolver.eigenvalues();
	VectorXi index=sort_vector(values);
	for(int i=0; i<index.size(); i++)
	{
		eval(i)=values(index(i));
		evec.col(i)=states.col(index(i));
	}
}