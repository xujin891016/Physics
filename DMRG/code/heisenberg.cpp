/* Jin Xu, a graduate student */
/* Department of Physics, University of Alberta, Edmonton, Alberta, T6G 2E1, Canada */
/* PHYS 580, the term project */

/* In this program, DMRG is used to study the Heisenberg chain */
/* The interations are between the nearest neighbours, with open boundary conditions */
/* The physical quantities measured are the ground state energy and the bond strengths */
/* The external library Eigen is used: http://eigen.tuxfamily.org/dox/ */

#include <Eigen/Dense>
using namespace Eigen;

#include <vector>
using namespace std;

#include <iostream>
using namespace std;

#include <iomanip>
using namespace std;

#include <cstdlib>
using namespace std;

#include "tensor.hpp"
#include "ED.hpp"
#include "RDM.hpp"
#include "WR.hpp"

int main(int argc, char **argv)
{
	if(argc !=4)
	{
		 cerr << "There should be three arguments!" << endl;
		 exit(1);
	}
	
	int size_of_chain=atoi(argv[1]); // The first argument is the size of the chain
	int m=atoi(argv[2]); // power(2,m) states will be kept after each truncation
	int halfsweeps=atoi(argv[3]); // The number of halfsweeps in the finite system DMRG part
	
	/* Initialize the infinite system DMRG */
	MatrixXd sigma_z(2,2);
	sigma_z << 0.5, 0, 0, -0.5;
	MatrixXd sigma_p(2,2);
	sigma_p << 0, 1, 0, 0;
	MatrixXd sigma_n(2,2);
	sigma_n << 0, 0, 1, 0; 
	MatrixXd I2(2,2);
	I2=MatrixXd::Identity(2,2);// The matrices from a single spin
	
	MatrixXd single_bond(4,4); // The matrices for two interating spins
	single_bond=tensor_product(sigma_z, sigma_z)+0.5*tensor_product(sigma_p, sigma_n)+0.5*tensor_product(sigma_n, sigma_p);
	
	vector<MatrixXd> correlation(m);
	int i;
	for(i=0; i< m; i++)
	{
		MatrixXd left=MatrixXd::Identity(power(2,i),power(2,i));
		MatrixXd right=MatrixXd::Identity(power(2,m-i-1), power(2, m-i-1));
		correlation[i]=tensor_product(left, tensor_product(single_bond, right));
	} // The matrix vector for bonds
	
	MatrixXd I_m=MatrixXd::Identity(power(2,m),power(2,m));
	MatrixXd I_M=MatrixXd::Identity(power(2,m+1),power(2,m+1));
	MatrixXd s_z=tensor_product(I_m, sigma_z);
	MatrixXd s_p=tensor_product(I_m, sigma_p);
	MatrixXd s_n=tensor_product(I_m, sigma_n); // The matrices for the single spin on the edge
	
	vector<MatrixXd> correlation_after(m);
	MatrixXd s_z_after(power(2,m), power(2,m));
	MatrixXd s_p_after(power(2,m), power(2,m));
	MatrixXd s_n_after(power(2,m), power(2,m)); // The operator matrices after each truncation
	
	vector<MatrixXd> nearest_correlation(2*m+1); // In the superblock, the bond matrices
	MatrixXd H(power(2, 2*m+2), power(2, 2*m+2)); // The superblock Hamiltonian matrix
	double groundstate_energy;
	VectorXd groundstate_vector(power(2, 2*m+2));
	
	MatrixXd density(power(2, m+1), power(2, m+1));
	MatrixXd TO(power(2,m+1), power(2,m));
	MatrixXd TT(power(2,m), power(2,m+1));
	
	int size_of_system=m+1;
	infinite_write(correlation, size_of_system);
	
	cout << "The infinite system DMRG part starts......" << endl;
	// Start the infinite system DMRG part 
	
	while(size_of_system<=size_of_chain/2)
	{
		nearest_correlation.resize(2*size_of_system-1);
		for(i=0;i<size_of_system-1;i++)
		{
			nearest_correlation[i]=tensor_product(correlation[i],I_M);
			nearest_correlation[2*size_of_system-2-i]=tensor_product(I_M, correlation[i]);
		}
		nearest_correlation[size_of_system-1]=tensor_product(s_z,s_z)+0.5*tensor_product(s_p,s_n)+0.5*tensor_product(s_n,s_p);
		H=MatrixXd::Zero(power(2, 2*m+2),power(2, 2*m+2));
		for(i=0;i<2*size_of_system-1;i++)
			H=H+nearest_correlation[i];// The superblock Hamiltonian matrix

		groundstate_energy=groundstate(H, groundstate_vector); // The ground state of the Hamiltonian
		cout << size_of_system << '\t' << size_of_system << '\t' << setprecision(12) << groundstate_energy/(2*size_of_system) << endl;
		density=reduced_density_matrix(groundstate_vector, m+1);
		TO=transformation_matrix(density, m);
		TT=TO.transpose(); // The transformation matrix
        
		correlation_after.resize(size_of_system-1);
		for(i=0;i<size_of_system-1;i++)
			correlation_after[i]=TT*correlation[i]*TO;
		s_z_after=TT*s_z*TO;
		s_p_after=TT*s_p*TO;
		s_n_after=TT*s_n*TO; // The new matrix representations of the operators

		size_of_system=size_of_system+1; // The new block Hamiltonian
		correlation.resize(size_of_system-1);
		for(i=0;i<size_of_system-2;i++)
			correlation[i]=tensor_product(correlation_after[i], I2);
		correlation[size_of_system-2]=tensor_product(s_z_after,sigma_z)+0.5*tensor_product(s_p_after,sigma_n)+0.5*tensor_product(s_n_after, sigma_p);
		infinite_write(correlation, size_of_system);
	} // End of the infinite system DMRG part

	cout << "The infinite system DMRG part has been finished!" << endl;


    cout << endl << "Then the finite system DMRG part starts......" << endl;
	// Start the finite system DMRG part
	int iter;
	size_of_system=size_of_chain/2;
	finite_read(correlation, size_of_system, 1); // Initialize the finite system DMRG part

	int size_of_environment;
	vector<MatrixXd> correlation_environment(size_of_chain/2);

	for(iter=0; iter<halfsweeps; iter++)
	{
		while(size_of_chain-size_of_system>m)
		{
			size_of_environment=size_of_chain-size_of_system;
			correlation_environment.resize(size_of_environment);
			finite_read(correlation_environment, size_of_environment, iter);
			// Read the environment block

			for(i=0;i<size_of_system-1;i++)
				nearest_correlation[i]=tensor_product(correlation[i], I_M);
			for(i=0;i<size_of_environment-1;i++)
				nearest_correlation[size_of_chain-2-i]=tensor_product(I_M, correlation_environment[i]);
			nearest_correlation[size_of_system-1]=tensor_product(s_z,s_z)+0.5*tensor_product(s_p,s_n)+0.5*tensor_product(s_n,s_p);
			H=MatrixXd::Zero(power(2, 2*m+2),power(2, 2*m+2));
		    for(i=0;i<size_of_chain-1;i++)
				H=H+nearest_correlation[i];// The superblock Hamiltonian matrix

			groundstate_energy=groundstate(H, groundstate_vector); // The ground state of the Hamiltonian
		    cout << size_of_system << '\t' << size_of_environment << '\t' << setprecision(12) <<groundstate_energy/(size_of_chain) << endl;
			printf_bonds(nearest_correlation, groundstate_vector, iter, size_of_system);

		    density=reduced_density_matrix(groundstate_vector, m+1);
		    TO=transformation_matrix(density, m);
		    TT=TO.transpose(); // The transformation matrix
			
			correlation_after.resize(size_of_system-1);
		    for(i=0;i<size_of_system-1;i++)
				correlation_after[i]=TT*correlation[i]*TO;
			s_z_after=TT*s_z*TO;
			s_p_after=TT*s_p*TO;
			s_n_after=TT*s_n*TO; // The new matrix representations of the operators

			size_of_system=size_of_system+1; // The new block Hamiltonian
			correlation.resize(size_of_system-1);
			for(i=0;i<size_of_system-2;i++)
				correlation[i]=tensor_product(correlation_after[i], I2);
			correlation[size_of_system-2]=tensor_product(s_z_after,sigma_z)+0.5*tensor_product(s_p_after,sigma_n)+0.5*tensor_product(s_n_after, sigma_p);
			finite_write(correlation, size_of_system, iter);
		}
		size_of_system=m+1;
		finite_read(correlation, size_of_system, iter);
	}// End of the finite system DMRG part
}

