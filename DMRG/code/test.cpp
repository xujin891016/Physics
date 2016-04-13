/* Jin Xu, a graduate student */
/* Department of Physics, University of Alberta, Edmonton, Alberta, T6G 2E1, Canada */
/* PHYS 580, the term project */

/* test.cpp is to test the functions I wrote */

#include <iostream>
using namespace std;

#include <Eigen/Dense>
using namespace Eigen;

#include <vector>
using namespace std;

#include "tensor.hpp"
#include "ED.hpp"
#include "RDM.hpp"
#include "WR.hpp"

void test_tensor()
{
	// It testes the functions in the library "tensor.hpp"
	cout << "power(2, 5) = " << power(2, 5) << endl;
	
	MatrixXd sigma_z(2, 2);
	sigma_z << 0.5, 0, 0, -0.5;
	MatrixXd sigma_p(2, 2);
	sigma_p << 0, 1, 0, 0;
	MatrixXd sigma_n(2, 2);
	sigma_n << 0, 0, 1, 0;
	
	MatrixXd H;
	H=tensor_product(sigma_z, sigma_z) + 0.5*tensor_product(sigma_p, sigma_n) + 0.5*tensor_product(sigma_n, sigma_p);
	
	cout << "H = " << endl;
	cout << H <<endl;
}

void test_ED()
{
	// It tests the functions in the library "ED.hpp"
	VectorXd V(6);
	V << 1.0, -1.0, -1.0, 1.0, 0.5, -0.5;
	cout << "The vector V = " << endl;
	cout << V <<endl;

	cout << "The index of the minimum element in the vector V is " << endl;
	cout << index_min(V) << endl;

	cout << "After the vector V is sorted, the order should be as follows: " << endl;
	cout << sort_vector(V) << endl;  
	
	cout << "After the sorting procedure, V = " << endl;
	cout << V << endl;

	MatrixXd sigma_z(2, 2);
	sigma_z << 0.5, 0, 0, -0.5;
	MatrixXd sigma_p(2, 2);
	sigma_p << 0, 1, 0, 0;
	MatrixXd sigma_n(2, 2);
	sigma_n << 0, 0, 1, 0;
	
	MatrixXd H;
	H=tensor_product(sigma_z, sigma_z) + 0.5*tensor_product(sigma_p, sigma_n) + 0.5*tensor_product(sigma_n, sigma_p);
	cout << "Before the exact diagonalization, H = " << endl;
	cout << H << endl;
	
	SelfAdjointEigenSolver <MatrixXd> eigensolver(H);
	cout << "The eigenvalues are " << endl;
	cout << eigensolver.eigenvalues() << endl;
	cout << "The eigenvectors are " << endl;
	cout << eigensolver.eigenvectors() << endl;

	VectorXd G_state(4);
	double groundstate_energy;
	groundstate_energy=groundstate(H, G_state);
	cout << "The groundstate energy is " << groundstate_energy << endl;
	cout << "The groundstate vector is " << endl;
	cout << G_state <<endl;

	VectorXd eval(4);
	MatrixXd evec(4,4);
	diagonalization(H, eval, evec);
	cout << "The eigenvalues are " << endl;
	cout << eval << endl;
	cout << "The eigenvectors are " << endl;
	cout << evec << endl;
	
	cout << "After the exact diagonalization, H = " << endl;
	cout << H << endl;
}

void test_RDM()
{
	// It tests the functions in the library "RDM.hpp"
	MatrixXd sigma_z(2, 2);
	sigma_z << 0.5, 0, 0, -0.5;
	MatrixXd sigma_p(2, 2);
	sigma_p << 0, 1, 0, 0;
	MatrixXd sigma_n(2, 2);
	sigma_n << 0, 0, 1, 0;
	
	MatrixXd H;
	H=tensor_product(sigma_z, sigma_z) + 0.5*tensor_product(sigma_p, sigma_n) + 0.5*tensor_product(sigma_n, sigma_p);
	cout << "The two-size block Hamiltonian matrix is " << endl;
	cout << H << endl;
	
	MatrixXd I2=MatrixXd::Identity(2,2);
	MatrixXd I4=MatrixXd::Identity(4,4);
	MatrixXd s_z=tensor_product(I2,sigma_z);
	MatrixXd s_p=tensor_product(I2,sigma_p);
	MatrixXd s_n=tensor_product(I2,sigma_n);
	
	MatrixXd superblock=tensor_product(H,I4)+tensor_product(I4,H)+tensor_product(s_z,s_z)+0.5*tensor_product(s_p,s_n)+0.5*tensor_product(s_n,s_p);
	cout << "The four-size superblock Hamiltonian matrix is " << endl;
	cout << superblock << endl;
	
	double groundstate_energy;
	VectorXd groundstate_vector(16);
	groundstate_energy=groundstate(superblock, groundstate_vector);
	cout << "The groundstate energy for the superblock Hamiltonian is " << endl;
	cout << groundstate_energy << endl;
	cout << "The vector representing the ground state of the superblock Hamiltonian is " << endl;
	cout << groundstate_vector << endl;
	
	MatrixXd density=reduced_density_matrix(groundstate_vector, 2);
	cout << "The reduced density matrix is " << endl;
	cout << density << endl;
	
	SelfAdjointEigenSolver<MatrixXd> eigensolver(density);
	cout << "The eigenvalues of the matrix density are " << endl;
	cout << eigensolver.eigenvalues() << endl;
	cout << "The eigenvectors of the matrix density are " << endl;
	cout << eigensolver.eigenvectors() << endl;
	
	MatrixXd T=transformation_matrix(density,1);
	cout << "The transformation matrix T = " << endl;
	cout << T << endl;
}

void test_WR()
{
	// It tests the functions in the library "WR.hpp"
	MatrixXd sigma_z(2, 2);
	sigma_z << 0.5, 0, 0, -0.5;
	MatrixXd sigma_p(2, 2);
	sigma_p << 0, 1, 0, 0;
	MatrixXd sigma_n(2, 2);
	sigma_n << 0, 0, 1, 0;
	MatrixXd I2=MatrixXd::Identity(2,2);
	
	MatrixXd H;
	H=tensor_product(sigma_z, sigma_z) + 0.5*tensor_product(sigma_p, sigma_n) + 0.5*tensor_product(sigma_n, sigma_p);
	MatrixXd I4=MatrixXd::Identity(4,4);
	MatrixXd s_z=tensor_product(I2,sigma_z);
	MatrixXd s_p=tensor_product(I2,sigma_p);
	MatrixXd s_n=tensor_product(I2,sigma_n);
	
	MatrixXd superblock=tensor_product(H,I4)+tensor_product(I4,H)+tensor_product(s_z,s_z)+0.5*tensor_product(s_p,s_n)+0.5*tensor_product(s_n,s_p);
	
	cout << "The superblock Hamiltonian matrix dumped to the file test.dat is " << endl;
	cout << superblock << endl;
	char filename1[9];
	filename1[0]='t';
	filename1[1]='e';
	filename1[2]='s';
	filename1[3]='t';
	filename1[4]='.';
	filename1[5]='d';
	filename1[6]='a';
	filename1[7]='t';
	filename1[8]='\0';
	save(filename1, superblock);
	
	MatrixXd superblock2;
	load(filename1, superblock2);
	cout << "The superblock Hamiltonian matrix read from the file test.dat is " << endl;
	cout << superblock2 << endl;
	
	vector<MatrixXd> correlation1(1);
	correlation1[0]=H;
	cout << "The block Hamiltonian matrix written to the disk is " << endl;
	cout << correlation1[0] << endl;
	infinite_write(correlation1, 2);
	finite_write(correlation1, 2, 0);
	
	vector<MatrixXd> correlation2(1);
	finite_read(correlation2, 2, 0);
	cout << "The block Hamiltonian matrix read from the disk is " << endl;
	cout << correlation2[0] << endl;
	
	vector<MatrixXd> nearest_correlation(3);
	nearest_correlation[0]=tensor_product(H, I4);
	nearest_correlation[1]=tensor_product(s_z,s_z)+0.5*tensor_product(s_p,s_n)+0.5*tensor_product(s_n,s_p);
	nearest_correlation[2]=tensor_product(I4, H);
	
	VectorXd groundstate_vector;
	double groundstate_energy=groundstate(superblock, groundstate_vector);
	printf_bonds(nearest_correlation, groundstate_vector, 0, 2);
}
	

int main(void)
{
	test_tensor();
	test_ED();
	test_RDM();
	test_WR();
	return 0;
}
	
	
	

