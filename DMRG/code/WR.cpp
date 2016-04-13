/* Jin Xu, a graduate student */
/* Department of Physics, University of Alberta, Edmonton, Alberta, T6G 2E1, Canada */
/* PHYS 580, the term project */

/* This library is for matrix write and matrix read */

#include <Eigen/Dense>
using namespace Eigen;

#include <iostream>
using namespace std;

#include <fstream>
using namespace std;

#include <vector>
using namespace std;

#include <cstdlib>
using namespace std;

#include <iomanip>
using namespace std;


void save(const char *filename, const MatrixXd& m)
{
	// Write the matrix m to a binary file
	ofstream f(filename, ios::binary);
	if(!f.is_open())
	{
		cerr << "ERROR: couldn't open " << filename << "!" << endl;
		exit(1);
	}
	int ROWS=m.rows();
	int COLS=m.cols();
	f.write((char *)&ROWS, sizeof(int));
	f.write((char *)&COLS, sizeof(int));
	f.write((char *)m.data(), sizeof(double)*ROWS*COLS);
	f.close();
}

void load(const char *filename, MatrixXd& m)
{
	// Read the matrix m from a binary file
	int rows, cols;
	ifstream f(filename, ios::binary);
	if(!f.is_open())
	{
		cerr << "ERROR: couldn't open " << filename << "!" << endl;
		exit(1);
	}
	f.read((char *)&rows, sizeof(int));
    f.read((char *)&cols, sizeof(int));
    m.resize(rows, cols);
	f.read((char *)m.data(), sizeof(double)*rows*cols);
	if (f.bad())
	{
		cerr << "Matrix Reading Error!" << endl;
		exit(1);
	}
	f.close();
} 

void infinite_write(vector<MatrixXd> correlation, int system_size)
{
	// In the infinite system DMRG part
	// Write the matrices representing H(i,i+1) to binary files
	// The vector<MatriXd> size should be system_size-1
	// The maximum size of the chain should be smaller than 1000
	
	char filename[12];
	filename[7]='.';
	filename[8]='d';
	filename[9]='a';
	filename[10]='t';
	filename[11]='\0';
	filename[1]=48+system_size/100;
	filename[2]=48+(system_size%100)/10;
	filename[3]=48+system_size%10; 
	
	int i;
	for(i=0;i<=system_size-2;i++)
	{
		filename[4]=48+i/100;
		filename[5]=48+(i%100)/10;
		filename[6]=48+i%10; 
		
		filename[0]='L';
		save(filename, correlation[i]);
		filename[0]='R';
		save(filename, correlation[i]);
	}
}

void finite_write(vector<MatrixXd> correlation, int system_size, int iter)
{
	// In the finite system DMRG part
	// Write the matrices representing H(i,i+1) to binary files
	// The vector<MatrixXd> size should be system_size-1
	// The maximum size of the chain should be smaller than 1000
	
	char filename[12];
	filename[7]='.';
	filename[8]='d';
	filename[9]='a';
	filename[10]='t';
	filename[11]='\0';
	filename[1]=48+system_size/100;
	filename[2]=48+(system_size%100)/10;
	filename[3]=48+system_size%10; 
	if(iter%2==0)
		filename[0]='L';
	else
		filename[0]='R';
	
	int i;
	for(i=0;i<=system_size-2;i++)
	{
		filename[4]=48+i/100;
		filename[5]=48+(i%100)/10;
		filename[6]=48+i%10; 
		save(filename, correlation[i]);
	}
}

void finite_read(vector<MatrixXd>& correlation, int system_size, int iter)
{
	// In the finite system DMRG part
	// Read the matrices representing H(i,i+1) from binary files
	// The vector<MatrixXd> size should be system_size-1
	// The maximum size of the chain should be smaller than 1000
	correlation.resize(system_size-1);
	char filename[12];
	filename[7]='.';
	filename[8]='d';
	filename[9]='a';
	filename[10]='t';
	filename[11]='\0';
	filename[1]=48+system_size/100;
	filename[2]=48+(system_size%100)/10;
	filename[3]=48+system_size%10; 
	if(iter%2==0)
		filename[0]='R';
	else
		filename[0]='L';
	
	int i;
	for(i=0;i<=system_size-2;i++)
	{
		filename[4]=48+i/100;
		filename[5]=48+(i%100)/10;
		filename[6]=48+i%10; 
		load(filename, correlation[i]);
	}
}


void printf_bonds(vector<MatrixXd> nearest_correlation, VectorXd evec, int iter, int system_size)
{
	// Calculate the bonds
	int L=nearest_correlation.size();
	
	char filename[12];
	filename[7]='.';
	filename[8]='t';
	filename[9]='x';
	filename[10]='t';
	filename[11]='\0';
	filename[0]='C';
	filename[1]=48+iter/100;
	filename[2]=48+(iter%100)/10;
	filename[3]=48+iter%10;
	filename[4]=48+system_size/100;
	filename[5]=48+(system_size%100)/10;
	filename[6]=48+system_size%10;
	
	VectorXd bond(L);
	for(int i=0; i<L; i++)
		bond(i)=evec.dot(nearest_correlation[i]*evec);
	
	ofstream f(filename);
	if(!f.is_open())
	{
		cerr << "couldn't open " << filename << "!" << endl;
		exit(1);
	}
	
	f << setprecision(12) << bond << endl;
	f.close();
}