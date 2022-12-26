#ifndef _wave_hpp_
#define _wave_hpp_

#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SparseLU>
#include <complex>
#include <cmath>



class wave	{

public:
	wave(int nn, int tt, double k0, double deltak, double xmin, double xmax, double deltat);
	void init_psi();
	void write_psi();
	void setMatrix();
	void updateR();
	void setV(int ver);
	void step();
	void simulationLoop();
	void printMatrix();
	
protected:	
	std::ofstream outputFile;
	
	Eigen::VectorXcd r, psi;
	Eigen::PartialPivLU<Eigen::MatrixXcd> AA;
	Eigen::MatrixXcd A;

	const double pi = 3.14159;
	double k_0, delta_k;
	double x_min, x_max, delta_x;
	double delta_t;
	double* V;

	int t;
	int n;
		
};


#endif
