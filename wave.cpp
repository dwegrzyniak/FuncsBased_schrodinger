


#include "wave.hpp"

using namespace std;

wave::wave(int nn, int tt, double k0, double deltak, double xmin, double xmax, double deltat){
	
	t = tt;
	n = nn;
	k_0 = k0;
	delta_k = deltak;
	x_min = xmin;
	x_max = xmax;
	delta_t = deltat;	
	delta_x = (x_max - x_min) / (n - 1);
	
	A = Eigen::MatrixXcd::Zero(n, n);
	r = Eigen::VectorXcd::Zero(n, 1);
	psi = Eigen::VectorXcd::Zero(n, 1); 
	
	outputFile.open("output_z_15.txt");
	
}

void wave::simulationLoop(){
	
	for (double i = 0; i<t; i+=delta_t){
		cout << i << endl;
		updateR();
		step();
		if (int(i / delta_t) % 3 == 0) write_psi();
	}
	
}

void wave::write_psi(){
	for (int it = 0; it < n; it++) {
		outputFile << x_min + it * delta_x << ",";	
	}
	outputFile << '\n' << endl;
	
	for (int it = 0; it < n; it++) {
		outputFile << norm(psi(it)) << "," ;
	}
	
	outputFile << '\n' << endl;
	
	for (int it = 0; it < n; it++) {
		outputFile << psi(it).real() << "," ;
	}
	
	outputFile << '\n' << endl;
	
	for (int it = 0; it < n; it++) {
		outputFile << psi(it).imag() <<",";
	}
	outputFile << '\n' << endl;
	
	for (int it = 0; it < n; it++) {
		outputFile << V[it] <<",";
	}
	outputFile << '\n' << endl;
}

void wave::init_psi(){
	
	for (int it = 0; it < n; it++) {
		
		if ( it == 0 || it == n - 1){
			psi(it) = 0;
		}
		
		else {
			double x = it * delta_x + x_min;
			psi(it) = exp(1i * k_0 * x) * exp(- 0.5 * pow(x * delta_k, 2 ) ) * pow(delta_k, 0.5) / pow(pi, 0.25);
			
		}
	}
	
}

void wave::printMatrix(){
	
	for (int i = 0; i < n ; i++){
		for (int j = 0; j < n; j++){
				cout << A(i, j);
			
		}
		cout << '\n' << endl;
	}
	/*
	for (int i = 0; i < n ; i++){
		cout << r(i)<< endl;
	}
	*/
}

void wave::setMatrix(){
	std::complex <double> ii(0, 1);
	for (int i = 0; i < n; i++){
		A(i, i) = 1. + 0.5 * delta_t * ii * (2 / (delta_x * delta_x) + V[i]);
		if (i != 0) A(i, i - 1) = -delta_t * ii * (1 / (2 *delta_x * delta_x)); 
		if (i !=  n-1) A(i, i + 1) = -delta_t * ii * (1 / (2 *delta_x * delta_x)); 
	}
	
	AA=A.lu();
}

void wave::updateR(){
	complex <double> ii(0, 1);
	r[0] = psi[0] + ii * delta_t * 0.5 * ((psi[1] - 2. * psi[0]) / (delta_x * delta_x) - V[0] * psi[0]);
	r[n - 1] = psi[0] + ii * delta_t * 0.5 * ((- 2. * psi[n - 1]) + psi[n - 2]/ (delta_x * delta_x) - V[n-1] * psi[n-1]);
	for (int i = 1; i < n - 1; i++){
		
		r[i] = psi[i] + ii * delta_t * 0.5 * ((psi[i - 1] - 2. * psi[i] + psi[i + 1])/(delta_x * delta_x) - V[i] * psi[i]);
		
	}
	
}

void::wave::step(){
	psi = AA.solve(r); 
	psi[0] = 0;
	psi[n - 1] = 0;
}

void::wave::setV(int ver){
	this -> V = new double[n - 2];
	if (ver == 0){
		for (int i = 0; i < n - 2; i++){
			V[i] = 0;
		}
		
	}
	if (ver == 1){
		double V0 = 105., sig = 0.5;
		for (int i = 0; i < n - 2; i++){
			V[i] = V0 * exp(-(i * delta_x + x_min - 10 )*(i * delta_x + x_min -10)/(sig*sig));
		}
		
	}
	
}
