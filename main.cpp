#include "wave.hpp"

/*
 * dk = 0.5 x_min = -5 x_max = 50, n = 1000, t = 4, dt = 0.002
 * dk = 0.75 x_min = -5 x_max = 25, n = 1000, t = 4, dt = 0.002
 * dk = 1. x_min = -5 x_max = 25, n = 1000, t = 4, dt = 0.002
 * dk = 1.5 x_min = -5 x_max = 25, n = 1000, t = 4, dt = 0.002
 * 
 */ 


int main() {
	wave wav(1000, 4, 10, 1.5, -5, 25, 0.002); //n(J), t, k0, deltak, xmin, xmax, deltat
	wav.setV(1);
	wav.setMatrix();
	wav.init_psi();
	//wav.printMatrix();
	wav.simulationLoop();


return 0;
}
