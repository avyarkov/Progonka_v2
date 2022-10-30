#pragma once
#include <iostream>

struct Coeff {
	const int in = 1;	// start index in the array; is arbitrary and can be 1 to INF, because in Progonka() index (in - 1) is used for convenience
	int num;				// number of (inner) spatial steps; количество (внутренних) €чеек
	int out;			// out := in + n; range [in, out + 1] is the range of used array elements; indexes [in + 1, out] correspond to inner steps of the grid
	int size;
	
	double* el;
	double* er;
	double* kl;
	double* kr;
	double* Rl;
	double* Rr;

	Coeff(int num);
	void print(std::ostream* f_out);

	void updateBoundaries(double A_in, double A_out);
};

