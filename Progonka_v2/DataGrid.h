#include <iostream>

#pragma once

struct DataGrid {
	const int in = 1;	// start index in the array; is arbitrary and can be 1 to INF, because in Progonka() index (in - 1) is used for convenience
	int num;				// number of (inner) spatial steps; количество (внутренних) €чеек
	int out;			// out := in + n; range [in, out + 1] is the range of used array elements; indexes [in + 1, out] correspond to inner steps of the grid
	int size;			// size := out + 2; required size of arrays

	double* r, * dtau;
	double* sigma_0, * sigma_1, * T_P, * T_M;
	double* sigma;

	double B_in, B_out;

	DataGrid(int num);
	static int sizeFromOut(int out);
	void print(std::ostream* f_out);

	void updateSigma();
	void updateDtau();

	DataGrid multiplied(int k);	// k > 0
	
	DataGrid withAddedSources(double* Add_P, double* Add_M);
	DataGrid toCharacteristics(double mu);
	DataGrid withDividedSigma0(double* D);
};

