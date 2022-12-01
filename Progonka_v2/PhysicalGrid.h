#pragma once
#include <iostream>

struct PhysicalGrid {
	const int in = 1;
	int num;
	int out;
	int size;

	double* SigmaA_l, * SigmaA_r, * SigmaS_l, * SigmaS_r;
	double* Q_l, * Q_r;

	double B_in, B_out;

	PhysicalGrid(int num);
	void print(std::ostream* f_out);
};

