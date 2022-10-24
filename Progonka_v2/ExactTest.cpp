#include <iostream>
#include "Progonka.h"

void runExactTest() {
	std::cout << "ExactTest:\n";
	int in = 1;
	int n = 2;
	int out = in + n;
	int size = out + 2;
	double* e = new double[size];
	double* k = new double[size];
	double* R = new double[size];
	e[1] = 0;
	e[2] = 0.5;
	e[3] = 0.5;
	e[4] = 0;
	for (int i = in; i <= out + 1; i++) {
		k[i] = 0.1;
		R[i] = 0.1;
	}
	R[1] = 0.2;

	Coeff coeff = Coeff(n);
	coeff.el = e;
	coeff.er = e;
	coeff.kl = k;
	coeff.kr = k;
	coeff.Rl = R;
	coeff.Rr = R;
	double* ct_res_P = new double[size];		// ExactTest results: Phi_P, Phi_M
	double* ct_res_M = new double[size];
	Progonka(in, out, coeff, ct_res_P, ct_res_M);
	for (int i = in; i <= out; i++) {
		std::cout << ct_res_P[i] << "\n";
	}
}