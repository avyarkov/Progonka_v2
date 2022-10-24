#include "Coeff.h"
/*
All incoming double* arrays must have (out + 2) elements.
Input arrarys: el, er, kl, kr, Rl, Rr
Results are written to: PHI_P, PHI_M
*/
void Progonka(int in, int out, Coeff coeff, double* PHI_P, double* PHI_M) {
	int size = out + 2;
	double* el = coeff.el, * er = coeff.er, * kl = coeff.kl, * kr = coeff.kr, * Rl = coeff.Rl, * Rr = coeff.Rr;
	
	double* alpha = new double[size];
	double* beta = new double[size];
	double* epsilon = new double[size];
	double* delta = new double[size];

	alpha[out + 1] = 0;
	beta[out + 1] = 0;
	double denom = 0;
	for (int i = out; i >= in; i--) {
		denom = er[i + 1] * (1 - alpha[i + 1]) + kr[i + 1] + el[i] + kl[i];
		alpha[i] = el[i] / denom;
		beta[i] = (er[i + 1] * beta[i + 1] + Rl[i] + Rr[i + 1]) / denom;
		// epsilon[i] = -alpha[i] * (el[i] + kl[i]) + el[i];	// "old" formula, no (INF - INF) resolution
		// delta[i] = beta[i] * (el[i] + kl[i]) - Rl[i];		// "old" formula, no (INF - INF) resolution
		epsilon[i] = alpha[i] * (er[i + 1] * (1 - alpha[i + 1]) + kr[i + 1]);
		delta[i] = ((el[i] + kl[i]) * (er[i + 1] * beta[i + 1] + Rr[i + 1]) - Rl[i] * (er[i + 1] * (1 - alpha[i + 1]) + kr[i + 1])) / denom;
	}
	PHI_P[in - 1] = 0;
	for (int i = in; i <= out; i++) {
		PHI_P[i] = alpha[i] * PHI_P[i - 1] + beta[i];
		PHI_M[i] = epsilon[i] * PHI_P[i - 1] - delta[i];
	}
}