#include <cmath>
#include "DataGrid.h"
#include "Coeff.h"

/*
All incoming double* arrays must have (out + 2) elements.
Input arrays: dtau, r, sigma_0, sigma_1, T_P, T_M, gamma;
Input boundary conditions: double B_in, double B_out; (considering A = 3/2 for now, maybe add later)
Resulting coefficients are written to: el, er, kl, kr, Rl, Rr
Spherical case, g = 1
*/
Coeff getSphericalCoefficients(DataGrid dg, double* gamma) {
	int g = 1; // considering g = 1: spherical case
	int in = dg.in, out = dg.out;
	dg.updateSigma();
	double* r = dg.r, * dtau = dg.dtau, * sigma_0 = dg.sigma_0, * sigma_1 = dg.sigma_1, * sigma = dg.sigma, * T_P = dg.T_P, * T_M = dg.T_M;
	double B_in = dg.B_in, B_out = dg.B_out;
	Coeff coeff = Coeff(dg.num);
	double* el = coeff.el, * er = coeff.er, * kl = coeff.kl, * kr = coeff.kr, * Rl = coeff.Rl, * Rr = coeff.Rr;

	for (int i = in + 1; i <= out; i++) {
		el[i] = sigma[i] / sigma_1[i] / sinh(dtau[i]);
		er[i] = sigma[i] / sigma_1[i] / sinh(dtau[i]);
		kl[i] = el[i] * (cosh(dtau[i]) - 1) - g * (1 - gamma[i]) / r[i] / sigma_1[i];
		kr[i] = er[i] * (cosh(dtau[i]) - 1) + g * (1 - gamma[i - 1]) / r[i - 1] / sigma_1[i];
		Rl[i] = el[i] * sigma_1[i] / pow(sigma[i], 2) * T_P[i] * (cosh(dtau[i]) - pow(r[i - 1] / r[i], g) - g * sinh(dtau[i]) / pow(r[i], g) / sigma[i]) +
			el[i] * T_M[i] / sigma[i] * (sinh(dtau[i]) - g * (2 - gamma[i]) * (cosh(dtau[i]) - 1) / pow(r[i], g) / sigma[i]);
		Rr[i] = er[i] * sigma_1[i] / pow(sigma[i], 2) * T_P[i] * (cosh(dtau[i]) - pow(r[i] / r[i - 1], g) + g * sinh(dtau[i]) / pow(r[i - 1], g) / sigma[i]) -
			er[i] * T_M[i] / sigma[i] * (sinh(dtau[i]) + g * (2 - gamma[i - 1]) * (cosh(dtau[i]) - 1) / pow(r[i - 1], g) / sigma[i]);
	}
	// TODO: check, add chi, use the full formula
	if (r[in] == 0) {
		el[in] = 0;
		kl[in] = 0;			// TODO symmetry condition PHI_M[in] = 0 in the center r = r_in = 0;
		Rl[in] = 0;
	}
	else {
		el[in] = 0;
		kl[in] = 3.0 / 2;			// considering A = 3/2
		Rl[in] = B_in;
	}
	er[out + 1] = 0;
	kr[out + 1] = 3.0 / 2;		// considering A = 3/2
	Rr[out + 1] = B_out;

	return coeff;
}

void leftSphericalSupplement(DataGrid dg, double* gamma, Coeff coeff, double* PHI_P, double* PHI_M) {
	int g = 1; // considering g = 1: spherical case
	int in = dg.in;
	dg.updateSigma();
	double* r = dg.r, * dtau = dg.dtau, * sigma_0 = dg.sigma_0, * sigma_1 = dg.sigma_1, * sigma = dg.sigma, * T_P = dg.T_P, * T_M = dg.T_M;
	double B_in = dg.B_in, B_out = dg.B_out;
	double* el = coeff.el, * er = coeff.er, * kl = coeff.kl, * kr = coeff.kr, * Rl = coeff.Rl, * Rr = coeff.Rr;
	
	double coeff_A = g * (1 - gamma[in]) / sigma_1[in + 1];
	double coeff_B = er[in + 1] * pow(r[in + 1], g);
	double sigma_in_plus1 = sqrt(abs(sigma_0[in + 1]) * sigma_1[in + 1]);
	double coeff_C = er[in + 1] * sigma_1[in + 1] / pow(sigma_in_plus1, 2) * T_P[in + 1] * (-pow(r[in + 1], g) + g * sinh(dtau[in + 1]) / sigma_in_plus1) -
		er[in + 1] * T_M[in + 1] / sigma_in_plus1 * g * (2 - gamma[in]) * (cosh(dtau[in + 1]) - 1) / sigma_in_plus1;
	PHI_P[in] = (coeff_B * PHI_P[in + 1] + coeff_C) / coeff_A;
	PHI_M[in] = 0;
}