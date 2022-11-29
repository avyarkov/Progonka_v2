#include <cmath>
#include "Coeff.h"
#include "OldDataGrid.h"
#include "DataGrid.h"

/*
All incoming double* arrays must have (out + 2) elements.
Input arrays: dtau, sigma_0, sigma_1, T_P, T_M;
Input boundary conditions: double B_in, double B_out; (considering A = 3/2 for now, maybe add later)
Resulting coefficients are written to: el, er, kl, kr, Rl, Rr
*/
Coeff getPlaneCoefficients(OldDataGrid dg) {
	// considering g = 0: plane case
	int in = dg.in, out = dg.out;
	dg.updateSigma();
	double* dtau = dg.dtau, * sigma_0 = dg.sigma_0, * sigma_1 = dg.sigma_1, * sigma = dg.sigma, * T_P = dg.T_P, * T_M = dg.T_M;
	double B_in = dg.B_in, B_out = dg.B_out;
	Coeff coeff = Coeff(dg.num);
	double* el = coeff.el, * er = coeff.er, * kl = coeff.kl, * kr = coeff.kr, * Rl = coeff.Rl, * Rr = coeff.Rr;

	for (int i = in + 1; i <= out; i++) {
		el[i] = sigma[i] / sigma_1[i] / sinh(dtau[i]);
		er[i] = sigma[i] / sigma_1[i] / sinh(dtau[i]);
		kl[i] = el[i] * (cosh(dtau[i]) - 1);
		kr[i] = er[i] * (cosh(dtau[i]) - 1);
		Rl[i] = el[i] * sigma_1[i] / pow(sigma[i], 2) * T_P[i] * (cosh(dtau[i]) - 1) +
			el[i] * T_M[i] / sigma[i] * sinh(dtau[i]);
		Rr[i] = er[i] * sigma_1[i] / pow(sigma[i], 2) * T_P[i] * (cosh(dtau[i]) - 1) -
			er[i] * T_M[i] / sigma[i] * sinh(dtau[i]);
	}

	// TODO: check, add chi, use the full formula
	el[in] = 0;
	kl[in] = 3.0 / 2;			// considering A = 3/2
	Rl[in] = B_in;
	er[out + 1] = 0;
	kr[out + 1] = 3.0 / 2;		// considering A = 3/2
	Rr[out + 1] = B_out;
	return coeff;
}

Coeff getPlaneCoefficients(OldDataGridNodeSources dg) {
	// considering g = 0: plain case
	int in = dg.in, out = dg.out;
	dg.updateSigma();
	double* dtau = dg.dtau, * sigma_0 = dg.sigma_0, * sigma_1 = dg.sigma_1, * sigma = dg.sigma, * T_P = dg.T_P, * T_M = dg.T_M;
	double B_in = dg.B_in, B_out = dg.B_out;
	Coeff coeff = Coeff(dg.num);
	double* el = coeff.el, * er = coeff.er, * kl = coeff.kl, * kr = coeff.kr, * Rl = coeff.Rl, * Rr = coeff.Rr;

	for (int i = in + 1; i <= out; i++) {
		el[i] = sigma[i] / sigma_1[i] / sinh(dtau[i]);
		er[i] = sigma[i] / sigma_1[i] / sinh(dtau[i]);
		kl[i] = el[i] * (cosh(dtau[i]) - 1);
		kr[i] = er[i] * (cosh(dtau[i]) - 1);
		Rl[i] = el[i] * sigma_1[i] / pow(sigma[i], 2) * T_P[i] * (cosh(dtau[i]) - 1) +
			el[i] * T_M[i] / sigma[i] * sinh(dtau[i]);
		Rr[i] = er[i] * sigma_1[i] / pow(sigma[i], 2) * T_P[i - 1] * (cosh(dtau[i]) - 1) -
			er[i] * T_M[i - 1] / sigma[i] * sinh(dtau[i]);
	}

	// TODO: check, add chi, use the full formula
	el[in] = 0;
	kl[in] = 3.0 / 2;			// considering A = 3/2
	Rl[in] = B_in;
	er[out + 1] = 0;
	kr[out + 1] = 3.0 / 2;		// considering A = 3/2
	Rr[out + 1] = B_out;
	return coeff;
}

Coeff getPlaneCoefficients(DataGrid dg) {
	// considering g = 0: plain case
	dg.updateSigma();
	dg.updateDtau();
	int in = dg.in, out = dg.out;
	double* dtau = dg.dtau;
	double* sigma0_l = dg.sigma0_l, * sigma0_r = dg.sigma0_r, * sigma1_l = dg.sigma1_l, * sigma1_r = dg.sigma1_r;
	double* sigma_l = dg.sigma_l, * sigma_r = dg.sigma_r;
	double* TP_l = dg.TP_l, * TP_r = dg.TP_r, * TM_l = dg.TM_l, * TM_r = dg.TM_r;
	double B_in = dg.B_in, B_out = dg.B_out;
	Coeff coeff = Coeff(dg.num);
	double* el = coeff.el, * er = coeff.er, * kl = coeff.kl, * kr = coeff.kr, * Rl = coeff.Rl, * Rr = coeff.Rr;

	for (int i = in + 1; i <= out; i++) {
		el[i] = sigma_l[i] / sigma1_l[i] / sinh(dtau[i]);
		er[i] = sigma_r[i - 1] / sigma1_r[i - 1] / sinh(dtau[i]);
		kl[i] = el[i] * (cosh(dtau[i]) - 1);
		kr[i] = er[i] * (cosh(dtau[i]) - 1);
		Rl[i] = el[i] * sigma1_l[i] / pow(sigma_l[i], 2) * TP_l[i] * (cosh(dtau[i]) - 1) +
			el[i] * TM_l[i] / sigma_l[i] * sinh(dtau[i]);
		Rr[i] = er[i] * sigma1_r[i - 1] / pow(sigma_r[i - 1], 2) * TP_r[i - 1] * (cosh(dtau[i]) - 1) -
			er[i] * TM_r[i - 1] / sigma_r[i - 1] * sinh(dtau[i]);
	}

	// TODO: check, add chi, use the full formula
	el[in] = 0;
	kl[in] = 3.0 / 2;			// considering A = 3/2
	Rl[in] = B_in;
	er[out + 1] = 0;
	kr[out + 1] = 3.0 / 2;		// considering A = 3/2
	Rr[out + 1] = B_out;
	return coeff;
}