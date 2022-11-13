#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "DataGrid.h"
#include "FileIO.h"
#include "Coeff.h"
#include "PlaneCase.h"
#include "Progonka.h"
#include "Quadrature.h"

using namespace std;

double* countGamma(double* D, int in, int out) {
	double* gamma = new double[DataGrid::sizeFromOut(out)];
	for (int i = in; i <= out; i++) {
		gamma[i] = 3.0 - 1 / D[i];
	}
	return gamma;
}

double calcF_P(double mu, double D, double PSI_0, double w_0, double w_2) {
	return (w_0 + 5.0 / 4 * w_2 * (3 * mu * mu - 1) * (3 * D - 1)) * PSI_0 / 4 / M_PI;
}

double calcF_M(double mu, double PSI_1, double w_1) {
	return 3 * w_1 / 4 / M_PI * mu * PSI_1;
}

void runQuasidiffusion() {
	ofstream fout("Quasidiffusion_output.txt");

	const double w_0 = 1.0, w_1 = 0.5, w_2 = 0.2;
	const double nu = 1;
	DataGrid quasiDG = readDataGrid("QuasidiffusionDataGrid.txt");
	quasiDG.updateSigma();
	quasiDG.updateDtau();
	//quasiDG.print(&fout);
	DataGrid eoDG = readDataGrid("EvenOddDataGrid.txt");
	eoDG.updateSigma();
	eoDG.updateDtau();

	int num = quasiDG.num;
	if (num != eoDG.num) { cout << "Incorrect input: quasiDG and eoDg are of different sizes!"; exit(239); }
	int size = quasiDG.size, in = quasiDG.in, out = quasiDG.out;
	double* Sigma_s = new double[size];
	for (int i = in + 1; i <= out; i++) { Sigma_s[i] = quasiDG.sigma_1[i] - quasiDG.sigma_0[i]; }
	double* D = new double[size];
	for (int i = in; i <= out; i++) { D[i] = 1.0 / 3; }
	double* gamma = countGamma(D, in, out);
	double A_in = 3.0 / 2, A_out = 3.0 / 2;
	int numberOfIterations = 5;
	Quadrature quadrature = Quadrature("Quadrature_Gauss–Legendre.txt");
	int order = 6;
	double* root = quadrature.rootsOfOrder(order);
	double* weight = quadrature.weightsOfOrder(order);
	for (int it = 0; it < numberOfIterations; it++) {
		// solve QuasidiffusionEquations
		DataGrid curQuasiDg = quasiDG.withDividedSigma0(D);
		curQuasiDg.updateSigma();
		curQuasiDg.updateDtau();
		Coeff quasiCoeff = getPlaneCoefficients(curQuasiDg);
		quasiCoeff.updateBoundaries(A_in, A_out);
		double* PSI_0 = new double[size], * PSI_1 = new double[size];
		Progonka(in, out, quasiCoeff, PSI_0, PSI_1);
		writeArrayEndl(PSI_0, in, out, &fout);

		double** F_P = new double* [order], ** F_M = new double* [order];
		double** phi_p = new double* [order], ** phi_m = new double* [order], ** phi = new double* [order];
		for (int k = 0; k < order; k++) {
			F_P[k] = new double[size];
			F_M[k] = new double[size];
			phi_p[k] = new double[size];
			phi_m[k] = new double[size];
			phi[k] = new double[size];
			for (int i = in; i <= out; i++) {
				F_P[k][i] = calcF_P(root[k], D[i], PSI_0[i], w_0, w_2);
				F_M[k][i] = calcF_M(root[k], PSI_1[i], w_1);
			}
			double* Add_P = new double[size], * Add_M = new double[size];
			for (int i = in + 1; i <= out; i++) {
				Add_P[i] = nu * Sigma_s[i] * F_P[k][i];
				Add_M[i] = nu * Sigma_s[i] * F_M[k][i];
			}
			// TODO check negative roots
			DataGrid curEODG = eoDG.toCharacteristics(root[k]).withAddedSources(Add_P, Add_M);
			curEODG.updateSigma();
			curEODG.updateDtau();
			Coeff curCoeff = getPlaneCoefficients(curEODG);
			curCoeff.updateBoundaries(1.0, 1.0);
			Progonka(in, out, curCoeff, phi_p[k], phi_m[k]);
			writeArrayEndl(phi_p[k], in, out, &fout);
			writeArrayEndl(phi_m[k], in, out, &fout);
			for (int i = in; i <= out; i++) { phi[k][i] = phi_p[k][i] + phi_m[k][i]; }
		}
		for (int i = in; i <= out; i++) {
			double* f2 = new double[order];
			double* curphi = new double[order];
			for (int k = 0; k < order; k++) {
				curphi[k] = phi_p[k][i];
				f2[k] = root[k] * root[k] * curphi[k];
			}
			double integral0 = quadrature.integrate(order, curphi);
			double integral2 = quadrature.integrate(order, f2);
			D[i] = integral2 / integral0;
			if (i == in) {
				double* f1 = new double[order];
				for (int k = 0; k < order; k++) { f1[k] = abs(root[k]) * curphi[k]; }
				double integral1 = quadrature.integrate(order, f1);
				A_in = integral1 / integral2;
			}
			if (i == out) {
				double* f1 = new double[order];
				for (int k = 0; k < order; k++) { f1[k] = abs(root[k]) * curphi[k]; }
				double integral1 = quadrature.integrate(order, f1);
				A_out = integral1 / integral2;
			}
		}
		writeArrayEndl(D, in, out, &fout);
		fout << A_in << " " << A_out << endl << "~~~~~~~~~~~" << endl;
	}
}