#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "OldDataGrid.h"
#include "FileIO.h"
#include "Coeff.h"
#include "PlaneCase.h"
#include "Progonka.h"
#include "Quadrature.h"

using namespace std;

double* countGamma(double* D, int in, int out) {
	double* gamma = new double[OldDataGrid::sizeFromOut(out)];
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
	DataGrid quasiDG = readDataGridCellSectionsNodeSources("QuasidiffusionDataGrid.txt");
	quasiDG.updateSigma();
	quasiDG.updateDtau();
	DataGrid eoDG = readDataGridCellSectionsNodeSources("EvenOddDataGrid.txt");
	eoDG.updateSigma();
	eoDG.updateDtau();

	int num = quasiDG.num;
	if (num != eoDG.num) { cout << "Incorrect input: quasiDG and eoDg are of different sizes!"; exit(239); }
	int size = quasiDG.size, in = quasiDG.in, out = quasiDG.out;
	double Sigma_s = quasiDG.sigma1_r[in] - quasiDG.sigma0_r[in];
	
	double* D = new double[size];
	for (int i = in; i <= out; i++) { D[i] = 1.0 / 3; }
	//double* gamma = countGamma(D, in, out);
	double A_in = 3.0 / 2, A_out = 3.0 / 2;
	
	Quadrature quadrature = Quadrature("Quadrature_Gauss–Legendre.txt");
	int order = 6;
	int qin = order / 2, qout = order - 1;
	double* root = quadrature.rootsOfOrder(order);
	double* weight = quadrature.weightsOfOrder(order);
	
	int numberOfIterations = 5;
	for (int it = 0; it < numberOfIterations; it++) {
		// solve Quasidiffusion equations
		DataGrid curQuasiDg = quasiDG.withDividedSigma0(D);
		Coeff quasiCoeff = getPlaneCoefficients(curQuasiDg);
		quasiCoeff.updateBoundaries(A_in, A_out);
		quasiCoeff.print(&cout);

		double* PSI_0 = new double[size], * PSI_1 = new double[size];
		Progonka(in, out, quasiCoeff, PSI_0, PSI_1);
		writeArrayEndl(PSI_0, in, out, &fout);

		// solve Even-Odd equations for positive quadrature roots
		double** F_P = new double* [order], ** F_M = new double* [order];
		double** phi_p = new double* [order], ** phi_m = new double* [order];
		for (int k = qin; k <= qout; k++) {
			F_P[k] = new double[size];
			F_M[k] = new double[size];
			phi_p[k] = new double[size];
			phi_m[k] = new double[size];
			for (int i = in; i <= out; i++) {
				F_P[k][i] = calcF_P(root[k], D[i], PSI_0[i], w_0, w_2);
				F_M[k][i] = calcF_M(root[k], PSI_1[i], w_1);
			}
			double* Add_P = new double[size], * Add_M = new double[size];
			for (int i = in; i <= out; i++) {
				Add_P[i] = nu * Sigma_s * F_P[k][i];
				Add_M[i] = nu * Sigma_s * F_M[k][i];
			}
			DataGrid curEODG = eoDG.toPlaneCharacteristics(root[k]).withAddedSources(Add_P, Add_M);
			Coeff curCoeff = getPlaneCoefficients(curEODG);
			curCoeff.updateBoundaries(1.0, 1.0);
			Progonka(in, out, curCoeff, phi_p[k], phi_m[k]);
			writeArrayEndl(phi_p[k], in, out, &fout);
			writeArrayEndl(phi_m[k], in, out, &fout);
		}
		// Calculate D[i], A_in, A_out
		for (int i = in; i <= out; i++) {
			double* f2 = new double[order];
			double* curphi = new double[order];
			for (int k = qin; k <= qout; k++) {
				curphi[k] = phi_p[k][i];
				f2[k] = root[k] * root[k] * curphi[k];
			}
			double integral0 = 2 * quadrature.integrateFromTo(order, qin, qout, curphi);
			double integral2 = 2 * quadrature.integrateFromTo(order, qin, qout, f2);
			D[i] = integral2 / integral0;
			if (i == in) {
				double* f1 = new double[order];
				for (int k = qin; k <= qout; k++) { f1[k] = abs(root[k]) * curphi[k]; }
				double integral1 = 2 * quadrature.integrateFromTo(order, qin, qout, f1);
				A_in = integral1 / integral2;
			}
			if (i == out) {
				double* f1 = new double[order];
				for (int k = qin; k <= qout; k++) { f1[k] = abs(root[k]) * curphi[k]; }
				double integral1 = 2 * quadrature.integrateFromTo(order, qin, qout, f1);
				A_out = integral1 / integral2;
			}
		}
		writeArrayEndl(D, in, out, &fout);
		fout << A_in << " " << A_out << endl << "~~~~~~~~~~~" << endl;
	}
}