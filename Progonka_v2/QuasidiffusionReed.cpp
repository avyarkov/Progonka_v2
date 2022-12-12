#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "OldDataGrid.h"
#include "FileIO.h"
#include "MyPlot.h"
#include "Coeff.h"
#include "PlaneCase.h"
#include "Progonka.h"
#include "Quadrature.h"
#include "Quasidiffusion.h"

using namespace std;

void runQuasidiffusion_Reed() {
	ofstream fout("Quasidiffusion_output_Reed.txt");

	const double w_0 = 1.0, w_1 = 0.5, w_2 = 0;
	const double nu = 1;
	const int mult = 20;
	DataGrid quasiDG = readDataGridCellSectionsCellSources("QuasidiffusionDataGrid_Reed.txt").multiplied(mult);
	quasiDG.updateSigma();
	quasiDG.updateDtau();
	DataGrid eoDG = readDataGridCellSectionsCellSources("EvenOddDataGrid_Reed.txt").multiplied(mult);
	eoDG.updateSigma();
	eoDG.updateDtau();

	//quasiDG.print(&cout);
	//eoDG.print(&cout);

	int num = quasiDG.num;
	if (num != eoDG.num) { cout << "Incorrect input: quasiDG and eoDg are of different sizes!"; exit(239); }
	int size = quasiDG.size, in = quasiDG.in, out = quasiDG.out;
	double* Sigma_s_l = new double[size], * Sigma_s_r = new double[size];
	for (int i = in; i <= out; i++) {
		Sigma_s_l[i] = quasiDG.sigma1_l[i] - quasiDG.sigma0_l[i];
		Sigma_s_r[i] = quasiDG.sigma1_r[i] - quasiDG.sigma0_r[i];
	}

	double* D = new double[size];
	for (int i = in; i <= out; i++) { D[i] = 1.0 / 3; }
	//double* gamma = countGamma(D, in, out);
	double A_in = 3.0 / 2, A_out = 3.0 / 2;

	Quadrature quadrature = Quadrature("Quadrature_Gauss–Legendre.txt");
	int order = 6;
	int qin = order / 2, qout = order - 1;
	double* root = quadrature.rootsOfOrder(order);
	double* weight = quadrature.weightsOfOrder(order);

	int numberOfIterations = 40;
	double** res_PSI_0 = new double* [numberOfIterations];
	for (int it = 0; it < numberOfIterations; it++) {
		// solve Quasidiffusion equations
		DataGrid curQuasiDg = quasiDG.withDividedSigma0(D);
		Coeff quasiCoeff = getPlaneCoefficients(curQuasiDg);
		quasiCoeff.updateBoundaries(A_in, A_out);
		//quasiCoeff.print(&cout);

		double* PSI_0 = new double[size], * PSI_1 = new double[size];
		Progonka(in, out, quasiCoeff, PSI_0, PSI_1);
		for (int i = in; i <= out; i++) { PSI_0[i] /= D[i]; }

		writeArrayEndl(PSI_0, in, out, &fout);
		res_PSI_0[it] = new double[size];
		for (int i = in; i <= out; i++) {
			res_PSI_0[it][i] = PSI_0[i];
		}
		//plot(curQuasiDg.r, PSI_0, in, out, "plot_Reed_" + std::to_string(it) + ".png");

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
			double* AddP_l = new double[size], * AddP_r = new double[size], * AddM_l = new double[size], * AddM_r = new double[size];
			for (int i = in; i <= out; i++) {
				AddP_l[i] = nu * Sigma_s_l[i] * F_P[k][i];
				AddP_r[i] = nu * Sigma_s_r[i] * F_P[k][i];
				AddM_l[i] = nu * Sigma_s_l[i] * F_M[k][i];
				AddM_r[i] = nu * Sigma_s_r[i] * F_M[k][i];
			}
			DataGrid curEODG = eoDG.toPlaneCharacteristics(root[k]).withAddedSources(AddP_l, AddP_r, AddM_l, AddM_r);
			Coeff curCoeff = getPlaneCoefficients(curEODG);
			curCoeff.updateBoundaries(1.0, 1.0);
			Progonka(in, out, curCoeff, phi_p[k], phi_m[k]);
			//writeArrayEndl(phi_p[k], in, out, &fout);
			//writeArrayEndl(phi_m[k], in, out, &fout);
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
		fout << "D: ";
		writeArrayEndl(D, in, out, &fout);
		fout << "A_in/A_out: " << A_in << " " << A_out << endl << "~~~~~~~~~~~" << endl;
	}

	int lastIteration = numberOfIterations - 1;
	double* measures = new double[numberOfIterations], * logMeasures = new double[numberOfIterations];
	for (int it = 0; it < numberOfIterations; it++) {
		for (int i = in; i <= out; i++) {
			double curMeasure = measure(in, out, res_PSI_0[it], res_PSI_0[lastIteration]);
			measures[it] = curMeasure;
			logMeasures[it] = log10(curMeasure);
		}
	}
	fout << "/////////////////" << endl;
	fout << "res_PSI_0 @ lastIteration: ";  writeArrayEndl(res_PSI_0[lastIteration], in, out, &fout);
	fout << "measures: "; writeArrayEndl(measures, 0, lastIteration, &fout);
	fout << "logMeasures: "; writeArrayEndl(logMeasures, 0, lastIteration, &fout);
	plot(range(0, lastIteration), measures, 0, lastIteration, "plot_Reed_measures.png");
	plot(range(0, lastIteration), logMeasures, 0, lastIteration - 1, "plot_Reed_logMeasures.png");
}