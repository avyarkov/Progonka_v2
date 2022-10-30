#include <iostream>
#include <cstdlib>
#include "DataGrid.h"
#include "FileIO.h"
#include "Coeff.h"
#include "PlaneCase.h"
#include "Progonka.h"

using namespace std;

double* countGamma(double* D, int in, int out) {
	double* gamma = new double[DataGrid::sizeFromOut(out)];
	for (int i = in; i <= out; i++) {
		gamma[i] = 3 - 1 / D[i];
	}
	return gamma;
}

void runQuasidiffusion() {
	DataGrid quasiDG = readDataGrid("QuasidiffusionDataGrid.txt");
	DataGrid eoDG = readDataGrid("EvenOddDataGrid.txt");
	int num = quasiDG.num;
	if (num != eoDG.num) {
		cout << "Incorrect input: quasiDG and eoDg are of different sizes!";
		exit(239);
	}
	int size = quasiDG.size, in = quasiDG.in, out = quasiDG.out;
	double* D = new double[size];
	for (int i = in; i <= out; i++) {
		D[i] = 1.0 / 3;
	}
	double* gamma = countGamma(D, in, out);
	double A_in = 3.0 / 2, A_out = 3.0 / 2;
	Coeff quasiCoeff = getPlaneCoefficients(quasiDG);
	int numberOfIterations = 5;
	for (int it = 0; it < numberOfIterations; it++) {
		// solve QuasidiffusionEquations
		quasiCoeff.updateBoundaries(A_in, A_out);
		double* Psi_0 = new double[size], * Psi_1 = new double[size];
		Progonka(in, out, quasiCoeff, Psi_0, Psi_1);
		quasiCoeff.print(&cout);
		writeArray(Psi_0, in, out, &cout);
		cout << endl;
		writeArray(Psi_1, in, out, &cout);
		cout << endl;
	}
}