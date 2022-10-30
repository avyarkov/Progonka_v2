#include <fstream>
#include <string>
#include "FileIO.h"
#include "Progonka.h"
#include "DataGrid.h"
#include "MyPlot.h"
#include "SphericalCase.h"

void runMultiplyingGridsTestSpherical() {
	std::ifstream f_in("MultiplyingGridsTestSpherical_input.txt");
	std::ofstream f_out("MultiplyingGridsTestSpherical_output.txt");
	std::string curToken;

	int num;
	f_in >> num;
	DataGrid dg0 = DataGrid(num);
	int in = dg0.in, out = dg0.out, size = dg0.size;
	f_in >> curToken;
	readArray(dg0.dtau, in + 1, out, &f_in);
	f_in >> curToken;
	readArray(dg0.sigma_0, in + 1, out, &f_in);
	f_in >> curToken;
	readArray(dg0.sigma_1, in + 1, out, &f_in);
	dg0.updateSigma();
	f_in >> curToken;
	readArray(dg0.T_P, in + 1, out, &f_in);
	f_in >> curToken;
	readArray(dg0.T_M, in + 1, out, &f_in);
	f_in >> curToken;
	readArray(dg0.r, in, out, &f_in);
	f_in >> curToken;
	f_in >> dg0.B_in >> dg0.B_out;
	f_in.close();

	int numIter = 15;
	for (int j = 1; j <= numIter; j++) {
		f_out << "j = " << j << ":\n";
		DataGrid dg = dg0.multiplied(j);
		int size = dg.size, in = dg.in, out = dg.out;
		dg.print(&f_out);
		double* gamma = new double[size];
		for (int i = 0; i < size; i++) {
			gamma[i] = 0;
		}
		Coeff coeff = getSphericalCoefficients(dg, gamma);
		coeff.print(&f_out);
		double* res_P = new double[size];	// MyDataTest results: Phi_P, Phi_M
		double* res_M = new double[size];
		if (dg.r[in] == 0) {
			Progonka(in + 1, out, coeff, res_P, res_M);
			leftSphericalSupplement(dg, gamma, coeff, res_P, res_M);
		}
		else {
			coeff.el[in] = 0;
			dg.B_in = 0;
			coeff.kl[in] = 0;
			Progonka(in, out, coeff, res_P, res_M);
		}
		f_out << "result Phi_P:" << "\n";
		writeArray(res_P, in, out, &f_out);
		f_out << "\n" << "result Phi_M:" << "\n";
		writeArray(res_M, in, out, &f_out);
		plot(dg.r, res_P, in, out, "plot" + std::to_string(j) + ".png");
		plot(dg.r, res_P, in, in + dg.num / 5, "plotBeg" + std::to_string(j) + ".png");
		f_out << "\n~~~~~~~\n";
	}

	f_out.close();

}