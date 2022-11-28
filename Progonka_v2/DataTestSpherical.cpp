#include <fstream>
#include "FileIO.h"
#include "SphericalCase.h"
#include "Progonka.h"
#include "OldDataGrid.h"
#include "MyPlot.h"

void runDataTestSpherical() {
	std::ifstream f_in("DataTestSpherical_input.txt");
	std::ofstream f_out("DataTestSpherical_output.txt");
	std::string curToken;

	int num;
	f_in >> num;
	OldDataGrid dg = OldDataGrid(num);
	int in = dg.in, out = dg.out, size = dg.size;
	double* gamma = new double[size];
	f_in >> curToken;
	readArray(dg.dtau, in + 1, out, &f_in);
	f_in >> curToken;
	readArray(dg.sigma_0, in + 1, out, &f_in);
	f_in >> curToken;
	readArray(dg.sigma_1, in + 1, out, &f_in);
	dg.updateSigma();
	f_in >> curToken;
	readArray(dg.T_P, in + 1, out, &f_in);
	f_in >> curToken;
	readArray(dg.T_M, in + 1, out, &f_in);
	f_in >> curToken;
	readArray(dg.r, in, out, &f_in);
	f_in >> curToken;
	readArray(gamma, in, out, &f_in);
	f_in >> curToken;
	f_in >> dg.B_in >> dg.B_out;
	f_in.close();

	Coeff coeff = getSphericalCoefficients(dg, gamma);

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
	f_out << "\n" << "-------" << "\n";
	f_out << "reference dtau:" << "\n";
	writeArray(dg.dtau, in + 1, out, &f_out);
	f_out << "\n" << "reference sigma_0:" << "\n";
	writeArray(dg.sigma_0, in + 1, out, &f_out);
	f_out << "\n" << "reference sigma_1:" << "\n";
	writeArray(dg.sigma_1, in + 1, out, &f_out);
	f_out << "\n" << "reference T_P:" << "\n";
	writeArray(dg.T_P, in + 1, out, &f_out);
	f_out << "\n" << "reference T_M:" << "\n";
	writeArray(dg.T_M, in + 1, out, &f_out);
	f_out << "\n" << "reference r:" << "\n";
	writeArray(dg.r, in, out, &f_out);
	f_out << "\n" << "reference gamma:" << "\n";
	writeArray(gamma, in, out, &f_out);
	f_out << "\n" << "reference B_in B_out:" << "\n" << dg.B_in << " " << dg.B_out << "\n";

	// Coefficients output
	bool ifPrintCoefficients = true;
	if (ifPrintCoefficients) {
		f_out << "-------\n";
		coeff.print(&f_out);
	}

	plot(range(in, out), res_P, in, out, "plot.png");

	f_out.close();

}