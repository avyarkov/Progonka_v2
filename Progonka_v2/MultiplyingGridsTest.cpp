#include <fstream>
#include <string>
#include "FileIO.h"
#include "PlaneCase.h"
#include "Progonka.h"
#include "DataGrid.h"
#include "MyPlot.h"

void runMultiplyingGridsTest() {
	std::ifstream f_in("MultiplyingGridsTest_input.txt");
	std::ofstream f_out("MultiplyingGridsTest_output.txt");
	std::string curToken;

	int num;
	f_in >> num;
	DataGrid dg0 = DataGrid(num); // g = 0 : plane case
	int in = dg0.in, out = dg0.out;
	f_in >> curToken;
	readArray(dg0.r, in, out, &f_in);
	f_in >> curToken;
	readArray(dg0.sigma_0, in + 1, out, &f_in);
	f_in >> curToken;
	readArray(dg0.sigma_1, in + 1, out, &f_in);
	dg0.updateSigma();
	dg0.updateDtau();
	f_in >> curToken;
	readArray(dg0.T_P, in + 1, out, &f_in);
	f_in >> curToken;
	readArray(dg0.T_M, in + 1, out, &f_in);
	f_in >> curToken;
	f_in >> dg0.B_in >> dg0.B_out;
	f_in.close();

	int numIter = 30;
	for (int j = 1; j <= numIter; j++) {
		f_out << "j = " << j << ":\n";
		DataGrid dg = dg0.multiplied(j);
		dg.print(&f_out);
		Coeff coeff = getPlaneCoefficients(dg);
		coeff.print(&f_out);
		int size = dg.size, in = dg.in, out = dg.out;
		double* res_P = new double[size];	// MyDataTest results: Phi_P, Phi_M
		double* res_M = new double[size];
		Progonka(in, out, coeff, res_P, res_M);
		f_out << "result Phi_P:" << "\n";
		writeArray(res_P, in, out, &f_out);
		f_out << "\n" << "result Phi_M:" << "\n";
		writeArray(res_M, in, out, &f_out);
		plot(range(in, out), res_P, in, out, "plot" + std::to_string(j) + ".png");
		f_out << "\n~~~~~~~\n";
	}

	f_out.close();

}