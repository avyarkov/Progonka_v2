#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include "FileIO.h"
#include "Progonka.h"
#include "PlaneCase.h"
#include "SphericalCase.h"
#include "Tests.h"

int main() {

	//runDataTestPlane();

	//runDataTestSpherical();

	//runExactTest();

	//runEquationTest();

	//runMultiplyingGridsTest();

	runMultiplyingGridsTestSpherical();

	std::ifstream f_in("DataGridTest_input.txt");
	std::ofstream f_out("DataGridTest_output.txt");
	std::string curToken;

	int num;
	f_in >> num;
	DataGrid dg = DataGrid(num); // g = 0 : plane case
	int in = dg.in, out = dg.out;
	f_in >> curToken;
	readArray(dg.r, in, out, &f_in);
	f_in >> curToken;
	readArray(dg.sigma_0, in + 1, out, &f_in);
	f_in >> curToken;
	readArray(dg.sigma_1, in + 1, out, &f_in);
	dg.updateSigma();
	dg.updateDtau();
	f_in >> curToken;
	readArray(dg.T_P, in + 1, out, &f_in);
	f_in >> curToken;
	readArray(dg.T_M, in + 1, out, &f_in);
	f_in >> curToken;
	f_in >> dg.B_in >> dg.B_out;
	f_in.close();

	dg.print(&f_out);
	dg.multiplied(100).print(&f_out);
}