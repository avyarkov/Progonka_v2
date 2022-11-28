#include <fstream>
#include "FileIO.h"
#include "PlaneCase.h"
#include "Progonka.h"

void runEquationTest() {
	std::ifstream f_in;
	std::ofstream f_out;
	std::string curToken;
	
	int num;
	f_in.open("EquationTest_input.txt");
	f_in >> num;
	OldDataGrid dg = OldDataGrid(num);
	int in = dg.in, out = dg.out, size = dg.size;
	double* Phi_P = new double[size], * Phi_M = new double[size];
	f_in >> curToken;
	readArray(dg.dtau, in + 1, out, &f_in);
	f_in >> curToken;
	readArray(dg.sigma_0, in + 1, out, &f_in);
	f_in >> curToken;
	readArray(dg.sigma_1, in + 1, out, &f_in);
	f_in >> curToken;
	readArray(Phi_P, in, out, &f_in);
	f_in >> curToken;
	readArray(Phi_M, in, out, &f_in);
	dg.updateSigma();
	for (int i = in + 1; i <= out; i++) {
		// T_P[i] = sigma[i] * (Phi_M[i] - Phi_M[i - 1]) / dtau[i] + sigma_0[i] * Phi_P[i];
		// T_M[i] = sigma[i] * (Phi_P[i] - Phi_P[i - 1]) / dtau[i] + sigma_1[i] * Phi_M[i];
		dg.T_P[i] = dg.sigma[i] * (Phi_M[i] - Phi_M[i - 1]) / dg.dtau[i] + dg.sigma_0[i] * (Phi_P[i - 1] + Phi_P[i]) / 2;
		dg.T_M[i] = dg.sigma[i] * (Phi_P[i] - Phi_P[i - 1]) / dg.dtau[i] + dg.sigma_1[i] * (Phi_M[i - 1] + Phi_M[i]) / 2;
	}
	f_in.close();
	dg.B_in = Phi_M[in] + 3.0 / 2 * Phi_P[in];
	dg.B_out = 3.0 / 2 * Phi_P[out] - Phi_M[out];

	Coeff coeff = getPlaneCoefficients(dg);

	double* et_res_P = new double[size];		// EquationTest results: Phi_P, Phi_M
	double* et_res_M = new double[size];
	Progonka(in, out, coeff, et_res_P, et_res_M);

	f_out.open("EquationTest_output.txt");
	f_out << num << "\n" << "results:\n";
	writeArrayWithTabs(et_res_P, in, out, &f_out);
	f_out << "\n";
	writeArrayWithTabs(et_res_M, in, out, &f_out);
	f_out << "\n" << "reference Phi:\n";
	writeArrayWithTabs(Phi_P, in, out, &f_out);
	f_out << "\n";
	writeArrayWithTabs(Phi_M, in, out, &f_out);
	f_out << "\n" << "reference sigma:\n";
	writeArrayWithTabs(dg.sigma_0, in + 1, out, &f_out);
	f_out << "\n";
	writeArrayWithTabs(dg.sigma_1, in + 1, out, &f_out);
	f_out.close();
}