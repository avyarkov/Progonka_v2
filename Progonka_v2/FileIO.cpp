#include <iostream>
#include <fstream>
#include <string>
#include "DataGrid.h"

// NB! Boundaries: [inclusive, inclusive]
void readArray(double* arr, int in, int out, std::istream* f_in) {
	for (int i = in; i <= out; i++) {
		*f_in >> arr[i];
	}
}

void writeArray(double* arr, int in, int out, std::ostream* f_out) {
	for (int i = in; i <= out; i++) {
		*f_out << arr[i] << " ";
	}
}

void writeArrayEndl(double* arr, int in, int out, std::ostream* f_out) {
	writeArray(arr, in, out, f_out);
	*f_out << std::endl;
}

void writeArrayWithTabs(double* arr, int in, int out, std::ostream* f_out) {
	for (int i = in; i <= out; i++) {
		*f_out << arr[i] << "\t";
	}
}

DataGrid readDataGrid(std::istream* f_in) {
	int num;
	*f_in >> num;
	DataGrid dg = DataGrid(num);
	int in = dg.in, out = dg.out;
	std::string curToken;
	*f_in >> curToken;
	readArray(dg.r, in, out, f_in);
	*f_in >> curToken;
	readArray(dg.sigma_0, in + 1, out, f_in);
	*f_in >> curToken;
	readArray(dg.sigma_1, in + 1, out, f_in);
	dg.updateSigma();
	dg.updateDtau();
	*f_in >> curToken;
	readArray(dg.T_P, in + 1, out, f_in);
	*f_in >> curToken;
	readArray(dg.T_M, in + 1, out, f_in);
	*f_in >> curToken;
	*f_in >> dg.B_in >> dg.B_out;
	return dg;
}

DataGrid readDataGrid(std::string fileName) {
	std::ifstream f_in(fileName);
	DataGrid dg = readDataGrid(&f_in);
	f_in.close();
	return dg;
}