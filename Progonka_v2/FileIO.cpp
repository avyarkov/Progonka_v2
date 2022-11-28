#include <iostream>
#include <fstream>
#include <string>
#include "OldDataGrid.h"

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

OldDataGrid readOldDataGrid(std::istream* f_in) {
	int num;
	*f_in >> num;
	OldDataGrid dg = OldDataGrid(num);
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

OldDataGrid readOldDataGrid(std::string fileName) {
	std::ifstream f_in(fileName);
	OldDataGrid dg = readOldDataGrid(&f_in);
	f_in.close();
	return dg;
}

OldDataGridNodeSources readOldDataGridNodeSources(std::istream* f_in) {
	int num;
	*f_in >> num;
	OldDataGridNodeSources dg = OldDataGridNodeSources(num);
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
	readArray(dg.T_P, in, out, f_in);
	*f_in >> curToken;
	readArray(dg.T_M, in, out, f_in);
	*f_in >> curToken;
	*f_in >> dg.B_in >> dg.B_out;
	return dg;
}

OldDataGridNodeSources readOldDataGridNodeSources(std::string fileName) {
	std::ifstream f_in(fileName);
	OldDataGridNodeSources dg = readOldDataGridNodeSources(&f_in);
	f_in.close();
	return dg;
}