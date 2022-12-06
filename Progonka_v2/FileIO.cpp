#include <iostream>
#include <fstream>
#include <string>
#include "OldDataGrid.h"
#include "FileIO.h"

// NB! Boundaries: [inclusive, inclusive]
void readArray(double* arr, int in, int out, std::istream* f_in) {
	for (int i = in; i <= out; i++) {
		*f_in >> arr[i];
	}
}

void readTokenArray(std::string* token, double* arr, int in, int out, std::istream* f_in) {
	*f_in >> *token;
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

DataGrid readDataGridCellSectionsNodeSources(std::istream* f_in) {
	std::string curToken;
	*f_in >> curToken;
	int num = std::stoi(curToken);
	DataGrid dg = DataGrid(num);
	int in = dg.in, out = dg.out, size = dg.size;
	readTokenArray(&curToken, dg.r, in , out, f_in);
	double* sigma0 = new double[size], * sigma1 = new double[size];
	readTokenArray(&curToken, sigma0, in + 1, out, f_in);
	readTokenArray(&curToken, sigma1, in + 1, out, f_in);
	for (int i = in + 1; i <= out; i++) {
		dg.sigma0_r[i - 1] = sigma0[i];
		dg.sigma0_l[i] = sigma0[i];
		dg.sigma1_r[i - 1] = sigma1[i];
		dg.sigma1_l[i] = sigma1[i];
	}
	dg.updateSigma();
	dg.updateDtau();
	double* TP = new double[size], * TM = new double[size];
	readTokenArray(&curToken, TP, in, out, f_in);
	readTokenArray(&curToken, TM, in, out, f_in);
	for (int i = in; i <= out; i++) {
		dg.TP_l[i] = TP[i];
		dg.TP_r[i] = TP[i];
		dg.TM_l[i] = TM[i];
		dg.TM_r[i] = TM[i];
	}
	*f_in >> curToken >> dg.B_in >> dg.B_out;
	return dg;
}

DataGrid readDataGridCellSectionsNodeSources(std::string fileName) {
	std::ifstream f_in(fileName);
	DataGrid dg = readDataGridCellSectionsNodeSources(&f_in);
	f_in.close();
	return dg;
}

DataGrid readDataGridCellSectionsCellSources(std::istream* f_in) {
	std::string curToken;
	*f_in >> curToken;
	int num = std::stoi(curToken);
	DataGrid dg = DataGrid(num);
	int in = dg.in, out = dg.out, size = dg.size;
	readTokenArray(&curToken, dg.r, in, out, f_in);
	double* sigma0 = new double[size], * sigma1 = new double[size];
	readTokenArray(&curToken, sigma0, in + 1, out, f_in);
	readTokenArray(&curToken, sigma1, in + 1, out, f_in);
	for (int i = in + 1; i <= out; i++) {
		dg.sigma0_r[i - 1] = sigma0[i];
		dg.sigma0_l[i] = sigma0[i];
		dg.sigma1_r[i - 1] = sigma1[i];
		dg.sigma1_l[i] = sigma1[i];
	}
	dg.updateSigma();
	dg.updateDtau();
	double* TP = new double[size], * TM = new double[size];
	readTokenArray(&curToken, TP, in + 1, out, f_in);
	readTokenArray(&curToken, TM, in + 1, out, f_in);
	for (int i = in + 1; i <= out; i++) {
		dg.TP_r[i - 1] = TP[i];
		dg.TP_l[i] = TP[i];
		dg.TM_r[i - 1] = TM[i];
		dg.TM_l[i] = TM[i];
	}
	*f_in >> curToken >> dg.B_in >> dg.B_out;
	return dg;
}

DataGrid readDataGridCellSectionsCellSources(std::string fileName) {
	std::ifstream f_in(fileName);
	DataGrid dg = readDataGridCellSectionsCellSources(&f_in);
	f_in.close();
	return dg;
}