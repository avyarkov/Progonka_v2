#include <iostream>
#include <cmath>
#include "DataGrid.h"
#include "FileIO.h"

DataGrid::DataGrid() {
}

DataGrid::DataGrid(int num) {
	this->num = num;
	out = in + num;
	size = sizeFromOut(out);
	r = new double[size];
	dtau = new double[size];
	sigma_0 = new double[size];
	sigma_1 = new double[size];
	sigma = nullptr;
	T_P = new double[size];
	T_M = new double[size];
	B_in = -566;
	B_out = -566;
}

int DataGrid::sizeFromOut(int out) {
	return out + 2;
}

void DataGrid::print(std::ostream* f_out) {
	*f_out << "DataGrid {\n";
	*f_out << "num=" << num << ", in=" << in << ", out=" << out << ";";
	*f_out << "\n" << "r:\t";
	writeArray(r, in, out, f_out);
	*f_out << "\n" << "dtau:\t";
	writeArray(dtau, in + 1, out, f_out);
	*f_out << "\n" << "sigma_0:\t";
	writeArray(sigma_0, in + 1, out, f_out);
	*f_out << "\n" << "sigma_1:\t";
	writeArray(sigma_1, in + 1, out, f_out);
	*f_out << "\n" << "sigma:\t";
	updateSigma();
	writeArray(sigma, in + 1, out, f_out);
	*f_out << "\n" << "T_P:\t";
	writeArray(T_P, in + 1, out, f_out);
	*f_out << "\n" << "T_M:\t";
	writeArray(T_M, in + 1, out, f_out);
	*f_out << "\n" << "B_in=" << B_in << ", B_out=" << B_out;
	*f_out << "\n}\n";
}

void DataGrid::updateSigma() {
	sigma = new double[size];
	for (int i = in + 1; i <= out; i++) {
		sigma[i] = sqrt(abs(sigma_0[i]) * sigma_1[i]);
	}
}

void DataGrid::updateDtau() {
	updateSigma();
	for (int i = in + 1; i <= out; i++) {
		dtau[i] = sigma[i] * (r[i] - r[i - 1]);
	}
}

DataGrid DataGrid::multiplied(int k) {
	updateSigma();
	int resNum = num * k;
	DataGrid res = DataGrid(resNum);
	res.sigma = new double[res.size];
	int start = res.in + 1;
	for (int i = start; i <= res.out; i++) {
		int j = (in + 1) + (i - start) / k;
		res.sigma_0[i] = sigma_0[j];
		res.sigma_1[i] = sigma_1[j];
		res.sigma[i] = sigma[j];
		res.T_P[i] = T_P[j];
		res.T_M[i] = T_M[j];
	}
	res.B_in = B_in;
	res.B_out = B_out;
	start = res.in;
	for (int i = start; i < res.out; i++) {
		int j = in + (i - start) / k;
		res.r[i] = r[j] + (r[j + 1] - r[j]) / k * ((i - start) % k);
	}
	res.r[res.out] = r[out];
	res.updateDtau();
	return res;
}

DataGrid DataGrid::withAddedSources(double* Add_P, double* Add_M) {
	updateSigma();
	updateDtau();
	DataGrid res = DataGrid(num);
	res.sigma = new double[res.size];
	for (int i = in; i <= out; i++) {
		res.r[i] = r[i];
		res.dtau[i] = dtau[i];
		res.sigma_0[i] = sigma_0[i];
		res.sigma_1[i] = sigma_1[i];
		res.sigma[i] = sigma[i];
		res.T_P[i] = T_P[i] + Add_P[i];
		res.T_M[i] = T_M[i] + Add_M[i];
	}
	res.B_in = B_in;
	res.B_out = B_out;
	return res;
}

DataGrid DataGrid::toCharacteristics(double mu) {
	if (mu < 0) { std::cout << "negative mu!"; exit(1); }
	updateSigma();
	updateDtau();
	DataGrid res = DataGrid(num);
	res.r[in] = r[in];
	res.T_P[in] = T_P[in];
	res.T_M[in] = T_M[in];
	for (int i = in + 1; i <= out; i++) {
		res.r[i] = res.r[i - 1] + (r[i] - r[i - 1]) / mu;
		res.sigma_0[i] = sigma_0[i];
		res.sigma_1[i] = sigma_1[i];
		res.T_P[i] = T_P[i];
		res.T_M[i] = T_M[i];
	}
	res.B_in = B_in;
	res.B_out = B_out;
	res.updateSigma();
	res.updateDtau();
	return res;
}

DataGrid DataGrid::withDividedSigma0(double* D) {
	updateSigma();
	updateDtau();
	DataGrid res = DataGrid(num);
	res.sigma = new double[res.size];
	for (int i = in; i <= out; i++) {
		res.r[i] = r[i];
		res.dtau[i] = dtau[i];
		res.sigma_0[i] = sigma_0[i] / D[i];
		res.sigma_1[i] = sigma_1[i];
		res.T_P[i] = T_P[i];
		res.T_M[i] = T_M[i];
	}
	res.updateSigma();
	res.updateDtau();
	res.B_in = B_in;
	res.B_out = B_out;
	return res;
}

DataGridNodeSources::DataGridNodeSources(int num) {
	this->num = num;
	out = in + num;
	size = sizeFromOut(out);
	r = new double[size];
	dtau = new double[size];
	sigma_0 = new double[size];
	sigma_1 = new double[size];
	sigma = nullptr;
	T_P = new double[size];
	T_M = new double[size];
	B_in = -566;
	B_out = -566;
}

void DataGridNodeSources::print(std::ostream* f_out) {
	*f_out << "DataGrid {\n";
	*f_out << "num=" << num << ", in=" << in << ", out=" << out << ";";
	*f_out << "\n" << "r:\t";
	writeArray(r, in, out, f_out);
	*f_out << "\n" << "dtau:\t";
	writeArray(dtau, in + 1, out, f_out);
	*f_out << "\n" << "sigma_0:\t";
	writeArray(sigma_0, in + 1, out, f_out);
	*f_out << "\n" << "sigma_1:\t";
	writeArray(sigma_1, in + 1, out, f_out);
	*f_out << "\n" << "sigma:\t";
	updateSigma();
	writeArray(sigma, in + 1, out, f_out);
	*f_out << "\n" << "T_P:\t";
	writeArray(T_P, in, out, f_out);
	*f_out << "\n" << "T_M:\t";
	writeArray(T_M, in, out, f_out);
	*f_out << "\n" << "B_in=" << B_in << ", B_out=" << B_out;
	*f_out << "\n}\n";
}

DataGridNodeSources DataGridNodeSources::multiplied(int k) {
	std::cout << "Can't multiply DataGridNodeSources";
	exit(1);
}

DataGridNodeSources DataGridNodeSources::withAddedSources(double* Add_P, double* Add_M) {
	updateSigma();
	updateDtau();
	DataGridNodeSources res = DataGridNodeSources(num);
	res.sigma = new double[res.size];
	for (int i = in; i <= out; i++) {
		res.r[i] = r[i];
		res.dtau[i] = dtau[i];
		res.sigma_0[i] = sigma_0[i];
		res.sigma_1[i] = sigma_1[i];
		res.sigma[i] = sigma[i];
		res.T_P[i] = T_P[i] + Add_P[i];
		res.T_M[i] = T_M[i] + Add_M[i];
	}
	res.B_in = B_in;
	res.B_out = B_out;
	return res;
}

DataGridNodeSources DataGridNodeSources::toCharacteristics(double mu) {
	if (mu < 0) { std::cout << "negative mu!"; exit(1); }
	updateSigma();
	updateDtau();
	DataGridNodeSources res = DataGridNodeSources(num);
	res.r[in] = r[in];
	res.T_P[in] = T_P[in];
	res.T_M[in] = T_M[in];
	for (int i = in + 1; i <= out; i++) {
		res.r[i] = res.r[i - 1] + (r[i] - r[i - 1]) / mu;
		res.sigma_0[i] = sigma_0[i];
		res.sigma_1[i] = sigma_1[i];
		res.T_P[i] = T_P[i];
		res.T_M[i] = T_M[i];
	}
	res.B_in = B_in;
	res.B_out = B_out;
	res.updateSigma();
	res.updateDtau();
	return res;
}

DataGridNodeSources DataGridNodeSources::withDividedSigma0(double* D) {
	updateSigma();
	updateDtau();
	DataGridNodeSources res = DataGridNodeSources(num);
	res.sigma = new double[res.size];
	for (int i = in; i <= out; i++) {
		res.r[i] = r[i];
		res.dtau[i] = dtau[i];
		res.sigma_0[i] = sigma_0[i] / D[i];
		res.sigma_1[i] = sigma_1[i];
		res.T_P[i] = T_P[i];
		res.T_M[i] = T_M[i];
	}
	res.updateSigma();
	res.updateDtau();
	res.B_in = B_in;
	res.B_out = B_out;
	return res;
}