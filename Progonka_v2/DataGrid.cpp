#include <fstream>
#include <cmath>
#include "DataGrid.h"
#include "FileIO.h"

DataGrid::DataGrid(int num) {
	this->num = num;
	out = in + num;
	size = out + 2;
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

void DataGrid::print(std::ofstream* f_out) {
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
	if (sigma == nullptr) {
		sigma = new double[size];
		for (int i = in + 1; i <= out; i++) {
			sigma[i] = sqrt(abs(sigma_0[i]) * sigma_1[i]);
		}
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
