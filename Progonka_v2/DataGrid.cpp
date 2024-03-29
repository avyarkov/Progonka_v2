#include <iostream>
#include <cmath>
#include "DataGrid.h"
#include "FileIO.h"

DataGrid::DataGrid(int num) {
	this->num = num;
	out = in + num;
	size = out + 2;
	r = new double[size];
	dtau = new double[size];
	sigma0_l = new double[size];
	sigma0_r = new double[size];
	sigma1_l = new double[size];
	sigma1_r = new double[size];
	sigma_l = new double[size];
	sigma_r = new double[size];
	TP_l = new double[size];
	TP_r = new double[size];
	TM_l = new double[size];
	TM_r = new double[size];
	B_in = -566;
	B_out = -566;
}

void DataGrid::print(std::ostream* f_out) {
	*f_out << "DataGrid {\n";
	*f_out << "num = " << num << ", in = " << in << ", out = " << out;
	*f_out << "\n" << "r:\t";
	writeArray(r, in, out, f_out);
	*f_out << "\n" << "dtau:\t";
	writeArray(dtau, in + 1, out, f_out);
	*f_out << "\n" << "sigma0_l:\t";
	writeArray(sigma0_l, in, out, f_out);
	*f_out << "\n" << "sigma0_r:\t";
	writeArray(sigma0_r, in, out, f_out);
	*f_out << "\n" << "sigma1_l;:\t";
	writeArray(sigma1_l, in, out, f_out);
	*f_out << "\n" << "sigma1_r:\t";
	writeArray(sigma1_r, in, out, f_out);
	updateSigma();
	*f_out << "\n" << "sigma_l:\t";
	writeArray(sigma_l, in, out, f_out);
	*f_out << "\n" << "sigma_r:\t";
	writeArray(sigma_r, in, out, f_out);
	*f_out << "\n" << "TP_l:\t";
	writeArray(TP_l, in, out, f_out);
	*f_out << "\n" << "TP_r:\t";
	writeArray(TP_r, in, out, f_out);
	*f_out << "\n" << "TM_l:\t";
	writeArray(TM_l, in, out, f_out);
	*f_out << "\n" << "TM_r:\t";
	writeArray(TM_r, in, out, f_out);
	*f_out << "\n" << "B_in=" << B_in << ", B_out=" << B_out;
	*f_out << "\n}\n";
}

void DataGrid::updateSigma() {
	sigma_l = new double[size];
	sigma_r = new double[size];
	sigma_r[in] = sqrt(abs(sigma0_r[in]) * sigma1_r[in]);
	sigma_l[out] = sqrt(abs(sigma0_l[out]) * sigma1_l[out]);
	for (int i = in + 1; i <= out - 1; i++) {
		sigma_l[i] = sqrt(abs(sigma0_l[i]) * sigma1_l[i]);
		sigma_r[i] = sqrt(abs(sigma0_r[i]) * sigma1_r[i]);
	}
}

void DataGrid::updateDtau() {
	updateSigma();
	for (int i = in + 1; i <= out; i++) {
		dtau[i] = (sigma_r[i - 1] + sigma_l[i]) / 2 * (r[i] - r[i - 1]);
	}
}

DataGrid DataGrid::clone() {
	DataGrid res = DataGrid(num);
	for (int i = in; i <= out; i++) {
		res.r[i] = r[i];
		res.dtau[i] = dtau[i];
		res.sigma0_l[i] = sigma0_l[i];
		res.sigma0_r[i] = sigma0_r[i];
		res.sigma1_l[i] = sigma1_l[i];
		res.sigma1_r[i] = sigma1_r[i];
		res.sigma_l[i] = sigma_l[i];
		res.sigma_r[i] = sigma_r[i];
		res.TP_l[i] = TP_l[i];
		res.TP_r[i] = TP_r[i];
		res.TM_l[i] = TM_l[i];
		res.TM_r[i] = TM_r[i];
	}
	res.B_in = B_in;
	res.B_out = B_out;
	return res;
}

double f(double yl, double yr, int m, int k) {
	return yl + (yr - yl) / k * m;
}

DataGrid DataGrid::multiplied(int k) {
	updateSigma();
	int resNum = num * k;
	DataGrid res = DataGrid(resNum);
	int start = res.in;
	for (int i = start; i < res.out; i++) {
		int j = in + (i - start) / k;
		int m = (i - start) % k;
		res.r[i] = f(r[j], r[j + 1], m, k);
	}
	res.r[res.out] = r[out];
	for (int i = start; i <= res.out; i++) {
		if ((i - start) % k == 0) {
			int j = in + (i - start) / k;
			res.sigma0_l[i] = sigma0_l[j];
			res.sigma0_r[i] = sigma0_r[j];
			res.sigma1_l[i] = sigma1_l[j];
			res.sigma1_r[i] = sigma1_r[j];
			res.sigma_l[i] = sigma_l[j];
			res.sigma_r[i] = sigma_r[j];
			res.TP_l[i] = TP_l[j];
			res.TP_r[i] = TP_r[j];
			res.TM_l[i] = TM_l[j];
			res.TM_r[i] = TM_r[j];
		} else {
			int j = in + (i - start) / k;
			int m = (i - start) % k;
			res.sigma0_l[i] = f(sigma0_r[j], sigma0_l[j + 1], m, k);
			res.sigma0_r[i] = res.sigma0_l[i];
			res.sigma1_l[i] = f(sigma1_r[j], sigma1_l[j + 1], m, k);
			res.sigma1_r[i] = res.sigma1_l[i];
			res.sigma_l[i] = f(sigma_r[j], sigma_l[j + 1], m, k);
			res.sigma_r[i] = res.sigma_l[i];
			res.TP_l[i] = f(TP_r[j], TP_l[j + 1], m, k);
			res.TP_r[i] = res.TP_l[i];
			res.TM_l[i] = f(TM_r[j], TM_l[j + 1], m, k);
			res.TM_r[i] = res.TM_l[i];
		}
	}
	res.updateDtau();
	res.B_in = B_in;
	res.B_out = B_out;
	return res;
}

DataGrid DataGrid::withAddedSources(double* AddP_l, double* AddP_r, double* AddM_l, double* AddM_r) {
	updateSigma();
	updateDtau();
	DataGrid res = this->clone();
	res.TP_r[in] += AddP_r[in];
	res.TM_r[in] += AddM_r[in];
	res.TP_l[out] += AddP_l[out];
	res.TM_l[out] += AddM_l[out];
	for (int i = in + 1; i <= out - 1; i++) {
		res.TP_l[i] += AddP_l[i];
		res.TP_r[i] += AddP_r[i];
		res.TM_l[i] += AddM_l[i];
		res.TM_r[i] += AddM_r[i];
	}
	return res;
}

DataGrid DataGrid::withAddedSources(double* Add_P, double* Add_M) {
	return DataGrid::withAddedSources(Add_P, Add_P, Add_M, Add_M);
}

DataGrid DataGrid::toPlaneCharacteristics(double mu) {
	if (mu < 0) { std::cout << "negative mu!"; exit(1); }
	updateSigma();
	updateDtau();
	DataGrid res = this->clone();
	res.r[in] = r[in];
	for (int i = in + 1; i <= out; i++) {
		res.r[i] = res.r[i - 1] + (r[i] - r[i - 1]) / mu;
	}
	res.updateDtau();
	return res;
}

DataGrid DataGrid::withDividedSigma0(double* D) {
	updateSigma();
	updateDtau();
	DataGrid res = this->clone();
	res.sigma0_r[in] /= D[in];
	res.sigma0_l[out] /= D[out];
	for (int i = in + 1; i <= out - 1; i++) {
		res.sigma0_r[i] /= D[i];
		res.sigma0_l[i] /= D[i];
	}
	res.updateSigma();
	res.updateDtau();
	return res;
}

