#include "Quadrature.h"
#include "FileIO.h"

using namespace std;

Quadrature::Quadrature(std::string fileName) {
	ifstream fin(fileName);
	string curToken;
	fin >> curToken >> minOrder >> maxOrder;
	int size = maxOrder + 1;
	root = new double* [size];
	weight = new double* [size];
	for (int i = minOrder; i <= maxOrder; i++) {
		fin >> curToken;
		root[i] = new double[i];
		weight[i] = new double[i];
		readArray(root[i], 0, i - 1, &fin);
		readArray(weight[i], 0, i - 1, &fin);
	}
}

double* Quadrature::rootsOfOrder(int order) {
	return root[order];
}

double* Quadrature::weightsOfOrder(int order) {
	return weight[order];
}

double Quadrature::integrate(int order, double* functionValues) {
	return integrateFromTo(order, 0, order - 1, functionValues);
}

double Quadrature::integrateFromTo(int order, int from, int to, double* functionValues) {
	double res = 0;
	for (int i = from; i <= to; i++) {
		res += weight[order][i] * functionValues[i];
	}
	return res;
}
