#pragma once
#include <string>
#include <fstream>

struct Quadrature {
	int minOrder, maxOrder;
	double** root;
	double** weight;

	Quadrature(std::string fileName);
	double* rootsOfOrder(int order);
	double* weightsOfOrder(int order);
	double integrate(int order, double* functionValues);
	double integrateFromTo(int order, int from, int to, double* functionValues);
};

