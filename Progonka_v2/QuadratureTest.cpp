#include "Quadrature.h"
#include "FileIO.h"
#include <iostream>

double f(double x) { return x * x; }

void runQuadratureTest() {
	Quadrature quadrature = Quadrature("Quadrature_Gauss–Legendre.txt");
	int order = 5;
	double* root = quadrature.rootsOfOrder(order);
	double* weight = quadrature.weightsOfOrder(order);
	for (int i = 0; i < order; i++) {
		std::cout << root[i] << " " << weight[i] << "\n";
	}
	std::cout << "-------\n";
	double* funcVal = new double[order];
	for (int i = 0; i < order; i++) {
		funcVal[i] = f(root[i]);
	}
	std::cout << quadrature.integrate(order, funcVal); // expecting result = 2/3 = 0.6666667
}