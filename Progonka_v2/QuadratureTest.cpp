#include <iostream>
#include <stdio.h>
#include "Quadrature.h"
#include "FileIO.h"

using namespace std;

double f(double x) { return x * x; }

void runQuadratureTest() {
	Quadrature quadrature = Quadrature("Quadrature_Gauss–Legendre.txt");
	for (int i = quadrature.minOrder; i <= quadrature.maxOrder; i++) {
		int order = i;
		double* root = quadrature.rootsOfOrder(order);
		double* weight = quadrature.weightsOfOrder(order);
		printf("order = %d\n", order);
		writeArrayEndl(root, 0, order - 1, &cout);
		writeArrayEndl(weight, 0, order - 1, &cout);
		double* funcVal = new double[order];
		for (int i = 0; i < order; i++) {
			funcVal[i] = f(root[i]);
		}
		double result = quadrature.integrate(order, funcVal); // expecting result = 2/3 = 0.6666667
		printf("result = %0.15lf", result);
		cout << "\n-----------\n";
	}
}