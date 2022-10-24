#pragma once
#include "DataGrid.h"
#include "Coeff.h"

Coeff getSphericalCoefficients(DataGrid dg, double* gamma);

void leftSphericalSupplement(DataGrid dg, double* gamma, Coeff coeff, double* PHI_P, double* PHI_M);