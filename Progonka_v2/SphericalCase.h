#pragma once
#include "OldDataGrid.h"
#include "Coeff.h"

Coeff getSphericalCoefficients(OldDataGrid dg, double* gamma);

void leftSphericalSupplement(OldDataGrid dg, double* gamma, Coeff coeff, double* PHI_P, double* PHI_M);