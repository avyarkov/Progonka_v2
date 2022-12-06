#pragma once

void runQuasidiffusion();

double* countGamma(double* D, int in, int out);

double calcF_P(double mu, double D, double PSI_0, double w_0, double w_2);

double calcF_M(double mu, double PSI_1, double w_1);
