#include <iostream>
#include "Coeff.h"
#include "FileIO.h"

Coeff::Coeff(int num) {
	this->num = num;
	out = in + num;
	size = out + 2;
	el = new double[size];
	er = new double[size];
	kl = new double[size];
	kr = new double[size];
	Rr = new double[size];
	Rl = new double[size];
}

void Coeff::print(std::ostream* f_out) {
	*f_out << "Coeff {\n";
	writeArray(el, in, out + 1, f_out);
	*f_out << "\n";
	writeArray(er, in, out + 1, f_out);
	*f_out << "\n";
	writeArray(kl, in, out + 1, f_out);
	*f_out << "\n";
	writeArray(kr, in, out + 1, f_out);
	*f_out << "\n";
	writeArray(Rl, in, out + 1, f_out);
	*f_out << "\n";
	writeArray(Rr, in, out + 1, f_out);
	*f_out << "\n}\n";
}

void Coeff::updateBoundaries(double A_in, double A_out) {
	kl[in] = A_in;
	kr[out + 1] = A_out;
}
