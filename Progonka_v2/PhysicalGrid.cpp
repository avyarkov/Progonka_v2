#include "PhysicalGrid.h"
#include "FileIO.h"

PhysicalGrid::PhysicalGrid(int num) {
	this->num = num;
	out = in + num;
	size = out + 2;
	SigmaA_l= new double[size];
	SigmaA_r = new double[size];
	SigmaS_l = new double[size];
	SigmaS_r = new double[size];
	Q_l = new double[size];
	Q_r = new double[size];
	B_in = -566;
	B_out = -566;
}

void PhysicalGrid::print(std::ostream* f_out) {
	*f_out << "PhysicalGrid {\n";
	*f_out << "num = " << num << ", in = " << in << ", out = " << out;
	*f_out << "\n" << "SigmaS_l:\t";
	writeArray(SigmaS_l, in, out, f_out);
	*f_out << "\n" << "SigmaS_r:\t";
	writeArray(SigmaS_r, in, out, f_out);
	*f_out << "\n" << "SigmaA_l:\t";
	writeArray(SigmaA_l, in, out, f_out);
	*f_out << "\n" << "SigmaA_r:\t";
	writeArray(SigmaA_r, in, out, f_out);
	*f_out << "\n" << "Q_l:\t";
	writeArray(Q_l, in, out, f_out);
	*f_out << "\n" << "Q_r:\t";
	writeArray(Q_r, in, out, f_out);
	*f_out << "\n" << "B_in = " << B_in << ", B_out = " << B_out;
	*f_out << "\n}\n";
}
