#include <fstream>

// NB! Boundaries: [inclusive, inclusive]
void readArray(double* arr, int in, int out, std::ifstream* f_in) {
	for (int i = in; i <= out; i++) {
		*f_in >> arr[i];
	}
}

void writeArray(double* arr, int in, int out, std::ofstream* f_out) {
	for (int i = in; i <= out; i++) {
		*f_out << arr[i] << " ";
	}
}

void writeArrayWithTabs(double* arr, int in, int out, std::ofstream* f_out) {
	for (int i = in; i <= out; i++) {
		*f_out << arr[i] << "\t";
	}
}