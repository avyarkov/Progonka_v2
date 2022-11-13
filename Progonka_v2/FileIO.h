#pragma once
#include <iostream>
#include <string>
#include "DataGrid.h"

void readArray(double* arr, int in, int out, std::istream* f_in);

void writeArray(double* arr, int in, int out, std::ostream* f_out);

void writeArrayEndl(double* arr, int in, int out, std::ostream* f_out);

void writeArrayWithTabs(double* arr, int in, int out, std::ostream* f_out);

DataGrid readDataGrid(std::istream* f_in);

DataGrid readDataGrid(std::string fileName);
