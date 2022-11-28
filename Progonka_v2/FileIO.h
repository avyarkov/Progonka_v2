#pragma once
#include <iostream>
#include <string>
#include "OldDataGrid.h"

void readArray(double* arr, int in, int out, std::istream* f_in);

void writeArray(double* arr, int in, int out, std::ostream* f_out);

void writeArrayEndl(double* arr, int in, int out, std::ostream* f_out);

void writeArrayWithTabs(double* arr, int in, int out, std::ostream* f_out);

OldDataGrid readOldDataGrid(std::istream* f_in);
OldDataGrid readOldDataGrid(std::string fileName);

OldDataGridNodeSources readOldDataGridNodeSources(std::istream* f_in);
OldDataGridNodeSources readOldDataGridNodeSources(std::string fileName);
