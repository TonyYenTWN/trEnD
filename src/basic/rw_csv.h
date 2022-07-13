// Header file for read and write csv
#pragma once

#ifndef RW_CSV
#define RW_CSV

#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "src/basic/Basic_Definitions.h"

std::vector <int> get_file_dim(std::string);
Eigen::MatrixXd read_file(int, int, std::string);
void write_file(Eigen::MatrixXd, std::string, std::vector<std::string>);

#endif
