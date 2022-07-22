// Header file for read and write csv
#pragma once

#ifndef RW_CSV
#define RW_CSV

#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "src/basic/basic_definitions.h"

namespace basic{
	std::vector <std::string> get_col_name(std::string, int);
	std::vector <int> get_file_dim(std::string);
	Eigen::MatrixXd read_file(int, int, std::string);
	void write_file(Eigen::MatrixXd, std::string, std::vector<std::string>);
}

#endif
