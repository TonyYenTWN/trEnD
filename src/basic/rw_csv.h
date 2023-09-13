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
	std::vector <int> get_file_dim(std::string, bool row_name = 0, bool col_name = 1);
	std::vector <std::string> get_col_name(std::string, int);
	std::vector <std::string> get_row_name(std::string, int);
	std::map<std::string, std::vector<std::string>> read_config_file(std::string);
	Eigen::MatrixXd read_file(int, int, std::string, bool row_name = 0);
	void write_file(Eigen::MatrixXd, std::string, std::vector<std::string>, int precision = 6);
}

#endif
