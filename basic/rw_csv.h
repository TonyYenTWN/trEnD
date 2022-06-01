// Geostat functions header file
#pragma once

#ifndef RW_CSV
#define RW_CSV

#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "../basic/Basic_Definitions.h"

Eigen::MatrixXd read_file(int, int, std::string);
void write_file(Eigen::MatrixXd, std::string, std::vector<std::string>);

#endif