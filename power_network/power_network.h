// Header file for information of the power system
#pragma once

#ifndef NETWORK_OBJECT
#define NETWORK_OBJECT

#include "../basic/Basic_Definitions.h"
#include "../basic/rw_csv.cpp"

struct points{
	Eigen::VectorXi node;
	Eigen::VectorXd population_density;
	Eigen::VectorXd x;
	Eigen::VectorXd y;
	Eigen::VectorXd lon;
	Eigen::VectorXd lat;
};

struct nodes{
	Eigen::VectorXi DSO;
	Eigen::VectorXi voltage_base;
	Eigen::VectorXd x;
	Eigen::VectorXd y;
	Eigen::VectorXd lon;
	Eigen::VectorXd lat;	
};

struct edges{
	Eigen::VectorXd from;
	Eigen::VectorXd to;
	Eigen::VectorXd type;
	Eigen::VectorXd voltage_base;
	Eigen::VectorXd distance;	
};

struct plants{
	std::string type;
	Eigen::VectorXd x;
	Eigen::VectorXd y;
	Eigen::VectorXd lon;
	Eigen::VectorXd lat;	
};

class network_inform{
	public:
		std::complex<double> x_trans_series(0, 5. * pow(10., -4.));		// Series impedence per meter of transmission line
		std::complex<double> x_trans_shunt(0, 0);						// Shunt impedence per meter of transmission line
		std::complex<double> x_distr_series(0, 7. * pow(10., -4.));		// Series impedence per meter of distribution line
		std::complex<double> x_distr_shunt(0, 0);						// Shunt impedence per meter of distribution line
		Eigen::VectorXi DSO_group;										// Which DSO group a DSO belongs to
		Eigen::VectorXi bidding_zone; 									// Which bidding zone a DSO group belongs to
		
	private:
};

#endif