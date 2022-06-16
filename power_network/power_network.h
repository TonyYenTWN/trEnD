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
	Eigen::VectorXd x;
	Eigen::VectorXd y;
	Eigen::VectorXd lon;
	Eigen::VectorXd lat;	
};

class network_inform{
	public:
		
	private:
};

#endif