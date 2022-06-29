// Header file for information of the power system
#pragma once

#ifndef NETWORK_OBJECT
#define NETWORK_OBJECT

#include "../basic/Basic_Definitions.h"
#include "../basic/rw_csv.h"

struct points{
	Eigen::MatrixXi coordinate_grid;
	Eigen::VectorXi bidding_zone;
	Eigen::VectorXi node;
	Eigen::VectorXd population_density;
	Eigen::VectorXd x;
	Eigen::VectorXd y;
	Eigen::VectorXd lon;
	Eigen::VectorXd lat;
};

struct nodes{
	Eigen::VectorXi bidding_zone;
	Eigen::VectorXi cluster;
	Eigen::VectorXi voltage_base;
	Eigen::VectorXd x;
	Eigen::VectorXd y;
	Eigen::VectorXd lon;
	Eigen::VectorXd lat;	
};

struct edges_orig{
	Eigen::VectorXi from;
	Eigen::VectorXi to;
	Eigen::VectorXi voltage_base;
	Eigen::VectorXd distance;	
};

struct edges_simp{
	Eigen::VectorXi from;
	Eigen::VectorXi to;
	Eigen::VectorXd conductance;	// non-dimensionalized into p.u.
};

struct plants_per_tech{
	Eigen::VectorXi node;
	Eigen::VectorXi type;
	Eigen::VectorXd cap;
	Eigen::VectorXd x;
	Eigen::VectorXd y;
	Eigen::VectorXd lon;
	Eigen::VectorXd lat;	
};

struct plants_all{
	plants_per_tech hydro;
	plants_per_tech wind;
};

struct DSO_cluster{
	std::vector <int> points_ID;
	std::vector <int> nodes_ID;
};

struct technical_parameters{
	std::complex<double> x_trans_series = std::complex<double> (0., 5. * pow(10., -4.));	// Series impedence per meter of transmission line
	std::complex<double> x_trans_shunt = std::complex<double>(0., 0.);						// Shunt impedence per meter of transmission line
	std::complex<double> x_distr_series = std::complex<double>(0., 7. * pow(10., -4.));		// Series impedence per meter of distribution line
	std::complex<double> x_distr_shunt = std::complex<double>(0., 0.);						// Shunt impedence per meter of distribution line
	std::complex<double> s_base = std::complex<double>(1000., 0.) * pow(3., .5);			// Reference value for non-dimensionalization of power into p.u.															// Reference values for non-dimensionalization into p.u.
};

struct network_inform{
	points points;
	nodes nodes;
	edges_orig edges_orig;
	edges_simp edges_simp;
	plants_all plants;
	std::vector <DSO_cluster> DSO_cluster;
};

#endif