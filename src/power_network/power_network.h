// Header file for information of the power system
#pragma once

#ifndef NETWORK_OBJECT
#define NETWORK_OBJECT

#include "src/basic/basic_definitions.h"
#include "src/basic/eigen_sparse.h"
#include "src/basic/rw_csv.h"

namespace power_network{
	// Power network objects
	struct points{
		double point_area = 100.; 		// square km
		double grid_length = 10000.; 	// meters
		Eigen::MatrixXi coordinate_grid;
		Eigen::MatrixXd distance;
		Eigen::SparseMatrix <double> covariance;
		Eigen::VectorXi bidding_zone;
		Eigen::VectorXi node;
		Eigen::VectorXi in_cluster_ID;
		Eigen::VectorXd population_density;
		Eigen::VectorXd x;
		Eigen::VectorXd y;
		Eigen::VectorXd lon;
		Eigen::VectorXd lat;
	};

	struct nodes{
		Eigen::VectorXi bidding_zone;
		Eigen::VectorXi cluster;
		Eigen::VectorXi in_cluster_ID;
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

    /** \brief Technical parameters of the power network.
     */
	struct technical_parameters{
		int voltage_cutoff_trans = 132;
		int voltage_cutoff_distr = 20;
		int line_num_distr = 124245;
		double line_density_distr;
		double fraction_dim_distr = 1.5;

		std::complex<double> z_trans_series = std::complex<double> (0., 5. * pow(10., -4.));	/**< Series impedence per meter of transmission line */
		std::complex<double> z_trans_shunt = std::complex<double> (0., 0.);							/**< Shunt impedence per meter of transmission line */
		std::complex<double> z_distr_series = std::complex<double> (0., 7. * pow(10., -4.));	/**< Series impedence per meter of distribution line */
		std::complex<double> z_distr_shunt = std::complex<double> (0., 0.);							/**< Shunt impedence per meter of distribution line */
		std::complex<double> s_base = std::complex<double> (1000., 0.) * pow(3., .5);			/**< Reference value for non-dimensionalization of power into p.u. */
	};

	struct network_inform{
		points points;
		nodes nodes;
		edges_orig edges_orig;
		edges_simp edges_simp;
		plants_all plants;
		std::vector <DSO_cluster> DSO_cluster;
		technical_parameters tech_parameters;

		// Set line density of distribution networks
		void set_line_density(){
			this->tech_parameters.line_density_distr = (double) this->tech_parameters.line_num_distr / (double) this->points.bidding_zone.size() / this->points.point_area;
		}
	};

	// Function for constructing distance and covariance matrix of points
	void point_distance_cov(points&, double);

	// Function for reading the files
	void power_network_input_process(network_inform&, std::string parent_dir);
}

#endif
