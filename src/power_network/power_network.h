// Header file for information of the power system
#pragma once

#ifndef NETWORK_OBJECT
#define NETWORK_OBJECT

#include "src/basic/basic_definitions.h"
#include "src/basic/eigen_sparse.h"
#include "src/basic/rw_csv.h"

namespace power_network{
    /** @brief Information of spatial points in the model.*/
	struct points{
		/*@{*/
		/** Area size of a point, in square kilometers.*/
		double point_area = 100.;
		/** Rough spacing distance between 2 adjacent points, in meters.*/
		double grid_length = 10000.;
		/*@{*/

		/*@{*/
		/** A matrix mapping point indexs onto x-y coordinates. Row indexs represent the x coordinates and column indexs represent the y coordinates.*/
		Eigen::MatrixXi coordinate_grid;
		/** Distance between the points.*/
		Eigen::MatrixXd distance;
		/** Covariance of the reference random fields between the points.*/
		Eigen::SparseMatrix <double> covariance;
		/*@{*/

		/**
		* @name hierarchy information
		*/
		/*@{*/
		/** Bidding zones where the points are located.*/
		Eigen::VectorXi bidding_zone;
		/** Transmission nodes where the power source / sink on the point will be injected to the transmission network.*/
		Eigen::VectorXi node;
		/** Index number of the point in the DSO cluster which it is located.*/
		Eigen::VectorXi in_cluster_ID;
		/*@{*/

		/**
		* @name features information
		*/
		/*@{*/
		/** Population density at the spatial point.*/
		Eigen::VectorXd population_density;
		/*@{*/

		/**
		* @name geometric information
		*/
		/*@{*/
		/** X-coordinate of the spatial point in European grid (ETRS89).*/
		Eigen::VectorXd x;
		/** Y-coordinate of the spatial point in European grid (ETRS89).*/
		Eigen::VectorXd y;
		/** Longitude of the spatial point.*/
		Eigen::VectorXd lon;
		/** Latitude of the spatial point.*/
		Eigen::VectorXd lat;
		/*@{*/
	};

	/** @brief Information of nodes on the transmission network.*/
	struct nodes{
		/**
		* @name hierarchy information
		*/
		/*@{*/
		/** Bidding zones where the nodes are located.*/
		Eigen::VectorXi bidding_zone;
		/** DSO clusters where the nodes are located.*/
		Eigen::VectorXi cluster;
		/** Index number of the node in the DSO cluster which it is located.*/
		Eigen::VectorXi in_cluster_ID;
		/*@{*/

		/**
		* @name features information
		*/
		/*@{*/
		/** Reference value for non-dimensionalization of voltage into p.u..*/
		Eigen::VectorXi voltage_base;
		/*@{*/

		/**
		* @name geometric information
		*/
		/*@{*/
		/** X-coordinate of the node in European grid (ETRS89).*/
		Eigen::VectorXd x;
		/** Y-coordinate of the node in European grid (ETRS89).*/
		Eigen::VectorXd y;
		/** Longitude of the node.*/
		Eigen::VectorXd lon;
		/** Latitude of the node.*/
		Eigen::VectorXd lat;
		/*@{*/
	};

	/** @brief Transmission network edges (power lines) from the original power network data.*/
	struct edges_orig{
		/**
		* @name topological information
		*/
		/*@{*/
		/** Starting node of the power line.*/
		Eigen::VectorXi from;
		/** Ending node of the power line.*/
		Eigen::VectorXi to;
		/** Ending node of the power line.*/
		/*@{*/

		/**
		* @name geometric information
		*/
		/*@{*/
		Eigen::VectorXd distance;
		/*@{*/

		/**
		* @name features information
		*/
		/*@{*/
		/** Reference value for non-dimensionalization of voltage into p.u..*/
		Eigen::VectorXi voltage_base;
		/*@{*/
	};

	/** @brief Transmission network edges (power lines) from simplified power network data.*/
	struct edges_simp{
		/**
		* @name topological information
		*/
		/*@{*/
		/** Starting node of the power line.*/
		Eigen::VectorXi from;
		/** Ending node of the power line.*/
		Eigen::VectorXi to;
		/*@{*/

		/**
		* @name features information
		*/
		/*@{*/
		/** Conductance of the edges in p.u..*/
		Eigen::VectorXd conductance;
		/*@{*/
	};

	/** @brief Power plants of a type of technology.*/
	struct plants_per_tech{
		Eigen::VectorXi node;

		/**
		* @name features information
		*/
		/*@{*/
		/** Type of the power plant.*/
		Eigen::VectorXi type;
		/** Capacity of the power plant.*/
		Eigen::VectorXd cap;
		/*@{*/

		/**
		* @name coordinates information
		*/
		/*@{*/
		/** X-coordinate of the power plant in European grid (ETRS89).*/
		Eigen::VectorXd x;
		/** Y-coordinate of the power plant in European grid (ETRS89).*/
		Eigen::VectorXd y;
		/** Longitude of the power plant.*/
		Eigen::VectorXd lon;
		/** Latitude of the power plant.*/
		Eigen::VectorXd lat;
		/*@{*/
	};

	/** @brief Power plants of all types of technologies.*/
	struct plants_all{
		/** Information of hydroelectric power plants.*/
		plants_per_tech hydro;
		/** Information of wind power plants.*/
		plants_per_tech wind;
	};

	/** @brief Clustered DSOs in the model.*/
	struct DSO_cluster{
		/** Indexes of spatial points included in the DSO cluster.*/
		std::vector <int> points_ID;
		/** Indexes of transmission nodes included in the DSO cluster.*/
		std::vector <int> nodes_ID;
	};

    /** @brief Technical parameters of the power network.
    *
    * Including:
    * - cutoff voltage level of the transmission and distribution network
    * - number and density of power lines of the distribution network
    * - fractional dimension of the power lines of the distribution network
    * - series and shunt impedance of the power lines
    */
	struct technical_parameters{
		/**
		* @name voltage on the power network
		*/
		/*@{*/
		/** Cutoff voltage level of the transmission network.*/
		int voltage_cutoff_trans = 132;
		/** Cutoff voltage level of the distribution network.*/
		int voltage_cutoff_distr = 20;
		/*@{*/

		/**
		* @name statistical parameters of power network
		*/
		/*@{*/
		/** Total number of power lines in the distribution network.*/
		int line_num_distr = 124245;
		/** Density of power lines per point in the distribution network.*/
		double line_density_distr;
		/** Fractional dimension of the distribution network.*/
		double fraction_dim_distr = 1.5;
		/*@{*/

		/**
		* @name physical parameters of power lines
		*/
		/*@{*/
		/** Series impedance per meter of transmission line.*/
		std::complex<double> z_trans_series = std::complex<double> (0., 5. * pow(10., -4.));
		/** Shunt impedance per meter of transmission line.*/
		std::complex<double> z_trans_shunt = std::complex<double> (0., 0.);
		/** Series impedance per meter of distribution line.*/
		std::complex<double> z_distr_series = std::complex<double> (0., 7. * pow(10., -4.));
		/** Shunt impedance per meter of distribution line.*/
		std::complex<double> z_distr_shunt = std::complex<double> (0., 0.);
		/*@{*/

		/**
		* @name non-dimensionalization parameters
		*/
		/*@{*/
		/** Reference value for non-dimensionalization of power into p.u..*/
		std::complex<double> s_base = std::complex<double> (1000., 0.) * pow(3., .5);
		/*@{*/
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
			this->tech_parameters.line_density_distr = (double) this->tech_parameters.line_num_distr / (double) this->points.bidding_zone.size();
		}
	};

	// Function for constructing distance and covariance matrix of points
	void point_distance_cov(points&, double);

	// Function for reading the files
	void power_network_input_process(network_inform&, std::string parent_dir);
}

#endif
