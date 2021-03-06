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
		/**
		* @name scale information
		*/
		/*@{*/
		/** Area size of a point in square kilometers.*/
		double point_area = 100.;
		/** Rough spacing distance between 2 adjacent points in meters.*/
		double grid_length = 10000.;
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
		/** A matrix mapping point indexs onto x-y coordinates. Row indexs represent the x coordinates and column indexs represent the y coordinates.*/
		Eigen::MatrixXi coordinate_grid;
		/** Distance between the points.*/
		Eigen::MatrixXd distance;
		/** Covariance of the reference random fields between the points.*/
		Eigen::SparseMatrix <double> covariance;
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

	/** @brief Transmission network edges (power lines) from the power network data.*/
	struct edges{
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
		/** Length of the power line in meters.*/
		Eigen::VectorXd distance;
		/*@{*/

		/**
		* @name features information
		*/
		/*@{*/
		/** Reference value of line-to-line voltage (in kV) for non-dimensionalization into p.u..*/
		Eigen::VectorXi voltage_base;
		/*@{*/
	};

	/** @brief Information of cross-border transmission with neighboring nations on the boundary.*/
	struct cbt{
		/** Names of the bidding zones included in the internationally-coupled market model.*/
		std::vector <std::string> bz_names;
		/** Nodes where power sources / sinks of bidding zones from neighbor nations are located.*/
		Eigen::MatrixXd entry_nodes;
		/***/
		Eigen::VectorXi entry_node_num;
		/** Constraint of flow exchange between bidding zones internationally-coupled market model, using NTC (net transmission capacity).
		* Rows: bidding zones exporting power; cols: bidding zones importing power.
		* Note the matrix can be asymmetric to capture the conditions of the bidding zones.*/
		Eigen::MatrixXd flow_constraint;
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
		* @name voltage boundaries of the power network
		*/
		/*@{*/
		/** Cutoff voltage level of the transmission network.*/
		int voltage_cutoff_trans = 132;
		/** Voltage level of connection between the transmission and the distribution network.*/
		int voltage_cutoff_connection = 66;
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
		/** Series impedance (ohm per meter) of transmission line.*/
		std::complex<double> z_trans_series = std::complex<double> (0., 5. * pow(10., -4.));
		/** Shunt impedance (ohm per meter) of transmission line.*/
		std::complex<double> z_trans_shunt = std::complex<double> (0., 0.);
		/** Series impedance (ohm per meter) of distribution line.*/
		std::complex<double> z_distr_series = std::complex<double> (0., 7. * pow(10., -4.));
		/** Shunt impedance (ohm per meter) of distribution line.*/
		std::complex<double> z_distr_shunt = std::complex<double> (0., 0.);
		/**Phase angle limits on a node.*/
		double theta_limit = boost::math::constants::pi<double>() / 18.;
		/**Hash table (mapping) of per phase power flow limits on an edge at different voltage base levels, in MW.*/
		std::map <int, double> power_limit;
		/*@{*/

		/**
		* @name non-dimensionalization parameters
		*/
		/*@{*/
		/** Reference value of power on the line (in MW, per phase) for non-dimensionalization into p.u..*/
		double s_base = 1.;
		/** Hash table (mapping) of line to line base voltage levels in kV.*/
		std::map <int, int> voltage_base_levels;
		/** Hash table (mapping) of per phase base impedance levels in Ohm.*/
		std::map <int, double> impedenace_base_levels;
		/*@{*/

		// Set the maps for levels of voltage and impedance base
		void set_level_maps(nodes &nodes){
			int voltage_min = nodes.voltage_base.minCoeff();
			int voltage_max = voltage_min;
			int level_count = 0;
			this->voltage_base_levels.insert(std::make_pair(voltage_max, level_count));
			this->impedenace_base_levels.insert(std::make_pair(voltage_max, (double) voltage_max * voltage_max / this->s_base / 3.));
			this->power_limit.insert(std::make_pair(voltage_max, (double) voltage_max));

			std::vector <int> voltage_base_sorted(nodes.voltage_base.data(), nodes.voltage_base.data() + nodes.voltage_base.size());
			std::sort(voltage_base_sorted.begin(), voltage_base_sorted.end());
			for(int node_iter = 0; node_iter < nodes.voltage_base.size(); ++ node_iter){
				if(voltage_base_sorted[node_iter] > voltage_max){
					voltage_max = voltage_base_sorted[node_iter];
					level_count += 1;
					this->voltage_base_levels.insert(std::make_pair(voltage_max, level_count));
					this->impedenace_base_levels.insert(std::make_pair(voltage_max, (double) voltage_max * voltage_max / this->s_base / 3.));
					this->power_limit.insert(std::make_pair(voltage_max, (double) voltage_max));
					//std::cout << voltage_base_levels[voltage_max] << "\t" <<  voltage_max << "\t" << impedenace_base_levels[level_count] << "\n";
				}
			}

//			// Set power flow limits for different base voltage levels
//			this->power_limit = Eigen::VectorXd(voltage_base_levels.size());
//			for(auto level_iter = ; level_iter < voltage_base_levels.size(); ++ level_iter){
//				this->power_limit(level_iter) = voltage_base_levels[level_iter];
//			}
		}
    };

	struct network_inform{
		points points;
		nodes nodes;
		edges edges;
		cbt cbt;
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
