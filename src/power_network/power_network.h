// Header file for information of the power system
#pragma once

#ifndef NETWORK_OBJECT
#define NETWORK_OBJECT

#include "src/basic/basic_definitions.h"
#include "src/basic/eigen_sparse.h"
#include "src/basic/rw_csv.h"
#include <iostream>

namespace power_network{
	namespace parameters{
		static inline double loss_factor(){
			double value = 0.;
			return value;
		}
	}

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
		/** Time series of mean electricity demand (kWh / person / km^2) field at each spatial points.*/
		Eigen::MatrixXd nominal_mean_demand_field;
		/** Time series of imbalance (MWh / MWh) field at each spatial points.*/
		Eigen::MatrixXd imbalance_field;
		/** Time series of onshore wind capacity factor field at each spatial points.*/
		Eigen::MatrixXd wind_on_cf;
		/** Time series of solar capacity factor field at each spatial points.*/
		Eigen::MatrixXd solar_cf;
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
		/** Number of nodes for an edge.*/
		Eigen::VectorXi entry_node_num;
		/** Corresponding bidding zone of an edge.*/
		Eigen::VectorXi entry_bz;
		/** Constraint of flow exchange between bidding zones internationally-coupled market model, using NTC (net transmission capacity).
		* Rows: bidding zones exporting power; cols: bidding zones importing power.
		* Note the matrix can be asymmetric to capture the conditions of the bidding zones.*/
		Eigen::MatrixXd flow_constraint;
	};

	/** @brief Information of weather stations.*/
	struct weather_stations{
		/** Names of the weather stations.*/
		std::vector <std::string> station_names;
		/** Points of the weather stations.*/
		Eigen::VectorXi point;

		/** X-coordinate of the node in European grid (ETRS89).*/
		Eigen::VectorXd x;
		/** Y-coordinate of the node in European grid (ETRS89).*/
		Eigen::VectorXd y;
		/** Longitude of the node.*/
		Eigen::VectorXd lon;
		/** Latitude of the node.*/
		Eigen::VectorXd lat;
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
		* @name voltage and power boundaries of the power network
		*/
		/*@{*/
		/** Cutoff voltage level of the transmission network.*/
		int voltage_cutoff_trans = 132;
		/** Voltage level of connection between the transmission and the distribution network.*/
		int voltage_cutoff_conn = 66;
		/** Cutoff voltage level of the distribution network.*/
		int voltage_cutoff_distr = 22;
		/** Power carrying capacity per voltage (MW / kV) of transmission lines.*/
		double power_limit_trans = 1.;
		/** Power carrying capacity per voltage (MW / kV) of HV distribution lines.*/
		double power_limit_conn = 1.;
		/** Power carrying capacity per voltage (MW / kV) of LV distribution lines.*/
		double power_limit_distr = .5;
		/*@{*/

		/**
		* @name statistical parameters of power network
		*/
		/*@{*/
		/** Density of power lines per point connecting distribution and transmission power network.*/
		double line_density_conn = 1.;
		/** Total number of power lines in the distribution network.*/
		int line_num_distr = 122782; // >=0: 126435; >= 40: 3653
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
		std::complex<double> z_trans_series = std::complex<double> (3. * pow(10., -5.), 2.5 * pow(10., -4.));
		/** Shunt admittance (ohm^(-1) per meter) of transmission line.*/
		std::complex<double> y_trans_shunt = std::complex<double> (0., 2 * pow(10., -9.));
		/** Series impedance (ohm per meter) of HV distribution line.*/
		std::complex<double> z_conn_series = std::complex<double> (3. * pow(10., -4.), 4. * pow(10., -4.));
		/** Shunt admittance (ohm^(-1) per meter) of HV distribution line.*/
		std::complex<double> y_conn_shunt = std::complex<double> (0., 2.5 * pow(10., -9.));
		/** Series impedance (ohm per meter) of MV distribution line.*/
		std::complex<double> z_distr_series = std::complex<double> (6. * pow(10., -4.), 6. * pow(10., -4.));
		/** Shunt admittance (ohm^(-1) per meter) of MV distribution line.*/
		std::complex<double> y_distr_shunt = std::complex<double> (0., 3 * pow(10., -9.));
		/**Phase angle limits on a transmission node.*/
		double theta_trans_limit = boost::math::constants::pi<double>() / 9.;
		/**Phase angle limits on a distribution node.*/
		double theta_distr_limit = boost::math::constants::pi<double>() / 18.;
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
		void set_level_maps(edges &edges){
			int voltage_min = edges.voltage_base.minCoeff();
			int voltage_max = voltage_min;
			int level_count = 0;
			this->voltage_base_levels.insert(std::make_pair(voltage_max, level_count));
			this->impedenace_base_levels.insert(std::make_pair(voltage_max, (double) voltage_max * voltage_max / this->s_base));
			this->power_limit.insert(std::make_pair(voltage_max, (double) this->power_limit_trans * voltage_max));
			//std::cout << voltage_base_levels[voltage_max] << "\t" <<  voltage_max << "\t" << impedenace_base_levels[voltage_max] << "\n";

			std::vector <int> voltage_base_sorted(edges.voltage_base.data(), edges.voltage_base.data() + edges.voltage_base.size());
			std::sort(voltage_base_sorted.begin(), voltage_base_sorted.end());
			for(int edge_iter = 0; edge_iter < edges.voltage_base.size(); ++ edge_iter){
				if(voltage_base_sorted[edge_iter] > voltage_max){
					voltage_max = voltage_base_sorted[edge_iter];
					level_count += 1;
					this->voltage_base_levels.insert(std::make_pair(voltage_max, level_count));
					this->impedenace_base_levels.insert(std::make_pair(voltage_max, (double) voltage_max * voltage_max / this->s_base));
					this->power_limit.insert(std::make_pair(voltage_max, (double) this->power_limit_trans * voltage_max));
					//std::cout << voltage_base_levels[voltage_max] << "\t" <<  voltage_max << "\t" << impedenace_base_levels[voltage_max] << "\n";
				}
			}
		}

		// Set z_base for distribution network
		double z_base_conn(){
			double value = this->voltage_cutoff_conn * this->voltage_cutoff_conn / this->s_base;
			return value;
		}
		double z_base_distr(){
			double value = this->voltage_cutoff_distr * this->voltage_cutoff_distr / this->s_base;
			return value;
		}
    };

     /** @brief Input parameters and results for power flow analysis.
    *
    * Including:
    * - nodal admittance matrix
    * - matrix for time series of voltage at each node
    * - matrix for time series of power source / sink at each node
    * - matrix for time series of power flow at each edge
    */
 	struct power_flow{
		Eigen::SparseMatrix <std::complex <double>> nodal_admittance;
		Eigen::SparseQR <Eigen::SparseMatrix <std::complex <double>>, Eigen::COLAMDOrdering <int>> solver_reg;
		Eigen::SparseQR <Eigen::SparseMatrix <std::complex <double>>, Eigen::COLAMDOrdering <int>> solver_hat;
		std::vector <int> PQ_bus;
		std::vector <int> PU_bus;
		std::vector <int> ref_bus;
		Eigen::MatrixXcd voltage;
		Eigen::MatrixXcd power_node;
		Eigen::MatrixXcd power_edge;
	};

	struct network_inform{
		points points;
		nodes nodes;
		edges edges;
		cbt cbt;
		weather_stations weather_stations;
		plants_all plants;
		std::vector <DSO_cluster> DSO_cluster;
		technical_parameters tech_parameters;
		power_flow power_flow;

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
