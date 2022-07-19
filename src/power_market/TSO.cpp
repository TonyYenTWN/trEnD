// Source file for re-dispatch and tertiary control reserve market clearing of TSO in Norway
#include <iostream>
//#include <chrono>
//#include "../basic/LP_gpa.h"
#include "src/basic/rw_csv.h"
#include "src/power_network/power_network.h"
#include "power_market.h"

// May be useful in the future (when calculating actual ac power flow)
//double impedence_conversion(Eigen::MatrixXd pu_dc_inform, double voltage){
//	int v_iter = 0;
//	while(pu_dc_inform(v_iter, 0) != voltage){
//		v_iter += 1;
//	}
//
//	return(pu_dc_inform(v_iter, 1));
//}

void power_market::TSO_Market_Set(market_inform &TSO_Market, power_network::network_inform &Power_network_inform, int Time){
//	double pi = boost::math::constants::pi<double>();

	// Input parameters of TSO market
	TSO_Market.num_zone = Power_network_inform.nodes.bidding_zone.size();
	TSO_Market.time_intervals = Time;
	TSO_Market.set_bidded_price();

	// Set node admittance matrix and line capacity matrix
	TSO_Market.network.num_vertice = TSO_Market.num_zone;
	Eigen::MatrixXd admittance = Eigen::MatrixXd::Zero(TSO_Market.network.num_vertice, TSO_Market.network.num_vertice);
	Eigen::MatrixXd capacity = Eigen::MatrixXd::Zero(TSO_Market.network.num_vertice, TSO_Market.network.num_vertice);
	for(int edge_iter = 0; edge_iter < Power_network_inform.edges.distance.size(); ++ edge_iter){
		int from_ID = Power_network_inform.edges.from(edge_iter );
		int to_ID = Power_network_inform.edges.to(edge_iter);
		int voltage = Power_network_inform.edges.voltage_base(edge_iter);
		double y_ij = 1. / Power_network_inform.edges.distance(edge_iter);
		y_ij /= Power_network_inform.tech_parameters.z_distr_series.imag();
		y_ij *= Power_network_inform.tech_parameters.impedenace_base_levels[voltage];
		admittance(from_ID, to_ID) += y_ij ;
		admittance(to_ID, from_ID) += y_ij ;
		capacity(from_ID, to_ID) += Power_network_inform.tech_parameters.power_limit[voltage];
	}

	// Set compact incidence matrix and edge admittance matrix
	double tol = 1E-6;
	TSO_Market.network.incidence.reserve(TSO_Market.network.num_vertice * TSO_Market.network.num_vertice);
	TSO_Market.network.admittance.reserve(TSO_Market.network.num_vertice * TSO_Market.network.num_vertice);
	std::vector <double> power_limit;
	power_limit.reserve(TSO_Market.network.num_vertice * TSO_Market.network.num_vertice);
	for(int row_iter = 0; row_iter < TSO_Market.network.num_vertice - 1; ++ row_iter){
		for(int col_iter = row_iter + 1; col_iter < TSO_Market.network.num_vertice; ++ col_iter){
			if(abs(admittance(row_iter , col_iter)) > tol){
				TSO_Market.network.incidence.push_back(Eigen::Vector2i(row_iter, col_iter));
				TSO_Market.network.admittance.push_back(admittance(row_iter , col_iter));
				power_limit.push_back(capacity(row_iter , col_iter));
			}
		}
	}
	TSO_Market.network.num_edges = TSO_Market.network.incidence.size();

//	TSO_Market.network.node_admittance_matrix = Eigen::SparseView(admittance, tol);
//	TSO_Market.network.line_capacity_matrix = Eigen::SparseView(capacity, tol);
//	TSO_Market.network.num_edges = TSO_Market.network.line_capacity_matrix.nonZeros() / 2;

	// Set compact incidence matrix and edge admittance matrix
//	TSO_Market.network.num_edges = Power_network_inform.edges_simp.from.size();
//	TSO_Market.network.incidence_matrix = Eigen::MatrixXi(TSO_Market.network.num_edges, 2);
//	TSO_Market.network.incidence_matrix.col(0) = Power_network_inform.edges_simp.from;
//	TSO_Market.network.incidence_matrix.col(1) = Power_network_inform.edges_simp.to;
//	TSO_Market.network.admittance_vector = Power_network_inform.edges_simp.conductance;

	// Set voltage and power constraints at each edge
	TSO_Market.network.voltage_constraint = Eigen::MatrixXd::Ones(TSO_Market.network.num_vertice, 2);
	TSO_Market.network.voltage_constraint.col(0) *= -Power_network_inform.tech_parameters.theta_limit;
	TSO_Market.network.voltage_constraint.col(1) *= Power_network_inform.tech_parameters.theta_limit;
	TSO_Market.network.power_constraint = Eigen::MatrixXd (TSO_Market.network.num_edges, 2);
	TSO_Market.network.power_constraint.col(1) = Eigen::Map <Eigen::VectorXd> (power_limit.data(), power_limit.size());
	TSO_Market.network.power_constraint.col(1) /= Power_network_inform.tech_parameters.s_base;
	TSO_Market.network.power_constraint.col(0) = -TSO_Market.network.power_constraint.col(1);

	// Initialization of process variables
	power_market::Market_Initialization(TSO_Market);

	// Initialization of output variables
	TSO_Market.confirmed_supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_price = Eigen::MatrixXd(Time, TSO_Market.num_zone);
	TSO_Market.network.confirmed_voltage = Eigen::MatrixXd(Time, TSO_Market.network.num_vertice);
	TSO_Market.network.confirmed_power = Eigen::MatrixXd(Time, TSO_Market.network.num_edges);
}
