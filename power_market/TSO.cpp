// Source file for re-dispatch and tertiary control reserve market clearing of TSO in Norway
#include <iostream>
//#include <chrono>
//#include "../basic/LP_gpa.h"
#include "../basic/rw_csv.h"
#include "../power_network/power_network.h"
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

void TSO_Market_Set(market_inform &TSO_Market, network_inform &Power_network_inform, int Time){
	double pi = boost::math::constants::pi<double>();
	
	// Input parameters of TSO market
	TSO_Market.num_zone = Power_network_inform.nodes.bidding_zone.size();
	TSO_Market.time_intervals = Time;
	TSO_Market.set_bidded_price();
	
	// Set compact incidence matrix and edge admittance matrix
	TSO_Market.network.num_vertice = TSO_Market.num_zone;
	TSO_Market.network.num_edges = Power_network_inform.edges_simp.from.size();
	TSO_Market.network.incidence_matrix = Eigen::MatrixXi(TSO_Market.network.num_edges, 2);
	TSO_Market.network.incidence_matrix.col(0) = Power_network_inform.edges_simp.from;
	TSO_Market.network.incidence_matrix.col(1) = Power_network_inform.edges_simp.to;
	TSO_Market.network.admittance_vector = Power_network_inform.edges_simp.conductance;

	// Set voltage and power constraints at each edge
	TSO_Market.network.voltage_constraint = Eigen::MatrixXd(TSO_Market.network.num_vertice, 2);
	TSO_Market.network.voltage_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, -pi / 18);
	TSO_Market.network.voltage_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, pi / 18);
	TSO_Market.network.power_constraint = Eigen::MatrixXd(TSO_Market.network.num_edges, 2);
	TSO_Market.network.power_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, -5.);
	TSO_Market.network.power_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, 5.);

	// Initialization of process variables
	Market_Initialization(TSO_Market);
	
	// Initialization of output variables
	TSO_Market.confirmed_supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_price = Eigen::MatrixXd(Time, TSO_Market.num_zone);
	TSO_Market.network.confirmed_power = Eigen::MatrixXd(Time, TSO_Market.network.num_edges);
}