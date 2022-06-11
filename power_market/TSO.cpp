// Source File for re-dispatch and tertiary control reserve market clearing of TSO in Norway
#include <iostream>
//#include <chrono>
#include "../basic/LP_gpa.cpp"
#include "power_market.cpp"

market_inform TSO_Market_Set(int Time){
	market_inform TSO_Market;
	
	// Input parameters of TSO market
	// A trivial test case with 2 nodes, the first bus being the reference 
	TSO_Market.num_zone = 2;
	TSO_Market.time_intervals = Time;
	TSO_Market.price_intervals = 600;
	TSO_Market.price_range_inflex << -500, 3000;
	TSO_Market.price_range_flex << -100, 500;
	TSO_Market.bidded_price = Eigen::VectorXd(TSO_Market.price_intervals + 2);
	TSO_Market.bidded_price(0) = TSO_Market.price_range_inflex(0);
	TSO_Market.bidded_price.array().tail(1) = TSO_Market.price_range_inflex(1);
	TSO_Market.bidded_price.array().segment(1, TSO_Market.price_intervals) = Eigen::VectorXd::LinSpaced(TSO_Market.price_intervals, TSO_Market.price_range_flex(0) + .5 * (TSO_Market.price_range_flex(1) - TSO_Market.price_range_flex(0)) / TSO_Market.price_intervals, TSO_Market.price_range_flex(1) - .5 * (TSO_Market.price_range_flex(1) - TSO_Market.price_range_flex(0)) / TSO_Market.price_intervals);
	
	// Set node admittance matrix Y_n
	TSO_Market.network.num_vertice = TSO_Market.num_zone;
	TSO_Market.network.num_edges = 1;
	TSO_Market.network.incidence_matrix = Eigen::MatrixXi(TSO_Market.network.num_edges, 2);
	TSO_Market.network.incidence_matrix.row(0) << 0, 1;
	TSO_Market.network.admittance_vector = Eigen::VectorXd(1);
	TSO_Market.network.admittance_vector << 1;
	TSO_Market.network.Y_n = Eigen::SparseMatrix <double> (TSO_Market.network.num_vertice, TSO_Market.network.num_vertice);
	Eigen::VectorXd Y_n_diag = Eigen::VectorXd::Zero(TSO_Market.network.num_vertice);
	std::vector<Trip> Y_n_trip;
	Y_n_trip.reserve(2 * TSO_Market.network.num_edges + TSO_Market.network.num_vertice);
	for(int edge_iter = 0; edge_iter < TSO_Market.network.num_edges; ++ edge_iter){
		Y_n_trip.push_back(Trip(TSO_Market.network.incidence_matrix(edge_iter, 0), TSO_Market.network.incidence_matrix(edge_iter, 1), -TSO_Market.network.admittance_vector(edge_iter)));
		Y_n_trip.push_back(Trip(TSO_Market.network.incidence_matrix(edge_iter, 1), TSO_Market.network.incidence_matrix(edge_iter, 0), -TSO_Market.network.admittance_vector(edge_iter)));
		Y_n_diag(TSO_Market.network.incidence_matrix(edge_iter, 0)) += TSO_Market.network.admittance_vector(edge_iter);
		Y_n_diag(TSO_Market.network.incidence_matrix(edge_iter, 1)) += TSO_Market.network.admittance_vector(edge_iter);
	}
	for(int node_iter = 0; node_iter < TSO_Market.network.num_vertice; ++ node_iter){
		Y_n_trip.push_back(Trip(node_iter, node_iter, Y_n_diag(node_iter)));
	}
	TSO_Market.network.Y_n.setFromTriplets(Y_n_trip.begin(), Y_n_trip.end());
	
	// Set voltage and power constraints at each edge
	TSO_Market.network.voltage_constraint = Eigen::MatrixXd(TSO_Market.network.num_vertice, 2);
	TSO_Market.network.voltage_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, -.2);
	TSO_Market.network.voltage_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, .2);
	TSO_Market.network.power_constraint = Eigen::MatrixXd(TSO_Market.network.num_edges, 2);
	TSO_Market.network.power_constraint.row(0) << -.1, .1;

	// Initialization of output variables
	TSO_Market.confirmed_supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_price = Eigen::MatrixXd(Time, TSO_Market.num_zone);
	TSO_Market.network.confirmed_power = Eigen::MatrixXd(Time, TSO_Market.network.num_edges);
	
	// For the trivial case only: initialize submitted supply and demand bids at each node
	TSO_Market.submitted_supply = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	TSO_Market.submitted_demand = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	TSO_Market.submitted_supply(0, 0) = 10;
	TSO_Market.submitted_supply(0, 1) = 20;
	TSO_Market.submitted_demand(TSO_Market.price_intervals + 1, 0) = 20;
	TSO_Market.submitted_demand(TSO_Market.price_intervals + 1, 1) = 10;	
	
	return(TSO_Market);
}

void TSO_LP_Set(market_inform &TSO_Market, LP_object &Problem){
	// Set dimension of the problem
	Problem.Constraints_eq_num = TSO_Market.network.num_edges + TSO_Market.network.num_vertice;
	Problem.Constraints_ie_num = 0;
	Problem.Variables_num = 2 * TSO_Market.network.num_edges + TSO_Market.network.num_vertice;

}

void TSO_Market_Optimization(int tick, market_inform &TSO_Market, bool print_result){
	Eigen::MatrixXd bidded_supply = TSO_Market.submitted_supply;
	Eigen::MatrixXd bidded_demand = TSO_Market.submitted_demand;
	
	// Initial market clearing within each nodes
	Eigen::VectorXi default_price_ID(TSO_Market.num_zone);
	Market_clearing_nodal(tick, TSO_Market, default_price_ID, bidded_supply, bidded_demand);
	std::cout << default_price_ID.transpose();
	
	// Initialization of process variables for the main optimization loop
	Eigen::MatrixXd bidded_supply_export = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	Eigen::MatrixXd bidded_demand_export = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	Eigen::MatrixXd bidded_supply_import = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	Eigen::MatrixXd bidded_demand_import = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	for(int zone_ID = 0; zone_ID < TSO_Market.num_zone; ++ zone_ID){
		// Supply curves
		bidded_supply_export.col(zone_ID).array().tail(TSO_Market.price_intervals + 1 - default_price_ID(zone_ID)) = TSO_Market.submitted_supply.col(zone_ID).array().tail(TSO_Market.price_intervals + 1 - default_price_ID(zone_ID));
		bidded_supply_export(default_price_ID(zone_ID), zone_ID) = bidded_supply(default_price_ID(zone_ID), zone_ID);
		bidded_supply_import.col(zone_ID).array().head(default_price_ID(zone_ID)) = TSO_Market.submitted_supply.col(zone_ID).array().head(default_price_ID(zone_ID));
		bidded_supply_import(default_price_ID(zone_ID), zone_ID) = TSO_Market.submitted_supply(default_price_ID(zone_ID), zone_ID) - bidded_supply(default_price_ID(zone_ID), zone_ID);
		
		// Demand curves
		bidded_demand_export.col(zone_ID).array().tail(TSO_Market.price_intervals + 1 - default_price_ID(zone_ID)) = TSO_Market.submitted_demand.col(zone_ID).array().tail(TSO_Market.price_intervals + 1 - default_price_ID(zone_ID));
		bidded_demand_export(default_price_ID(zone_ID), zone_ID) = TSO_Market.submitted_demand(default_price_ID(zone_ID), zone_ID) - bidded_demand(default_price_ID(zone_ID), zone_ID);
		bidded_demand_import.col(zone_ID).array().head(default_price_ID(zone_ID)) = TSO_Market.submitted_demand.col(zone_ID).array().head(default_price_ID(zone_ID));
		bidded_demand_import(default_price_ID(zone_ID), zone_ID) = bidded_demand(default_price_ID(zone_ID), zone_ID);		
	}
	
	// Set the LP problem
	LP_object TSO_Problem;
	TSO_LP_Set(TSO_Market, TSO_Problem);
}

int main(){
	market_inform TSO_Market = TSO_Market_Set(1);
	TSO_Market_Optimization(0, TSO_Market, 1);
}