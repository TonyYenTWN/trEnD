// Source file for re-dispatch and tertiary control reserve market clearing of TSO in Norway
#include <iostream>
#include <iomanip>
#include <chrono>
#include "../basic/rw_csv.h"
#include "../basic/alglib/optimization.h"
#include "../power_network/power_network.h"
#include "../power_market/power_market.h"

void TSO_Market_Set_Test_1(market_inform &TSO_Market, int Time){
	// Input parameters of TSO market
	// A trivial test case with 3 nodes, the first bus being the reference connecting to the remaining 2
	TSO_Market.num_zone = 3;
	TSO_Market.time_intervals = Time;
	TSO_Market.set_bidded_price();

	// Set compact node admittance matrix Y_n
	TSO_Market.network.num_vertice = TSO_Market.num_zone;
	TSO_Market.network.num_edges = 3;
	TSO_Market.network.incidence_matrix = Eigen::MatrixXi(TSO_Market.network.num_edges, 2);
	TSO_Market.network.incidence_matrix.row(0) << 0, 1;
	TSO_Market.network.incidence_matrix.row(1) << 0, 2;
	TSO_Market.network.incidence_matrix.row(2) << 1, 2;
	TSO_Market.network.admittance_vector = Eigen::VectorXd(3);
	TSO_Market.network.admittance_vector << 100., 100., 100.;

	// Set voltage and power constraints at each edge
	TSO_Market.network.voltage_constraint = Eigen::MatrixXd(TSO_Market.network.num_vertice, 2);
	TSO_Market.network.voltage_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, -.2);
	TSO_Market.network.voltage_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, .2);
	TSO_Market.network.power_constraint = Eigen::MatrixXd(TSO_Market.network.num_edges, 2);
	TSO_Market.network.power_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, -50.);
	TSO_Market.network.power_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, 50.);

	// Initialization of output variables
	TSO_Market.confirmed_supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_price = Eigen::MatrixXd(Time, TSO_Market.num_zone);
	TSO_Market.network.confirmed_power = Eigen::MatrixXd(Time, TSO_Market.network.num_edges);

	// For the trivial case only: initialize submitted supply and demand bids at each node
	Market_Initialization(TSO_Market);
	TSO_Market.submitted_supply(0, 0) = 10.;
	TSO_Market.submitted_supply(0, 1) = 20.;
	TSO_Market.submitted_supply(0, 2) = 20.;
	TSO_Market.submitted_supply(1, 0) = 5.;
	TSO_Market.submitted_supply(1, 1) = 15.;
	TSO_Market.submitted_supply(1, 2) = 25.;
	TSO_Market.submitted_demand(TSO_Market.price_intervals + 1, 0) = 25.;
	TSO_Market.submitted_demand(TSO_Market.price_intervals + 1, 1) = 15.;
	TSO_Market.submitted_demand(TSO_Market.price_intervals + 1, 2) = 5.;
	TSO_Market.submitted_demand(TSO_Market.price_intervals, 0) = 10.;
	TSO_Market.submitted_demand(TSO_Market.price_intervals, 1) = 7.5;
	TSO_Market.submitted_demand(TSO_Market.price_intervals, 2) = 5.;
	//S = {15, 35, 45}; D = {35, 22.5, 10}
}

void TSO_Market_Set_Test_2(market_inform &TSO_Market, int Time){
	// Input parameters of TSO market
	// A trivial test case with 20 nodes connected as a radial line
	TSO_Market.num_zone = 5000;
	TSO_Market.time_intervals = Time;
	TSO_Market.set_bidded_price();

	// Set node admittance matrix Y_n
	TSO_Market.network.num_vertice = TSO_Market.num_zone;
	TSO_Market.network.num_edges = TSO_Market.num_zone - 1;
	TSO_Market.network.incidence_matrix = Eigen::MatrixXi(TSO_Market.network.num_edges, 2);
	for(int edge_iter = 0; edge_iter < TSO_Market.network.num_edges; ++ edge_iter){
		TSO_Market.network.incidence_matrix.row(edge_iter) << edge_iter, edge_iter + 1;
	}
	TSO_Market.network.admittance_vector = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, 5000.);

	// Set voltage and power constraints at each edge
	TSO_Market.network.voltage_constraint = Eigen::MatrixXd(TSO_Market.network.num_vertice, 2);
	TSO_Market.network.voltage_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, -.1);
	TSO_Market.network.voltage_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, .1);
	TSO_Market.network.power_constraint = Eigen::MatrixXd(TSO_Market.network.num_edges, 2);
	TSO_Market.network.power_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, -150.);
	TSO_Market.network.power_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, 150.);

	// Initialization of output variables
	TSO_Market.confirmed_supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_price = Eigen::MatrixXd(Time, TSO_Market.num_zone);
	TSO_Market.network.confirmed_power = Eigen::MatrixXd(Time, TSO_Market.network.num_edges);

	// For the trivial case only: initialize submitted supply and demand bids at each node
	Market_Initialization(TSO_Market);
	TSO_Market.submitted_supply.leftCols(TSO_Market.num_zone / 2) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 1.);
	TSO_Market.submitted_demand.rightCols(TSO_Market.num_zone / 2) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 1.);
}

int main(){
	network_inform Power_network_inform;
	power_network_input_process(Power_network_inform, "../power_network/");

	market_inform TSO_Market;
	TSO_Market_Set_Test_2(TSO_Market, 1);
	alglib::minlpstate TSO_Problem;

	auto start = std::chrono::high_resolution_clock::now();
	Flow_Based_Market_LP_Set(TSO_Market, TSO_Problem);
	//Flow_Based_Market_LP_Set_Test(TSO_Market, TSO_Problem);
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast <std::chrono::microseconds> (stop - start);
	std::cout << "Set time: " << duration.count() << " microseconds" << "\n\n";
	start = std::chrono::high_resolution_clock::now();
	Flow_Based_Market_Optimization(0, TSO_Market, TSO_Problem);
	//Flow_Based_Market_Optimization_Test(0, TSO_Market, TSO_Problem);
	stop = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast <std::chrono::microseconds> (stop - start);
	std::cout << "Optimization time: " << duration.count() << " microseconds" << "\n\n";

	// Check issues of memory allocation
	TSO_Market_Set_Test_1(TSO_Market, 1);

	start = std::chrono::high_resolution_clock::now();
	Flow_Based_Market_LP_Set(TSO_Market, TSO_Problem);
	//Flow_Based_Market_LP_Set_Test(TSO_Market, TSO_Problem);
	stop = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast <std::chrono::microseconds> (stop - start);
	std::cout << "Set time: " << duration.count() << " microseconds" << "\n\n";
	start = std::chrono::high_resolution_clock::now();
	Flow_Based_Market_Optimization(0, TSO_Market, TSO_Problem);
	//Flow_Based_Market_Optimization_Test(0, TSO_Market, TSO_Problem);
	stop = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast <std::chrono::microseconds> (stop - start);
	std::cout << "Optimization time: " << duration.count() << " microseconds" << "\n\n";
}
