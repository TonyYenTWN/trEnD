// Main Source File
#include <iostream>
#include <chrono>
#include <omp.h>
#include <boost/math/distributions/normal.hpp>
#include "../basic/rw_csv.cpp"

boost::math::normal norm_dist(0.0, 1.0);

struct network_graph{
	// Input parameters
	int num_vertices;
	int num_edges;
	
	// Compact Incidence Matrix
	// 0th col: start; 1st col: end
	Eigen::MatrixXi incidence_matrix;
	
	// Power flow constrint
	// 0th col: from start to end; 1st col: from end to start
	Eigen::MatrixXd power_constraint;

	// Process variables
	// Check if an edge is still available (flow constraint not reached)
	Eigen::MatrixXi available_edge;
};

struct market_inform{
	
	// Input parameters
	int num_zone;							// Can be the actual bidding zones, or just a node / spatial element
	int time_intervals;
	int price_intervals;
	std::vector<std::string> zone_names;
	Eigen::Vector2d price_range_inflex;
	Eigen::Vector2d price_range_flex;
	Eigen::VectorXd bidded_price;
	Eigen::MatrixXd merit_order_curve;
	network_graph network;
	
	// Process variables
	Eigen::MatrixXd bidded_supply;
	Eigen::MatrixXd bidded_demand;
	
	// Output results
	Eigen::MatrixXd confirmed_supply;
	Eigen::MatrixXd confirmed_demand;
	Eigen::MatrixXd confirmed_price;
	
};

int main(){
	int Time = 8760;
	int prob_intervals = 100;
	market_inform International_Market;
	
	// Parameters of international market
	International_Market.num_zone = 13;
	International_Market.time_intervals = Time;
	International_Market.price_intervals = 600;
	International_Market.zone_names = {"NO1","NO2","NO3","NO4","NO5","DE-LU","DK1","FI","GB","NL","SE1","SE2","SE3"};
	International_Market.price_range_inflex << -500, 3000;
	International_Market.price_range_flex << -100, 500;
	International_Market.network.num_vertices = International_Market.num_zone;
	International_Market.network.num_edges = 15;
	International_Market.network.incidence_matrix = Eigen::MatrixXi(International_Market.network.num_edges, 2);
	International_Market.network.incidence_matrix.col(0) << 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3;
	International_Market.network.incidence_matrix.col(1) << 1, 2, 4, 12, 4, 5, 6, 8, 9, 3, 4, 11, 7, 10, 11;
	International_Market.network.power_constraint = Eigen::MatrixXd(International_Market.network.num_edges, 2);
	International_Market.network.power_constraint.col(0) << 1900, 100, 500, 2130, 300, 1400, 1680, 720, 720, 300, 500, 600, 0, 650, 200;
	International_Market.network.power_constraint.col(1) << 3400, 350, 3900, 2095, 500, 1400, 1150, 720, 720, 1100, 450, 1000, 0, 600, 250;
	International_Market.network.available_edge = Eigen::MatrixXi::Zero(International_Market.network.num_edges, 2);

	International_Market.bidded_price = Eigen::VectorXd(International_Market.price_intervals + 2);
	International_Market.bidded_price(0) = International_Market.price_range_inflex(0);
	International_Market.bidded_price.array().tail(1) = International_Market.price_range_inflex(1);
	International_Market.bidded_price.array().segment(1, International_Market.price_intervals) = Eigen::VectorXd::LinSpaced(International_Market.price_intervals, International_Market.price_range_flex(0) + .5, International_Market.price_range_flex(1) - .5);
	
	// Quantity density at each price
	// Read infered merit order curve data
	std::string fin_name;
	int num_row = International_Market.price_intervals + 2; 
	int num_col = International_Market.num_zone;
	fin_name = "input/merit_order_curve_q_unified.csv";
	Eigen::MatrixXd merit_order_curve_q = read_file(num_row, num_col, fin_name);
	Eigen::MatrixXd diff_merit_order_curve_q = merit_order_curve_q.bottomRows(num_row - 1) - merit_order_curve_q.topRows(num_row - 1);
	International_Market.merit_order_curve = merit_order_curve_q;
	International_Market.merit_order_curve.bottomRows(num_row - 1) = diff_merit_order_curve_q;
	International_Market.merit_order_curve = International_Market.merit_order_curve.array() * International_Market.merit_order_curve.array().max(0);
	
	std::cout << International_Market.merit_order_curve.topRows(200) << std::endl;
}