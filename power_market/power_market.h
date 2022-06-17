// Header file for market object and functions for market clearing
#pragma once

#ifndef MARKET_OBJECT
#define MARKET_OBJECT

#include "../basic/Basic_Definitions.h"
#include "../basic/Eigen_Sparse.h"
#include "../basic/rw_csv.cpp"

// Power market objects
struct network_graph{
	// Input parameters
	int num_vertice;
	int num_edges;
	Eigen::MatrixXi incidence_matrix;		// Compact Incidence Matrix (0th col: start; 1st col: end)
	Eigen::VectorXd admittance_vector;		// (Imaginary part of) admittance of each edge; used in TSO and DSO
	Eigen::MatrixXd voltage_constraint;		// Voltage constraint at each node; used in TSO and DSO
	Eigen::MatrixXd power_constraint;  		// Power flow constraint (0th col: from start to end; 1st col: from end to start)

	// Output variables
	Eigen::MatrixXd confirmed_power;		// Power flow across each edge; positive indicates flowing from start to end
	Eigen::MatrixXd confirmed_voltage;		// Voltage at each vertice
};

struct market_inform{
	// Input parameters
	int num_zone;							// Can be the actual bidding zones, or just a node / spatial element
	int cross_border_zone_start;			// Indicate the index, after whose corresponding bidding zones are on the boundary and only act by cross-border flow
	int time_intervals;
	int price_intervals;
	std::vector<std::string> zone_names;
	Eigen::Vector2d price_range_inflex;
	Eigen::Vector2d price_range_flex;
	Eigen::VectorXd bidded_price;
	Eigen::MatrixXd merit_order_curve;
	Eigen::MatrixXd demand_default;			// Default demand profiles of the bidding zones; in later runs demand bids from Norway should come from lower level markets

	// Process Variables
	Eigen::MatrixXd submitted_supply;		// Supply bid submitted in the bidding zones
	Eigen::MatrixXd submitted_demand;		// Demand bid submitted in the bidding zones

	// Output variables
	Eigen::MatrixXd confirmed_supply;		// Confirmed supply quantity from MO or TSO
	Eigen::MatrixXd confirmed_demand;		// Confirmed demand quantity from MO or TSO
	Eigen::MatrixXd confirmed_price;		// Confirmed market clearing price from MO or TSO
	//Eigen::MatrixXd confirmed_bid_supply; // Confirmed supply bids from DSOs forwarding to wholesale electricity market
	//Eigen::MatrixXd confirmed_bid_demand;	// Confirmed demand bids from DSOs forwarding to wholesale electricity market

	// Mixed Substructure
	network_graph network;
};

#endif

#ifndef POWER_MARKET
#define POWER_MARKET

void Market_Initialization(market_inform&);

#endif

#ifndef IMO
#define IMO

market_inform International_Market_Set(int, std::string, std::string);
void International_Market_Optimization(int, market_inform&, bool print_result = 1);

#endif
