// Header file for market object and functions for market clearing
#pragma once

#ifndef MARKET_OBJECT
#define MARKET_OBJECT

#include "src/basic/Basic_Definitions.h"
#include "src/basic/Eigen_Sparse.h"
#include "src/basic/alglib/optimization.h"
#include "src/power_network/power_network.h"

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
	int price_intervals = 600;
	std::vector<std::string> zone_names;
	Eigen::Vector2d price_range_inflex = Eigen::Vector2d(-500., 3000.);
	Eigen::Vector2d price_range_flex = Eigen::Vector2d(-100., 500.);
	Eigen::VectorXd bidded_price;
	Eigen::MatrixXd merit_order_curve;
	Eigen::MatrixXd demand_default;			// Default demand profiles of the bidding zones; in later runs demand bids from Norway should come from lower level markets

	// Set bidded price
	void set_bidded_price(){
		this->bidded_price = Eigen::VectorXd(this->price_intervals + 2);
		this->bidded_price(0) = this->price_range_inflex(0);
		this->bidded_price.array().tail(1) = this->price_range_inflex(1);
		this->bidded_price.segment(1, this->price_intervals) = Eigen::VectorXd::LinSpaced(this->price_intervals, this->price_range_flex(0) + .5 * (this->price_range_flex(1) - this->price_range_flex(0)) / this->price_intervals, this->price_range_flex(1) - .5 * (this->price_range_flex(1) - this->price_range_flex(0)) / this->price_intervals);
	}

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

struct DSO_Markets{
	std::vector <market_inform> markets;
	//std::vector <LP_object> problem;
};

#endif

// generic functions for power markets
#ifndef POWER_MARKET
#define POWER_MARKET

void Market_Initialization(market_inform&);
void Market_clearing_nodal(int, market_inform&, Eigen::VectorXi&, Eigen::MatrixXd&, Eigen::MatrixXd&);
void Submitted_bid_calculation(int, DSO_Markets&, market_inform&, market_inform&, network_inform&, std::string);
void Flow_Based_Market_LP_Set(market_inform&, alglib::minlpstate&);
void Flow_Based_Market_Optimization(int, market_inform&, alglib::minlpstate&);

#endif

// Functions for the internationally-coupled market operator
#ifndef IMO
#define IMO

void International_Market_Set(market_inform&, int, std::string, std::string);
void International_Market_Optimization(int, market_inform&, bool print_result = 1);
void International_Market_Output(market_inform&);

#endif

// Functions for the transmission system operator
#ifndef TSO
#define TSO

void TSO_Market_Set(market_inform&, network_inform&, int);

#endif

// Functions for the distribution system operators
#ifndef DSO
#define DSO

void DSO_Markets_Set(DSO_Markets&, network_inform&, int);

#endif
