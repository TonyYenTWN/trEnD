// Header file for market object and functions for market clearing
#pragma once

#ifndef MARKET_OBJECT
#define MARKET_OBJECT

#include "src/alglib/optimization.h"
#include "src/basic/basic_definitions.h"
#include "src/basic/eigen_sparse.h"
#include "src/power_network/power_network.h"

namespace power_market{
	// Power market objects
	/**Information of the corresponding power network of the market.*/
	struct network_graph{
		// Input parameters
		/**Number of vertices (nodes) in the power network of the market.*/
		int num_vertice;
		/**Number of edges (lines) in the power network of the market.*/
		int num_edges;
		/**Compact form of the incidence matrix (0th col: start; 1st col: end).*/
		Eigen::MatrixXi incidence_matrix;
		/**(Imaginary part of) admittance of each edge; used in TSO and DSOs markets.*/
		Eigen::VectorXd admittance_vector;
		/**Voltage constraints at each node; used in TSO and DSOs markets.*/
		Eigen::MatrixXd voltage_constraint;
		/**Power flow constraints at each edge (0th col: from start to end; 1st col: from end to start).*/
		Eigen::MatrixXd power_constraint;

		// Output variables
		/**Power flow across each edge; a positive value indicates flowing from start to end.*/
		Eigen::MatrixXd confirmed_power;
		/**Voltage at each vertex.*/
		Eigen::MatrixXd confirmed_voltage;
	};

	/**Information of the power market.*/
	struct market_inform{
		// Input parameters
		/**Number of bidding zones / nodes in the market.*/
		int num_zone;
		/**Indicate the index, after whose corresponding bidding zones are on the boundary and will only act as cross-border flow (for the IMO market only).*/
		int cross_border_zone_start;
		/**Total time intervals of the model.*/
		int time_intervals;
		/**Total price intervals for flexible supply and demand in the model.*/
		int price_intervals = 600;
		/**Name of bidding zones (for the IMO market only).*/
		std::vector<std::string> zone_names;
		/**Range of lowest and highest possible bidding prices.*/
		Eigen::Vector2d price_range_inflex = Eigen::Vector2d(-500., 3000.);
		/**Range of bidding prices for flexible supply and demand in the model.*/
		Eigen::Vector2d price_range_flex = Eigen::Vector2d(-100., 500.);
		/**Bidding prices in the model.*/
		Eigen::VectorXd bidded_price;
		/**Merit order curve of the supply bids (for the IMO market only).*/
		Eigen::MatrixXd merit_order_curve;
		/**Default demand of biddings zones without cross-border flows (for the IMO market only).*/
		Eigen::MatrixXd demand_default;			// Default demand profiles of the bidding zones; in later runs demand bids from Norway should come from lower level markets

		// Process Variables
		/**Supply bid submitted in the bidding zones.*/
		Eigen::MatrixXd submitted_supply;
		/**Demand bid submitted in the bidding zones.*/
		Eigen::MatrixXd submitted_demand;
		/**Member object that stores the linear programming problem (for TSO and DSO markets).*/
		alglib::minlpstate Problem;

		// Output Variables
		/**Confirmed supply quantity of the market.*/
		Eigen::MatrixXd confirmed_supply;
		/**Confirmed demand quantity of the market.*/
		Eigen::MatrixXd confirmed_demand;
		/**Confirmed market clearing price of the market.*/
		Eigen::MatrixXd confirmed_price;
		/**Filtered supply bids from DSOs forwarding to TSO before redispatch.*/
		Eigen::MatrixXd filtered_supply;
		/**Filtered demand bids from DSOs forwarding to TSO before redispatch.*/
		Eigen::MatrixXd filtered_demand;

		// Mixed Substructure
		/**Information of the corresponding power network of the market.*/
		network_graph network;

		// Functions
		/**Set the bidding prices of the power market.*/
		void set_bidded_price(){
			this->bidded_price = Eigen::VectorXd(this->price_intervals + 2);
			this->bidded_price(0) = this->price_range_inflex(0);
			this->bidded_price.array().tail(1) = this->price_range_inflex(1);
			this->bidded_price.segment(1, this->price_intervals) = Eigen::VectorXd::LinSpaced(this->price_intervals, this->price_range_flex(0) + .5 * (this->price_range_flex(1) - this->price_range_flex(0)) / this->price_intervals, this->price_range_flex(1) - .5 * (this->price_range_flex(1) - this->price_range_flex(0)) / this->price_intervals);
		}
	};

	/**A vector of power markets. This type is used in setting the DSO markets.*/
	typedef std::vector <market_inform> markets_inform;
}

#endif

// generic functions for power markets
#ifndef POWER_MARKET
#define POWER_MARKET

namespace power_market{
	void Market_Initialization(market_inform&);
	void Market_clearing_nodal(int, market_inform&, Eigen::VectorXi&, Eigen::MatrixXd&, Eigen::MatrixXd&);
	void Submitted_bid_calculation(int, markets_inform&, market_inform&, market_inform&, power_network::network_inform&, std::string);
	void Flow_Based_Market_LP_Set(market_inform&, alglib::minlpstate &);
	void Flow_Based_Market_Optimization(int, market_inform&, alglib::minlpstate&);
}

#endif

// Functions for the internationally-coupled market operator
#ifndef IMO
#define IMO

namespace power_market{
	void International_Market_Set(market_inform&, int, std::string, std::string);
	void International_Market_Optimization(int, market_inform&, bool print_result = 1);
	void International_Market_Output(market_inform&);
}

#endif

// Functions for the transmission system operator
#ifndef TSO
#define TSO

namespace power_market{
	void TSO_Market_Set(market_inform&, power_network::network_inform&, int);
}

#endif

// Functions for the distribution system operators
#ifndef DSO
#define DSO

namespace power_market{
	void DSO_Markets_Set(markets_inform&, power_network::network_inform&, int);
}

#endif
