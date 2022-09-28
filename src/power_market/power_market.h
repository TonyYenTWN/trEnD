// Header file for market object and functions for market clearing
#pragma once

#ifndef MARKET_OBJECT
#define MARKET_OBJECT

#include "src/agent/end-user.h"
#include "src/agent/industrial.h"
#include "src/alglib/optimization.h"
#include "src/basic/basic_definitions.h"
#include "src/basic/eigen_sparse.h"
#include "src/power_network/power_network.h"

namespace power_market{
	namespace parameters{
		static inline int Time(){
			int value = 8760;
			return value;
		}

		static inline int price_interval(){
			int value = 600;
			return value;
		}
	}

	// Power market objects
	/** Information of the corresponding power network of the market.*/
	struct network_graph{
		/**
		* @name input parameters
		*/
		/*@{*/
		/**Number of vertices (nodes) in the power network of the market.*/
		int num_vertice;
		/**Number of edges (lines) in the power network of the market.*/
		int num_edges;
		/**Line capacity matrix for IMO. Term (i, j) represents the maximum power flow capacity from node #i to #j.*/
		Eigen::MatrixXd line_capacity_matrix;
		/**Compact form of the incidence matrix (0th index: start; 1st index: end).*/
		std::vector <Eigen::Vector2i> incidence;
		/**(Imaginary part of) admittance of each edge; used in TSO and DSOs markets.*/
		std::vector <double> admittance;
		/**Power flow constraints at each edge.*/
		std::vector <double> power_limit;

		/**Voltage constraints at each node; used in TSO and DSOs markets.*/
		Eigen::MatrixXd voltage_constraint;
		/**Power flow constraints at each edge.
		* (For IMO) 0th col: from start to end; 1st col: from end to start.
		* (For flow-based markets) 0th col: lower bound; 1st col: upper bound.
		*/
		Eigen::MatrixXd power_constraint;
		/*@{*/

		/**
		* @name output parameters
		*/
		/*@{*/
		/**Power flow across each edge; a positive value indicates flowing from start to end.*/
		Eigen::MatrixXd confirmed_power;
		/**Voltage at each vertex.*/
		Eigen::MatrixXd confirmed_voltage;
		/*@{*/
	};

	/**Information of the control reserve market*/
	struct control_reserve_inform{
		Eigen::MatrixXd activated_positive;
		Eigen::MatrixXd activated_negative;

		Eigen::MatrixXd available_ratio_supply;
		Eigen::MatrixXd available_ratio_demand;
		Eigen::MatrixXd submitted_positive_supply;
		Eigen::MatrixXd submitted_positive_demand;
		Eigen::MatrixXd submitted_negative_supply;
		Eigen::MatrixXd submitted_negative_demand;
	};

	/**Information of a power market.*/
	struct market_inform{
		// Input parameters
		/**Number of bidding zones / nodes in the market.*/
		int num_zone;
		/**Indicate the index, after whose corresponding bidding zones are on the boundary and will only act as cross-border flow (for the IMO market only).*/
		int cross_border_zone_start;
		/**Total time intervals of the model.*/
		int time_intervals;
		/**Total price intervals for flexible supply and demand in the model.*/
		int price_intervals = parameters::price_interval();
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
		/**Filtered supply bids from DSOs forwarding to TSO before redispatch.*/
		Eigen::MatrixXd filtered_supply;
		/**Filtered demand bids from DSOs forwarding to TSO before redispatch.*/
		Eigen::MatrixXd filtered_demand;

		// Output Variables
		/**Highest price of filtered supply in the DSOs.*/
		Eigen::MatrixXd filtered_price_supply;
		/**Lowest price of filtered demand in the DSOs.*/
		Eigen::MatrixXd filtered_price_demand;
		/**Ratio of filtered supply to submitted supply at the highest price in the DSOs.*/
		Eigen::MatrixXd filtered_ratio_supply;
		/**Ratio of filtered demand to submitted demand at the lowest price in the DSOs.*/
		Eigen::MatrixXd filtered_ratio_demand;
		/**Confirmed supply quantity of the market.*/
		Eigen::MatrixXd confirmed_supply;
		/**Confirmed demand quantity of the market.*/
		Eigen::MatrixXd confirmed_demand;
		/**Confirmed market clearing price of the market.*/
		Eigen::MatrixXd confirmed_price;
		/**Ratio of supply or demand confirmed at marginal price:
		* -1: 100% demand; 1: 100% supply.
		*/
		Eigen::MatrixXd confirmed_price_ratio;
		/**Actual supply quantity (after real time control reserve activation) of the market.*/
		Eigen::MatrixXd actual_supply;
		/**Actual demand quantity (after real time control reserve activation) of the market.*/
		Eigen::MatrixXd actual_demand;
		/** Market clearing price after control reserve activation of the market.*/
		Eigen::MatrixXd actual_price;
		/**Ratio of supply or demand confirmed at marginal price:
		* -1: 100% demand; 1: 100% supply.
		*/
		Eigen::MatrixXd actual_price_ratio;

		// Mixed Substructure
		/**Information of the corresponding power network of the market.*/
		network_graph network;
		control_reserve_inform control_reserve;

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

	/**A vector of LP problems of power markets. This type is used in setting the LP for DSO markets.*/
	typedef std::vector <alglib::minlpstate> Problems;

	/**Information of the entire power market landscape.*/
	struct market_whole_inform{
		market_inform International_Market;
		alglib::minlpstate IMO_Problem;
		market_inform TSO_Market;
		alglib::minlpstate TSO_Problem;
		markets_inform DSO_Markets;
		Problems DSO_Problems;
		agent::end_user::profiles end_user_profiles;
		Eigen::MatrixXd industrial_submitted_demand;
	};
}

#endif

// generic functions for power markets
#ifndef POWER_MARKET
#define POWER_MARKET

namespace power_market{
	void Market_Initialization(market_inform&);
	void Submitted_bid_calculation(int, market_whole_inform&, power_network::network_inform&, bool);
	void TSO_boundary_update(int, market_inform&, market_inform&, power_network::network_inform&);
	void Flow_Based_Market_LP_Set(market_inform&, alglib::minlpstate&);
	void Flow_Based_Market_Optimization(market_inform&, alglib::minlpstate&);
	void Filtered_bid_calculation(int, markets_inform&, market_inform&, power_network::network_inform&, std::vector <alglib::minlpstate>&);
	void default_demand_set(power_network::network_inform&, market_whole_inform&);
	void power_market_process_set(power_network::network_inform&, market_whole_inform&, bool, bool);
	void power_market_process_update(power_network::network_inform&, market_whole_inform&, bool, bool);
}

#endif

// Functions for the internationally-coupled market operator
#ifndef IMO
#define IMO

namespace power_market{
	struct fin_market{
		std::string dir;
		std::string moc;
		std::string demand;
		std::string cbt;
		std::string wind_on;
		std::string wind_off;
		std::string solar;
	};

	void International_Market_Set(market_inform&, alglib::minlpstate&, power_network::network_inform&, int, fin_market);
	void International_Market_Optimization(int, market_inform&, alglib::minlpstate&);
	void International_Market_Output(market_inform&);
	void International_Market_Price_Estimation(int, market_inform&, alglib::minlpstate&, power_network::network_inform&);
	std::vector <agent::sorted_vector> International_Market_Price_Sorted(int,  market_inform&);
}

#endif

// Functions for the transmission system operator
#ifndef TSO
#define TSO

namespace power_market{
	void TSO_Market_Set(market_inform&, power_network::network_inform&, int);
	void TSO_Market_Results_Get(int, market_inform&, alglib::minlpstate&);
	void TSO_Market_control_reserve(int, market_whole_inform&, power_network::network_inform&, bool);
}

#endif

// Functions for the distribution system operators
#ifndef DSO
#define DSO

namespace power_market{
	void DSO_Markets_Set(markets_inform&, power_network::network_inform&, int);
	agent::end_user::profiles DSO_agents_set(market_inform&, power_network::network_inform&);
	void DSO_agents_update(int, agent::end_user::profiles&, market_inform&, market_inform&, power_network::network_inform&);
	void Source_Node_Set(market_inform&, power_network::DSO_cluster&);
	void Sink_Node_Set(market_inform&, power_network::DSO_cluster&);
	void DSO_Market_Results_Get(int, market_inform&, alglib::minlpstate&, power_network::DSO_cluster&, bool supply = 1);
}

#endif
