// Header file for market object and functions for market clearing
#pragma once

#ifndef MARKET_OBJECT
#define MARKET_OBJECT

#include "src/agent/aggregator.h"
#include "src/agent/end-user.h"
#include "src/agent/industrial.h"
#include "src/agent/power_supplier.h"
#include "src/alglib/optimization.h"
#include "src/basic/basic_definitions.h"
#include "src/basic/eigen_sparse.h"
#include "src/power_network/power_network.h"

namespace power_market{
	namespace parameters{
		struct price_ID_bimap{
			Eigen::VectorXd bidded_price;
			std::map <double, int> price_ID;
		};

		static inline int Time(){
			int value = 8760;
			return value;
		}

		static inline int price_interval(){
			int value = 600;
			return value;
		}

		static inline double redispatch_price_max(){
			double value = 3000.;
			return value;
		}

		void bidded_price(price_ID_bimap&);
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
		/*@{*/
	};

	/**Information of the settlement of a market.*/
	struct settlement{
		/**Upward supply quantity of the market.*/
		Eigen::MatrixXd supply_up;
		/**Downward supply quantity of the market.*/
		Eigen::MatrixXd supply_down;
		/**Upward demand quantity of the market.*/
		Eigen::MatrixXd demand_up;
		/**Downward demand quantity of the market.*/
		Eigen::MatrixXd demand_down;
		/**Supply cost of the market*/
		Eigen::MatrixXd cost_supply;
		/**Demand cost of the market*/
		Eigen::MatrixXd cost_demand;
	};

	/**Information of a power market.*/
	struct market_inform{
		/**
		* @name input parameters
		*/
		/*@{*/
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
		/**Bidding prices in the model.*/
		parameters::price_ID_bimap bidded_price_map;
		/**Merit order curve of the supply bids (for the IMO market only).*/
		Eigen::MatrixXd merit_order_curve;
		/**Default demand of biddings zones without cross-border flows (for the IMO market only).*/
		Eigen::MatrixXd demand_default;			// Default demand profiles of the bidding zones; in later runs demand bids from Norway should come from lower level markets
		/*@{*/

		/**
		* @name process variables
		*/
		/*@{*/
		/**Supply bid submitted in the bidding zones.*/
		Eigen::MatrixXd submitted_supply;
		/**Demand bid submitted in the bidding zones.*/
		Eigen::MatrixXd submitted_demand;
		/**Reference prices for each zone.*/
		Eigen::VectorXd reference_price;
		/**Ratio of supply confirmed at marginal price.*/
		Eigen::VectorXd confirmed_ratio_supply;
		/**Ratio of demand confirmed at marginal price.*/
		Eigen::VectorXd confirmed_ratio_demand;
		/**Ratio of actual supply at marginal price.*/
		Eigen::VectorXd actual_ratio_supply;
		/**Ratio of actual demand at marginal price.*/
		Eigen::VectorXd actual_ratio_demand;
		/*@{*/

		/**
		* @name output variables
		*/
		/*@{*/
		/**Confirmed supply quantity of the market.*/
		Eigen::MatrixXd confirmed_supply;
		/**Confirmed demand quantity of the market.*/
		Eigen::MatrixXd confirmed_demand;
		/**Confirmed market clearing price of the market.*/
		Eigen::MatrixXd confirmed_price;
		/**Settlement for redispatch.*/
		settlement redispatch;
		/**Settlement for balancing.*/
		settlement balancing;
		/**Actual supply quantity (after real time control reserve activation) of the market.*/
		Eigen::MatrixXd actual_supply;
		/**Actual demand quantity (after real time control reserve activation) of the market.*/
		Eigen::MatrixXd actual_demand;
		/** Market clearing price after control reserve activation of the market.*/
		Eigen::MatrixXd actual_price;
		/*@{*/

		/**
		* @name mixed substructure
		*/
		/*@{*/
		/**Information of the corresponding power network of the market.*/
		network_graph network;
		/*@{*/
	};

	/**A vector of power markets. This type is used in setting the DSO markets.*/
	typedef std::vector <market_inform> markets_inform;

	/**A vector of LP problems of power markets. This type is used in setting the LP for DSO markets.*/
	typedef std::vector <alglib::minlpstate> Problems;

	/**Information of agents*/
	struct agent_profiles{
		agent::aggregator::profiles aggregators;
		agent::end_user::profiles end_users;
		agent::industrial::profiles industrial;
		agent::power_supplier::profiles power_supplier;
	};

	/**Information of the entire power market landscape.*/
	struct market_whole_inform{
		parameters::price_ID_bimap price_map;
		market_inform International_Market;
		alglib::minlpstate IMO_Problem;
		market_inform TSO_Market;
		alglib::minlpstate TSO_Problem;
		markets_inform DSO_Markets;
		Problems DSO_Problems;
		agent_profiles agent_profiles;
	};
}

#endif

// generic functions for power markets
#ifndef POWER_MARKET
#define POWER_MARKET

namespace power_market{
	void Market_Initialization(market_inform&);
	void TSO_boundary_update(int, market_inform&, market_inform&, power_network::network_inform&);
	void Flow_Based_Market_LP_Set(market_inform&, alglib::minlpstate&);
	void Flow_Based_Market_Optimization(market_inform&, alglib::minlpstate&);
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
	void Submitted_bid_calculation(market_whole_inform&, power_network::network_inform&);
	void International_Market_Optimization(int, market_inform&, alglib::minlpstate&);
	void International_Market_Output(market_inform&);
	void International_Market_Price_Estimation(int, market_inform&, alglib::minlpstate&, power_network::network_inform&);
}

#endif

// Functions for the transmission system operator
#ifndef TSO
#define TSO

namespace power_market{
	void TSO_Market_Set(market_inform&, power_network::network_inform&, int);
	void Confirmed_bid_calculation(int, market_whole_inform&, power_network::network_inform&);
	void TSO_Market_Scheduled_Results_Get(int, market_inform&, alglib::minlpstate&);
	void Balancing_bid_calculation(int, market_whole_inform&, power_network::network_inform&);
	void TSO_Market_Actual_Results_Get(int, market_inform&, alglib::minlpstate&);
}

#endif

// Functions for the distribution system operators
#ifndef DSO
#define DSO

namespace power_market{
	void DSO_Markets_Set(markets_inform&, power_network::network_inform&, int);
	//void DSO_agents_update(int, agent::end_user::profiles&, market_inform&, market_inform&, power_network::network_inform&);
	void Filtered_bid_demand_calculation(int, market_whole_inform&, power_network::network_inform&);
	void Filtered_bid_supply_calculation(int, market_whole_inform&, power_network::network_inform&);
	void DSO_Market_Results_Get(int, market_inform&, alglib::minlpstate&, power_network::DSO_cluster&, bool supply = 1);
}

#endif
