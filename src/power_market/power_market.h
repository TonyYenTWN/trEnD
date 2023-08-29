// Header file for market object and functions for market clearing
#pragma once

#ifndef MARKET_OBJECT
#define MARKET_OBJECT

#include "src/agent/aggregator.h"
#include "src/agent/cross_border.h"
#include "src/agent/end-user.h"
#include "src/agent/industrial.h"
#include "src/agent/power_supplier.h"
#include "src/alglib/optimization.h"
#include "src/basic/basic_definitions.h"
#include "src/basic/eigen_sparse.h"
#include "src/configuration/configuration.h"
#include "src/power_network/power_network.h"

namespace power_market{
	namespace parameters{
		struct price_ID_bimap{
			Eigen::VectorXd bidded_price;
			std::map <double, int> price_ID;
		};

		static inline int price_interval(){
			int value = 600;
			return value;
		}

		static inline double redispatch_price_max(){
			double value = 1000.;
			return value;
		}

		static inline double balancing_price_max(){
			double value = 1000.;
			return value;
		}

		void bidded_price(price_ID_bimap&);
	}

	// Power market objects
	/** Information of the state of charge limits (for contingency analysis).*/
	struct soc_range{
        Eigen::MatrixXd soc_min;
        Eigen::MatrixXd soc_max;
        Eigen::VectorXd soc_current;
        Eigen::VectorXd capacity_max;
	};

	/** Information of the flexibility situation at each market node (for contingency analysis).*/
	struct flexibility_status{
	    // Processed Variables for power market simulation
	    // (input variables for contingency analysis)
        Eigen::MatrixXd demand_inflex;         // Currently = default demand - smart appliance of end-users + inflexible industrial demand
        Eigen::MatrixXd demand_flex;            // Currently = flexible industrial demand
        Eigen::MatrixXd demand_shiftable;    // Currently = smart appliance of end-users
        Eigen::MatrixXd supply_inflex;           // Currently = actual supply from active end-users
        Eigen::MatrixXd supply_flex;              // Currently = total submitted supply from power supplier + positive redispatch from slack gas

        // Input variables for contingency analysis but does not exist in power market simulation
        soc_range BESS_soc;
        soc_range EV_soc;

        // Process variables for contingency analysis
        Eigen::MatrixXd unfulfilled_demand;
	};

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
		/**Confirmed power flow across each edge; a positive value indicates flowing from start to end.*/
		Eigen::MatrixXd confirmed_power;
		/**Actual power flow across each edge; a positive value indicates flowing from start to end.*/
		Eigen::MatrixXd actual_power;
		/*@{*/
	};

	/**Information of the results of a market.*/
	struct results{
		/**Supply quantity of the market.*/
		Eigen::MatrixXd supply;
		/**Demand quantity of the market.*/
		Eigen::MatrixXd demand;
		/**Market clearing price of the market.*/
		Eigen::MatrixXd price;
		/**Ratio of supply at marginal price (process variable).*/
		Eigen::VectorXd ratio_supply;
		/**Ratio of demand at marginal price (process variable).*/
		Eigen::VectorXd ratio_demand;
	};

	/**Information of the settlement of a market.*/
	struct settlement{
		/**Upward supply quantity of the process.*/
		Eigen::MatrixXd supply_up;
		/**Downward supply quantity of the process.*/
		Eigen::MatrixXd supply_down;
		/**Upward demand quantity of the process.*/
		Eigen::MatrixXd demand_up;
		/**Downward demand quantity of the process.*/
		Eigen::MatrixXd demand_down;
		/**Supply cost of the process*/
		Eigen::MatrixXd price_supply;
		/**Demand cost of the process*/
		Eigen::MatrixXd price_demand;
		/**Upward flexibility cost of the process*/
		Eigen::MatrixXd price_up;
		/**Downward flexibility cost of the process*/
		Eigen::MatrixXd price_down;
		/**Total cost of the process*/
		Eigen::MatrixXd cost;
		/**Total utility of the process*/
		Eigen::MatrixXd utility;
	};

	struct schedule{
		Eigen::MatrixXd EOM;
		Eigen::MatrixXd redispatch;
		Eigen::MatrixXd balancing;
	};

	/**Information of the settlement of a market.*/
	struct schedules{
		schedule cross_border;
		schedule end_user;
		schedule industrial;
		schedule hydro;
		schedule wind;
		schedule pump_storage;
		schedule slack;
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
		/**Name of bidding zones (for the IMO market and TSO market).*/
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
		/** LP object*/
		alglib::minlpstate Problem;
		/*@{*/

		/**
		* @name output variables
		*/
		/*@{*/
		/**Confirmed scheduled result.*/
		results confirmed;
		/**Actual result.*/
		results actual;
		/**Settlement for EOM.*/
		settlement EOM;
		/**Settlement for redispatch.*/
		settlement redispatch;
		/**Settlement for imbalance.*/
		settlement imbalance;
		/**Settlement for balancing.*/
		settlement balancing;
		/**Actual power flow on the corresponding network.*/
		power_network::power_flow_struct power_flow;
		/**Operational schedule.*/
		schedules operation;
		/**Flexibility status.*/
		flexibility_status flex_stat;
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
	struct agent_profiles_struct{
		agent::aggregator::profiles aggregators;
		agent::cross_border::edge_profiles cross_border;
		agent::end_user::profiles end_users;
		Eigen::MatrixXd end_user_type;
		agent::industrial::profiles industrial;
		agent::power_supplier::profiles power_supplier;
	};

	/**Information of the entire power market landscape.*/
	struct market_whole_inform{
		parameters::price_ID_bimap price_map;
		market_inform International_Market;
		market_inform TSO_Market;
		markets_inform DSO_Markets;
		agent_profiles_struct agent_profiles;
	};
}

#endif

// generic functions for power markets
#ifndef POWER_MARKET
#define POWER_MARKET

namespace power_market{
	void Market_Initialization(market_inform&);
	void Operation_Initialization(market_inform&, int);
	void Flow_Based_Market_LP_Set(market_inform&);
	void Flow_Based_Market_Optimization(market_inform&);
	void default_demand_set(power_network::network_inform&, market_whole_inform&, configuration::process_config&);
	void power_market_process_set(power_network::network_inform&, market_whole_inform&, configuration::process_config&);
	void power_market_process_update(power_network::network_inform&, market_whole_inform&, configuration::process_config&);
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

	void International_Market_Set(market_inform&, power_network::network_inform&, int, fin_market);
	void Submitted_bid_calculation(market_whole_inform&, power_network::network_inform&);
	void International_Market_Optimization(int, market_inform&);
	void International_Market_Price_Estimation(int, market_inform&, power_network::network_inform&, configuration::process_config&);
}

#endif

// Functions for the transmission system operator
#ifndef TSO
#define TSO

namespace power_market{
	void TSO_Market_Set(market_inform&, power_network::network_inform&, int);
	void Confirmed_bid_calculation(int, market_whole_inform&, power_network::network_inform&);
	void TSO_Market_Scheduled_Results_Get(int, market_inform&);
	void Balancing_bid_calculation(int, market_whole_inform&, power_network::network_inform&);
	void TSO_Market_Actual_Results_Get(int, market_inform&);
}

#endif

// Functions for the distribution system operators
#ifndef DSO
#define DSO

namespace power_market{
	void DSO_Markets_Set(markets_inform&, power_network::network_inform&, int);
	void Filtered_bid_demand_calculation(int, market_whole_inform&, power_network::network_inform&);
	void Filtered_bid_supply_calculation(int, market_whole_inform&, power_network::network_inform&);
	void DSO_Market_Results_Get(int, market_inform&, power_network::DSO_cluster_struct&, bool supply = 1);
}

#endif

// Functions for output results
#ifndef MARKET_OUTPUT
#define MARKET_OUTPUT

namespace power_market{
	void Markets_results_print(market_whole_inform&, configuration::process_config&);
	void Simplified_network_print(market_whole_inform&, configuration::process_config&);
}

#endif
