// Header file for agents
#pragma once
#include "src/alglib/optimization.h"
#include "src/basic/basic_definitions.h"
#include "src/basic/eigen_sparse.h"

namespace agent{
	namespace aggregator{
		namespace parameters{
			static inline int foresight_time(){
				int value = 24;
				return value;
			}
		}
	}

	namespace parameters{
		static inline double residential_ratio(){
			//double value = 1.;
			double value = .4;
			return value;
		}
	}

	struct bids{
		// Submitted bids into the IMO
		Eigen::VectorXd submitted_supply_inflex;
		Eigen::VectorXd submitted_demand_inflex;
		Eigen::VectorXd submitted_supply_flex;
		Eigen::VectorXd submitted_demand_flex;

		// Equivalent bidding price when considering redispatch
		Eigen::VectorXd redispatch_supply_inflex;
		Eigen::VectorXd redispatch_demand_inflex;
		Eigen::VectorXd redispatch_supply_flex;
		Eigen::VectorXd redispatch_demand_flex;

		// Filtered bids into redisatch of TSO
		Eigen::VectorXd filtered_supply_inflex;
		Eigen::VectorXd filtered_demand_inflex;
		Eigen::VectorXd filtered_supply_flex;
		Eigen::VectorXd filtered_demand_flex;

		// Confirmed bids after TSO redispatch
		Eigen::VectorXd confirmed_supply_inflex;
		Eigen::VectorXd confirmed_demand_inflex;
		Eigen::VectorXd confirmed_supply_flex;
		Eigen::VectorXd confirmed_demand_flex;

		// Equivalent bidding price when control reserve is activated
		Eigen::VectorXd actual_supply_inflex;
		Eigen::VectorXd actual_demand_inflex;
		Eigen::VectorXd actual_supply_flex;
		Eigen::VectorXd actual_demand_flex;
	};

	struct sorted_vector{
		Eigen::VectorXi id;
		Eigen::VectorXd value;
	};

	// Functions
	sorted_vector sort(Eigen::VectorXd);
}
