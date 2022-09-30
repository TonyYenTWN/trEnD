// Header file for agents
#pragma once
#include "src/alglib/optimization.h"
#include "src/basic/basic_definitions.h"
#include "src/basic/eigen_sparse.h"

namespace agent{
	namespace parameters{
		static inline double residential_ratio(){
			//double value = 1.;
			double value = .4;
			return value;
		}
	}

	struct bids{
		// Inflex: bids with dispatch priority
		// Flex: bids that can be redispatched or provide control reserve

		// Submitted bids into the IMO
		Eigen::VectorXd submitted_supply_inflex;
		Eigen::VectorXd submitted_demand_inflex;
		Eigen::VectorXd submitted_supply_flex;
		Eigen::VectorXd submitted_demand_flex;

		// Equivalent bidding price when considering redispatch
		Eigen::VectorXd redispatch_supply;
		Eigen::VectorXd redispatch_demand;

		// Filtered bids into redisatch of TSO (for distributed resources)
		Eigen::VectorXd filtered_supply;
		Eigen::VectorXd filtered_demand;

		// Confirmed bids after TSO redispatch
		Eigen::VectorXd accepted_supply;
		Eigen::VectorXd accepted_demand;

		// Equivalent bidding price when control reserve is activated
		Eigen::VectorXd balancing_supply;
		Eigen::VectorXd balancing_demand;
	};

	struct results{
		// Confirmed schedule after TSO redispatch
		Eigen::VectorXd confirmed_supply;
		Eigen::VectorXd confirmed_demand;

		// Actual supply / demand profile after control reserve is activated
		Eigen::VectorXd actual_supply;
		Eigen::VectorXd actual_demand;
	};

	struct sorted_vector{
		Eigen::VectorXi id;
		Eigen::VectorXd value;
	};

	// Functions
	sorted_vector sort(Eigen::VectorXd);
}
