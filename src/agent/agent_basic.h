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
		Eigen::VectorXd filter_supply;
		Eigen::VectorXd filter_demand;

		// Equivalent bidding price when control reserve is activated
		Eigen::VectorXd balancing_supply;
		Eigen::VectorXd balancing_demand;
	};

	struct results{
		// Clearing results in the IMO
		double cleared_supply;
		double cleared_demand;

		// Confirmed schedule after TSO redispatch
		double confirmed_supply;
		double confirmed_demand;

		// Actual supply / demand profile after control reserve is activated
		double actual_supply;
		double actual_demand;
	};

	struct settlement_process{
		double EOM;
		double redispatch;
		double balancing;
	};

	struct settlement{
		settlement_process volume_supply;
		settlement_process volume_demand;
		// Cost of supplying service
		settlement_process cost;
		// Price of using service
		settlement_process price;
		// Utility of using service
		settlement_process utility;
	};
}
