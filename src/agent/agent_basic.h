// Header file for agents
#pragma once
#include "src/alglib/optimization.h"
#include "src/basic/basic_definitions.h"
#include "src/basic/eigen_sparse.h"

namespace agent{
	namespace parameters{
		static inline double residential_ratio(){
		    // default is .4
			//double value = 0.;
			double value = .4;
			return value;
		}
	}

	struct bids_struct{
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

		// Equivalent bids when control reserve is activated
		Eigen::VectorXd imbalance_demand;
		Eigen::VectorXd imbalance_supply;
		Eigen::VectorXd balancing_supply;
		Eigen::VectorXd balancing_demand;
	};

	struct results_struct{
		// Clearing results in the IMO
		double cleared_supply;
		double cleared_demand;

		// Confirmed schedule after TSO redispatch
		double confirmed_supply;
		double confirmed_demand;

		// Imbalance in real time
		double imbalance_supply;
		double imbalance_demand;

		// Actual supply / demand profile after control reserve is activated
		double actual_supply;
		double actual_demand;
	};

	struct settlement_process{
		double EOM;
		double redispatch;
		double imbalance;
		double balancing;
		double BESS;
	};

	struct settlement_struct{
		settlement_process volume_supply;
		settlement_process volume_demand;
		settlement_process volume_supply_up;
		settlement_process volume_supply_down;
		settlement_process volume_demand_up;
		settlement_process volume_demand_down;
		// Cost of supplying service
		settlement_process cost_supply;
		settlement_process cost_demand;
		// Price of using service
		settlement_process price;
		// Reimburse of using service
		settlement_process reimburse;
		// Utility of using service
		settlement_process utility_supply;
		settlement_process utility_demand;
	};
}
