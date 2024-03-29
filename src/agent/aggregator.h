// Header file for aggregators
#pragma once
#include "agent_basic.h"

namespace agent{
	namespace aggregator{
		namespace parameters{
			static inline int foresight_time(){
				int value = 24;
				return value;
			}

			static inline double arbitrage_demand(){
				double value = 0.; // Temporary change to non-zero value for sensitivity test
				return value;
			}

			static inline double arbitrage_supply(){
				double value = 0.; // Temporary change to non-zero value for sensitivity test
				return value;
			}
		}

		struct profile{
			int point_ID;
			double arbitrage_demand;
			double arbitrage_supply;
			Eigen::VectorXd price_expected_profile;
			Eigen::VectorXd price_demand_profile;
			Eigen::VectorXd price_supply_profile;
			bids_struct bids;
			results_struct results;
			settlement_struct settlement;
		};

		typedef std::vector <profile> profiles;
	}
}
