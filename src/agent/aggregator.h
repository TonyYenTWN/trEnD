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
		}

		struct profile{
			int point_ID;
			Eigen::VectorXd price_expected_profile;
			Eigen::VectorXd price_demand_profile;
			Eigen::VectorXd price_supply_profile;
		};

		typedef std::vector <profile> profiles;
	}
}
