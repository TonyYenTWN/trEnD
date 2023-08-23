// Header file for operation and investment of industrial consumers
#pragma once
#include "agent_basic.h"

namespace agent{
	namespace industrial{
		namespace parameters{
			static inline double flexible_ratio(){
				double value = .1;
				return value;
			}

			static inline double power_factor(){
				double value = 1.;
				return value;
			}
		}

		struct profile{
			int point_ID;
			bids_struct bids;
			results_struct results;
			settlement_struct settlement;
		};

		struct profiles{
			std::vector <profile> HV;
			std::vector <profile> LV;
		};
	}
}
