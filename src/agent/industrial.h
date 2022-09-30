// Header file for operation and investment of industrial consumers
#pragma once
#include "agent_basic.h"

namespace agent{
	namespace industrial{
		static inline double flexible_ratio(){
			double value = .05;
			return value;
		}

		struct profile{
			int point_ID;
			bids bids;
		};

		struct profiles{
			std::vector <profile> HV;
			std::vector <profile> LV;
		};
	}
}
