// Header file for operation and investment of industrial consumers
#pragma once
#include "agent.h"

namespace agent{
	namespace industrial{
		static inline double flexible_ratio(){
			double value = .05;
			return value;
		}

		struct profile{
			int point_ID;
			int node_ID;
			Eigen::VectorXd submitted_demand;
		};

		struct profiles{
			std::vector <profile> HV;
			std::vector <profile> LV;
		};
	}
}
