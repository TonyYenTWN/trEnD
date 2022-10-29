// Header file for cross-border flow (model as an inflexible agent)
#pragma once
#include "agent_basic.h"

namespace agent{
	namespace cross_border{
		struct profile{
			int node_ID;
			bids bids;
			results results;
			settlement settlement;
		};
		typedef std::vector <profile> profiles;

		struct zonal_profile{
			int exchange_zone_ID;
			int node_num;
			profiles profiles;
		};
		typedef std::vector <zonal_profile> zonal_profiles;
	}
}
