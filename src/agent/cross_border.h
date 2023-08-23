// Header file for cross-border flow (model as an inflexible agent)
#pragma once
#include "agent_basic.h"

namespace agent{
	namespace cross_border{
		struct profile{
			int node_ID;
			bids_struct bids;
			results_struct results;
			settlement_struct settlement;
		};
		typedef std::vector <profile> profiles_typedef;

		struct edge_profile{
			int node_num;
			int entry_bz_ID;
            profiles_typedef profiles;
		};
		typedef std::vector <edge_profile> edge_profiles;
	}
}
