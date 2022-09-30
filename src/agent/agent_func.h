// Header file for generic functions of agents
#pragma once
#include "src/power_market/power_market.h"
#include "aggregator.h"
#include "end-user.h"
#include "industrial.h"

namespace agent{
	void agents_set(power_market::market_whole_inform&, power_network::network_inform&);
}
