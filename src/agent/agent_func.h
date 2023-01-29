// Header file for generic functions of agents
#pragma once
#include "src/power_market/power_market.h"
#include "aggregator.h"
#include "cross_border.h"
#include "end-user.h"
#include "industrial.h"
#include "power_supplier.h"

namespace agent{
	void agents_set(int, power_market::market_whole_inform&, power_network::network_inform&, std::string);
	void agents_redispatch_update(int, power_market::market_whole_inform&, power_network::network_inform&);
	void agents_filter_demand_update(int, power_market::market_whole_inform&, power_network::network_inform&);
	void agents_filter_supply_update(int, power_market::market_whole_inform&, power_network::network_inform&);
	void agents_balancing_update(int, power_market::market_whole_inform&, power_network::network_inform&);
	void agents_status_update(int, power_market::market_whole_inform&, power_network::network_inform&, bool);
	void agents_submit_update(int, power_market::market_whole_inform&, power_network::network_inform&);
	void agents_results_print(power_market::market_whole_inform&, power_network::network_inform&);
}
