// Source file for printing results from simulation and settlement of the power market
#include <iostream>
#include "src/basic/rw_csv.h"
#include "power_market.h"

// ------------------------------------------------------------------------------------------------
// Generic functions for all market types
// ------------------------------------------------------------------------------------------------
namespace{
	void Market_results_print(power_market::market_inform Market, std::string name, bool nodal = 1){
		std::string fout_name;

		// Confirmed results
		fout_name = "csv/output/" + name + "_confirmed_price.csv";
		basic::write_file(Market.confirmed.price, fout_name, Market.zone_names);
		fout_name = "csv/output/" + name + "_confirmed_demand.csv";
		basic::write_file(Market.confirmed.demand, fout_name, Market.zone_names);
		fout_name = "csv/output/" + name + "_confirmed_supply.csv";
		basic::write_file(Market.confirmed.supply, fout_name, Market.zone_names);

		if(!nodal){
			fout_name = "csv/output/" + name + "_confirmed_power.csv";
			basic::write_file(Market.network.confirmed_power, fout_name, Market.zone_names);
			return;
		}

		// EOM settlement
		fout_name = "csv/output/" + name + "_EOM_cost.csv";
		basic::write_file(Market.EOM.cost, fout_name, Market.zone_names);
		fout_name = "csv/output/" + name + "_EOM_utility.csv";
		basic::write_file(Market.EOM.utility, fout_name, Market.zone_names);

		// Redispatch settlement
		fout_name = "csv/output/" + name + "_redispatch_cost.csv";
		basic::write_file(Market.redispatch.cost, fout_name, Market.zone_names);
		fout_name = "csv/output/" + name + "_redispatch_utility.csv";
		basic::write_file(Market.redispatch.utility, fout_name, Market.zone_names);
		fout_name = "csv/output/" + name + "_redispatch_demand.csv";
		basic::write_file(Market.redispatch.demand_up - Market.redispatch.demand_down, fout_name, Market.zone_names);
		fout_name = "csv/output/" + name + "_redispatch_supply.csv";
		basic::write_file(Market.redispatch.supply_up - Market.redispatch.supply_down, fout_name, Market.zone_names);

		// Actual results
		fout_name = "csv/output/" + name + "_actual_price.csv";
		basic::write_file(Market.actual.price, fout_name, Market.zone_names);
		fout_name = "csv/output/" + name + "_actual_demand.csv";
		basic::write_file(Market.actual.demand, fout_name, Market.zone_names);
		fout_name = "csv/output/" + name + "_actual_supply.csv";
		basic::write_file(Market.actual.supply, fout_name, Market.zone_names);

		// Balancing settlement
		fout_name = "csv/output/" + name + "_balancing_cost.csv";
		basic::write_file(Market.balancing.cost, fout_name, Market.zone_names);
		fout_name = "csv/output/" + name + "_balancing_utility.csv";
		basic::write_file(Market.balancing.utility, fout_name, Market.zone_names);
		fout_name = "csv/output/" + name + "_balancing_demand.csv";
		basic::write_file(Market.balancing.demand_up - Market.balancing.demand_down, fout_name, Market.zone_names);
		fout_name = "csv/output/" + name + "_balancing_supply.csv";
		basic::write_file(Market.balancing.supply_up - Market.balancing.supply_down, fout_name, Market.zone_names);
	}
}

void power_market::Markets_results_print(market_whole_inform &Power_market_inform){
	Market_results_print(Power_market_inform.International_Market, "IMO", 0);
	Market_results_print(Power_market_inform.TSO_Market, "TSO");
}
