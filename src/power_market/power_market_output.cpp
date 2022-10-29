// Source file for printing results from simulation and settlement of the power market
#include <iostream>
#include "src/basic/rw_csv.h"
#include "power_market.h"

// ------------------------------------------------------------------------------------------------
// Generic functions for all market types
// ------------------------------------------------------------------------------------------------
namespace{
	void Market_results_print(power_market::market_inform Market, std::string name){
		std::string fout_name = "csv/output/" + name + "_confirmed_price.csv";
		basic::write_file(Market.confirmed.price, fout_name, Market.zone_names);
	}
}

void power_market::Markets_results_print(market_whole_inform &Power_market_inform){
	Market_results_print(Power_market_inform.International_Market, "IMO");
	Market_results_print(Power_market_inform.TSO_Market, "TSO");
}

