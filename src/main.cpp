// Source file for the main procedure of the power market clearing
#include <iostream>
//#include <chrono>
//#include "../basic/LP_gpa.h"
#include "basic/rw_csv.h"
#include "alglib/optimization.h"
#include "power_network/power_network.h"
#include "power_market/power_market.h"

int main(){
	// Initialization of power network information
	power_network::network_inform Power_network_inform;
	power_network::power_network_input_process(Power_network_inform, "csv/input/power_network/");
	Power_network_inform.set_line_density();

	// Initialization of the IMO
	int Time = 8760;
	std::string fin_name_moc = "csv/input/power_market/merit_order_curve_q_assimilated_2021.csv";
	std::string fin_name_demand = "csv/input/power_market/residual_load_default_forecast_2021.csv";
	power_market::market_inform International_Market;
	power_market::International_Market_Set(International_Market, Time, fin_name_moc, fin_name_demand);

	// Initialization of the TSO
	power_market::market_inform TSO_Market;
	power_market::TSO_Market_Set(TSO_Market, Power_network_inform, Time);
	alglib::minlpstate TSO_Problem;
	power_market::Flow_Based_Market_LP_Set(TSO_Market, TSO_Problem);

	// Initialization of the DSO
	power_market::DSO_Markets DSO_Markets;
	power_market::DSO_Markets_Set(DSO_Markets, Power_network_inform, Time);
	alglib::minlpstate DSO_Problems;//[DSO_Markets.markets.size()];
	power_market::Flow_Based_Market_LP_Set(DSO_Markets.markets[6], DSO_Problems);
	// create a minlpstate object in market_inform structure to avoid collision
//	//for(int DSO_iter = 0; DSO_iter < DSO_Markets.markets.size(); ++ DSO_iter){
//	for(int DSO_iter = 0; DSO_iter < 10; ++ DSO_iter){
//		if(DSO_Markets.markets[DSO_iter].num_zone > 1){
//		}
//	}

	// Re-initialization of submitted bids
	std::string fin_point_demand = "csv/processed/spatial_field/nominal_mean_demand_field_10km_annual_mean.csv";
	power_market::Submitted_bid_calculation(0, DSO_Markets, TSO_Market, International_Market, Power_network_inform, fin_point_demand);

	// Ideal market clearing in IMO

	// Bid-filtering in DSOs
	//Flow_Based_Market_Optimization(0, DSO_Markets.markets[6], DSO_Problems);

	// Re-dispatch + tertiary control reserve in TSO
	power_market::Flow_Based_Market_Optimization(0, TSO_Market, TSO_Problem);
}
