// Source file for the main procedure of the power market clearing
#include <iostream>
//#include <chrono>
//#include "../basic/LP_gpa.h"
#include "../basic/rw_csv.h"
#include "../basic/alglib/optimization.h"
#include "../power_network/power_network.h"
#include "power_market.h"

int main(){
	// Initialization of power network information
	network_inform Power_network_inform;
	Power_network_inform.set_line_density();
	power_network_input_process(Power_network_inform);
	
	// Initialization of the IMO
	int Time = 8760;
	std::string fin_name_moc = "input/merit_order_curve_q_assimilated_2021.csv";
	std::string fin_name_demand = "input/residual_load_default_forecast_2021.csv";
	market_inform International_Market;
	International_Market_Set(International_Market, Time, fin_name_moc, fin_name_demand);
	
	// Initialization of the TSO
	market_inform TSO_Market;
	alglib::minlpstate TSO_Problem; 
	TSO_Market_Set(TSO_Market, Power_network_inform, 1);
	Flow_Based_Market_LP_Set(TSO_Market, TSO_Problem);		
	
	// Initialization of the DSO
	DSO_Markets DSO_Markets;
	DSO_Markets_Set(DSO_Markets, Power_network_inform, 1);

	// Re-initialization of submitted bids
	std::string fin_point_demand = "../spatial_field/processed/nominal_mean_demand_field_10km_annual_mean.csv";
	Submitted_bid_calculation(0, DSO_Markets, TSO_Market, International_Market, Power_network_inform, fin_point_demand);
	
	// Bid-filtering in DSOs
	
	// Ideal market clearing in IMO
	
	// Re-dispatch in TSO
	Flow_Based_Market_Optimization(0, TSO_Market, TSO_Problem);
}