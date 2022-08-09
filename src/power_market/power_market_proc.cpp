// Source file for the main procedure of the power market clearing
#include "src/basic/rw_csv.h"
#include "src/alglib/optimization.h"
#include "src/power_network/power_network.h"
#include "src/power_market/power_market.h"

void power_market::power_market_process_set(power_network::network_inform &Power_network_inform, market_whole_inform &Power_market_inform, bool DSO_filter_flag){
	int Time = parameters::Time();

	// Initialization of power network information
	power_network::power_network_input_process(Power_network_inform, "csv/input/power_network/");
	spatial_field::spatial_field_store(Power_network_inform, "csv/processed/spatial_field/nominal_mean_demand_field_10km_ts_", Time);

	// Initialization of the IMO
	std::string fin_name_moc = "csv/input/power_market/merit_order_curve_q_assimilated_2021.csv";
	std::string fin_name_demand = "csv/input/power_market/residual_load_default_forecast_2021.csv";
	International_Market_Set(Power_market_inform.International_Market, Power_network_inform, Time, fin_name_moc, fin_name_demand);

	// Initialization of the TSO
	TSO_Market_Set(Power_market_inform.TSO_Market, Power_network_inform, Time);
	Flow_Based_Market_LP_Set(Power_market_inform.TSO_Market, Power_market_inform.TSO_Problem);

	// Initialization of the DSO
	DSO_Markets_Set(Power_market_inform.DSO_Markets, Power_network_inform, Time);
	Power_market_inform.DSO_Problems = Problems (Power_market_inform.DSO_Markets.size());
	for(int DSO_iter = 0; DSO_iter < Power_market_inform.DSO_Markets.size(); ++ DSO_iter){
		Flow_Based_Market_LP_Set(Power_market_inform.DSO_Markets[DSO_iter], Power_market_inform.DSO_Problems[DSO_iter]);
	}

	// Initial estimation of market clearing price in the IMO
	International_Market_Price_Estimation(0, Power_market_inform.International_Market, Power_network_inform);

	// Bidding strategies of end-users
	Power_market_inform.end_user_profiles = DSO_agents_set(Power_market_inform.International_Market, Power_network_inform);

	// Initialization of submitted bids in DSOs and TSOs
	Submitted_bid_calculation(Power_market_inform.end_user_profiles, Power_market_inform.DSO_Markets, Power_market_inform.TSO_Market, Power_market_inform.International_Market, Power_network_inform, DSO_filter_flag);

	// Ideal market clearing in IMO
	International_Market_Optimization(0, Power_market_inform.International_Market, 0);

	// Set cross-border transmission as boundary conditions of TSO
	TSO_boundary_update(0, Power_market_inform.TSO_Market, Power_market_inform.International_Market, Power_network_inform);

	// Bid-filtering in DSOs
	if(DSO_filter_flag){
		Filtered_bid_calculation(Power_market_inform.DSO_Markets, Power_market_inform.TSO_Market, Power_network_inform, Power_market_inform.DSO_Problems);
	}

	// Re-dispatch + tertiary control reserve in TSO
	Flow_Based_Market_Optimization(Power_market_inform.TSO_Market, Power_market_inform.TSO_Problem);
	TSO_Market_Results_Get(0, Power_market_inform.TSO_Market, Power_market_inform.TSO_Problem);
	TSO_Market_control_reserve(0, Power_market_inform.TSO_Market, Power_market_inform.TSO_Problem);
}

void power_market::power_market_process_update(power_network::network_inform &Power_network_inform, market_whole_inform &Power_market_inform, bool DSO_filter_flag){
	int Time = parameters::Time();

	int tick = 1;
	International_Market_Price_Estimation(tick, Power_market_inform.International_Market, Power_network_inform);
	DSO_agents_update(tick, Power_market_inform.end_user_profiles, Power_market_inform.TSO_Market, Power_market_inform.International_Market, Power_network_inform);
	Submitted_bid_calculation(Power_market_inform.end_user_profiles, Power_market_inform.DSO_Markets, Power_market_inform.TSO_Market, Power_market_inform.International_Market, Power_network_inform, DSO_filter_flag);
	International_Market_Optimization(tick, Power_market_inform.International_Market, 0);
	TSO_boundary_update(tick, Power_market_inform.TSO_Market, Power_market_inform.International_Market, Power_network_inform);

	if(DSO_filter_flag){
		Filtered_bid_calculation(Power_market_inform.DSO_Markets, Power_market_inform.TSO_Market, Power_network_inform, Power_market_inform.DSO_Problems);
	}

	Flow_Based_Market_Optimization(Power_market_inform.TSO_Market, Power_market_inform.TSO_Problem);
	TSO_Market_Results_Get(tick, Power_market_inform.TSO_Market, Power_market_inform.TSO_Problem);
	TSO_Market_control_reserve(tick, Power_market_inform.TSO_Market, Power_market_inform.TSO_Problem);
}
