// Source file for the main procedure of the power market clearing
#include "src/agent/agent_func.h"
#include "src/power_market/power_market.h"
#include "src/spatial_field/spatial_field.h"

void power_market::default_demand_set(power_network::network_inform &Power_network_inform, market_whole_inform &Power_market_inform){
	int Time = parameters::Time();

	// Initialization of the IMO
	fin_market fin_market;
	fin_market.dir = "csv/input/power_market/";
	fin_market.moc = fin_market.dir + "merit_order_curve_q_assimilated_2021.csv";
	fin_market.demand = fin_market.dir + "generation_total_forecast_2021.csv";
	fin_market.cbt = fin_market.dir + "cbt_forecast_2021.csv";
	fin_market.solar= fin_market.dir + "generation_solar_forecast_2021.csv";
	fin_market.wind_on = fin_market.dir + "generation_wind_onshore_forecast_2021.csv";
	fin_market.wind_off = fin_market.dir + "generation_wind_offshore_forecast_2021.csv";
	International_Market_Set(Power_market_inform.International_Market, Power_market_inform.IMO_Problem, Power_network_inform, Time, fin_market);
}

void power_market::power_market_process_set(power_network::network_inform &Power_network_inform, market_whole_inform &Power_market_inform, bool DSO_filter_flag, bool control_reserve_flag){
	int Time = parameters::Time();

	// Initialization of processed spatial fields
	spatial_field::fin_field fin_field_processed;
	fin_field_processed.dir = "csv/processed/spatial_field/";
	fin_field_processed.demand = fin_field_processed.dir + "nominal_mean_demand_field_10km_ts_";
	fin_field_processed.imbalance = fin_field_processed.dir + "imbalance_field_10km_ts_";
	fin_field_processed.solar = fin_field_processed.dir + "solar_radiation_field_10km_ts_";
	fin_field_processed.wind_on = fin_field_processed.dir + "wind_onshore_cf_field_10km_ts_";
	spatial_field::spatial_field_store(Power_network_inform, fin_field_processed, Time);

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
	International_Market_Price_Estimation(0, Power_market_inform.International_Market, Power_market_inform.IMO_Problem, Power_network_inform);

	// Bidding strategies of agents
	agent::agents_set(Power_market_inform, Power_network_inform);

	// Ideal market clearing in IMO
	Submitted_bid_calculation(Power_market_inform, Power_network_inform);
	International_Market_Optimization(0, Power_market_inform.International_Market, Power_market_inform.IMO_Problem);

	// Equivalent redispatch bids of agents
	agent::agents_redispatch_update(0, Power_market_inform, Power_network_inform);

	// Redispatch in DSO
	if(DSO_filter_flag){
		power_market::Filtered_bid_demand_calculation(0, Power_market_inform, Power_network_inform);
		agent::agents_filter_demand_update(0, Power_market_inform, Power_network_inform);
		power_market::Filtered_bid_supply_calculation(0, Power_market_inform, Power_network_inform);
		agent::agents_filter_supply_update(0, Power_market_inform, Power_network_inform);
	}

	// Redispatch in TSO
	Confirmed_bid_calculation(0, Power_market_inform, Power_network_inform, DSO_filter_flag);
	Flow_Based_Market_Optimization(Power_market_inform.TSO_Market, Power_market_inform.TSO_Problem);
	TSO_Market_Scheduled_Results_Get(0, Power_market_inform.TSO_Market, Power_market_inform.TSO_Problem);

	// Equivalent balancing bids of agents
	agent::agents_balancing_update(0, Power_market_inform, Power_network_inform, DSO_filter_flag);

	// Control reserve activation in TSO
	if(control_reserve_flag){
		Balancing_bid_calculation(0, Power_market_inform, Power_network_inform);
		Flow_Based_Market_Optimization(Power_market_inform.TSO_Market, Power_market_inform.TSO_Problem);
		TSO_Market_Actual_Results_Get(0, Power_market_inform.TSO_Market, Power_market_inform.TSO_Problem);
	}

	// Update state variables of agents
	agent::agents_status_update(0, Power_market_inform, Power_network_inform, control_reserve_flag);
}

void power_market::power_market_process_update(power_network::network_inform &Power_network_inform, market_whole_inform &Power_market_inform, bool DSO_filter_flag, bool control_reserve_flag){
	int Time = parameters::Time();
	int tick = 1;

	International_Market_Price_Estimation(tick, Power_market_inform.International_Market, Power_market_inform.IMO_Problem, Power_network_inform);

//	International_Market_Price_Estimation(tick, Power_market_inform.International_Market, Power_market_inform.IMO_Problem, Power_network_inform);
//	DSO_agents_update(tick, Power_market_inform.end_user_profiles, Power_market_inform.TSO_Market, Power_market_inform.International_Market, Power_network_inform);
//	Submitted_bid_calculation(tick, Power_market_inform, Power_network_inform, DSO_filter_flag);
//	International_Market_Optimization(tick, Power_market_inform.International_Market, Power_market_inform.IMO_Problem);
//	TSO_boundary_update(tick, Power_market_inform.TSO_Market, Power_market_inform.International_Market, Power_network_inform);
//
//	if(DSO_filter_flag){
//		Filtered_bid_calculation(tick, Power_market_inform.DSO_Markets, Power_market_inform.TSO_Market, Power_network_inform, Power_market_inform.DSO_Problems);
//	}
//
//	Flow_Based_Market_Optimization(Power_market_inform.TSO_Market, Power_market_inform.TSO_Problem);
//	TSO_Market_Results_Get(tick, Power_market_inform.TSO_Market, Power_market_inform.TSO_Problem);
//	if(control_reserve_flag){
//		TSO_Market_control_reserve(tick, Power_market_inform, Power_network_inform, DSO_filter_flag);
//	}
}
