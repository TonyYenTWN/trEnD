// Source file for dispatch filtering of DSOs in Norway
#include <iostream>
//#include <chrono>
#include "../basic/LP_gpa.h"
#include "../basic/rw_csv.cpp"
#include "../power_network/power_network_input.cpp"
#include "power_market.h"
#include "power_market_func.cpp"
#include "IMO.cpp"
#include "TSO.cpp"

void DSO_Markets_Set(DSO_Markets &DSO_Markets, network_inform &Power_network_inform, int Time){
	DSO_Markets.markets.clear();
	DSO_Markets.markets = std::vector <market_inform> (Power_network_inform.DSO_cluster.size());
	
	// Initialize markets at each clustered DSO
	for(int DSO_iter = 0; DSO_iter < DSO_Markets.markets.size(); ++ DSO_iter){
		// Input parameters of DSO market
		DSO_Markets.markets[DSO_iter].num_zone = Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() + Power_network_inform.DSO_cluster[DSO_iter].nodes_ID.size();
		DSO_Markets.markets[DSO_iter].time_intervals = Time;
		DSO_Markets.markets[DSO_iter].set_bidded_price();

		// Set compact incidence matrix and edge admittance matrix
		// Use fractional Laplacian here
		DSO_Markets.markets[DSO_iter].network.num_vertice = DSO_Markets.markets[DSO_iter].num_zone;
		DSO_Markets.markets[DSO_iter].network.num_edges = DSO_Markets.markets[DSO_iter].num_zone * (DSO_Markets.markets[DSO_iter].num_zone - 1) / 2;
		
		// Initialization of output variables
		DSO_Markets.markets[DSO_iter].confirmed_supply = Eigen::MatrixXd::Zero(Time, DSO_Markets.markets[DSO_iter].num_zone);
		DSO_Markets.markets[DSO_iter].confirmed_demand = Eigen::MatrixXd::Zero(Time, DSO_Markets.markets[DSO_iter].num_zone);
		DSO_Markets.markets[DSO_iter].confirmed_price = Eigen::MatrixXd(Time, DSO_Markets.markets[DSO_iter].num_zone);
		DSO_Markets.markets[DSO_iter].network.confirmed_power = Eigen::MatrixXd(Time, DSO_Markets.markets[DSO_iter].network.num_edges);		
	}
}

int main(){
	network_inform Power_network_inform;
	power_network_input_process(Power_network_inform);
	
	int Time = 8760;
	std::string fin_name_moc = "input/merit_order_curve_q_assimilated_2021.csv";
	std::string fin_name_demand = "input/residual_load_default_forecast_2021.csv";
	market_inform International_Market;
	International_Market_Set(International_Market, Time, fin_name_moc, fin_name_demand);
	
	market_inform TSO_Market;
	LP_object TSO_Problem; 
	TSO_Market_Set(TSO_Market, Power_network_inform, 1);
	Flow_Based_Market_LP_Set(TSO_Market, TSO_Problem);		
	
	std::string fin_point_demand = "../area_deamand/processed/nominal_mean_demand_field_10km_annual_mean.csv";
	DSO_Markets DSO_Markets;
	DSO_Markets_Set(DSO_Markets, Power_network_inform, 1);
	Submitted_bid_calculation(0, DSO_Markets, TSO_Market, International_Market, Power_network_inform, fin_point_demand);
}