// Source file for dispatch filtering of DSOs in Norway
#include <iostream>
//#include <chrono>
#include "../basic/LP_gpa.h"
#include "power_market.cpp"
#include "../power_network/power_network_input.cpp"
#include "../basic/rw_csv.cpp"

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

void DSO_Markets_submitted_bid_calculation(int tick, DSO_Markets &DSO_Markets, network_inform &Power_network_inform, std::string fin_point_demand){
	auto fin_point_demand_dim = get_file_dim(fin_point_demand);
	auto point_demand_inform = read_file(fin_point_demand_dim[0], fin_point_demand_dim[1], fin_point_demand);
	
	// Initialize submit bids of clustered DSOs
	for(int DSO_iter = 0; DSO_iter < DSO_Markets.markets.size(); ++ DSO_iter){
		Market_Initialization(DSO_Markets.markets[DSO_iter]);
	}
	
	// Trivial case: demand at each point are 100% inflexible
	int DSO_ID;
	for(int point_iter = 0; point_iter < Power_network_inform.points.bidding_zone.size(); ++ point_iter){
		DSO_ID = Power_network_inform.nodes.cluster(Power_network_inform.points.node(point_iter));
		DSO_Markets.markets[DSO_ID].submitted_demand(DSO_Markets.markets[DSO_ID].price_intervals + 1, Power_network_inform.points.in_cluster_ID(point_iter))
			= point_demand_inform(point_iter, 0) * Power_network_inform.points.population_density(point_iter); //* Power_network_inform.points.point_area;
			// nominal demand currently wrong in processed files, should change them and then multiply area of a point later
	}
}

int main(){
	network_inform Power_network_inform;
	power_network_input_process(Power_network_inform);	
	
	std::string fin_point_demand = "../area_deamand/processed/nominal_mean_demand_field_10km_annual_mean.csv";
	DSO_Markets DSO_Markets;
	DSO_Markets_Set(DSO_Markets, Power_network_inform, 1);
	DSO_Markets_submitted_bid_calculation(0, DSO_Markets, Power_network_inform, fin_point_demand);
}