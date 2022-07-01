// Source file for dispatch filtering of DSOs in Norway
#include <iostream>
//#include <chrono>
#include "../basic/LP_gpa.h"
#include "../basic/rw_csv.h"
#include "power_market.h"
#include "IMO.cpp"
#include "TSO.cpp"
#include "../power_network/power_network_input.cpp"

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

void Submitted_bid_calculation(int tick, DSO_Markets &DSO_Markets, market_inform &TSO_Market, market_inform &International_Market, network_inform &Power_network_inform, std::string fin_point_demand){
	auto fin_point_demand_dim = get_file_dim(fin_point_demand);
	auto point_demand_inform = read_file(fin_point_demand_dim[0], fin_point_demand_dim[1], fin_point_demand);
	
	// Initialize submit bids of markets
	Market_Initialization(TSO_Market);
	Market_Initialization(International_Market);
	for(int DSO_iter = 0; DSO_iter < DSO_Markets.markets.size(); ++ DSO_iter){
		Market_Initialization(DSO_Markets.markets[DSO_iter]);
	}

	int bz_ID;
	int DSO_ID;	
	int node_ID;
	int x_ID;
	int y_ID;	
	int point_ID;
	double bid_quan;
	Eigen::VectorXd bid_vec;
	// Trivial case: demand at each point are 100% inflexible
	for(int point_iter = 0; point_iter < Power_network_inform.points.bidding_zone.size(); ++ point_iter){
		node_ID = Power_network_inform.points.node(point_iter);
		DSO_ID = Power_network_inform.nodes.cluster(node_ID);
		bz_ID = Power_network_inform.points.bidding_zone(point_iter);
		bid_quan = point_demand_inform(point_iter, 0) * Power_network_inform.points.population_density(point_iter); //* Power_network_inform.points.point_area;
		// nominal demand currently wrong in processed files, should change them and then multiply area of a point later
		
		DSO_Markets.markets[DSO_ID].submitted_demand(DSO_Markets.markets[DSO_ID].price_intervals + 1, Power_network_inform.points.in_cluster_ID(point_iter)) = bid_quan;
		TSO_Market.submitted_demand(DSO_Markets.markets[DSO_ID].price_intervals + 1, node_ID) += bid_quan;
		International_Market.submitted_demand(DSO_Markets.markets[DSO_ID].price_intervals + 1, bz_ID) += bid_quan;
	}
	
	// Supply at each point (LV power plants) / node (HV power plants)
	for(int hydro_iter = 0; hydro_iter < Power_network_inform.plants.hydro.node.size(); ++ hydro_iter){
		if(Power_network_inform.plants.hydro.cap(hydro_iter) >= 20.){
			node_ID = Power_network_inform.plants.hydro.node(hydro_iter);
			DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);
			bid_vec = International_Market.merit_order_curve.col(Power_network_inform.nodes.bidding_zone(node_ID)) 
				* Power_network_inform.plants.hydro.cap(hydro_iter) / International_Market.merit_order_curve.col(Power_network_inform.nodes.bidding_zone(node_ID)).sum();	
			
			DSO_Markets.markets[DSO_ID].submitted_supply.col(Power_network_inform.DSO_cluster[DSO_ID].points_ID.size() + Power_network_inform.nodes.in_cluster_ID(node_ID)) += bid_vec;
		}
		else{
			x_ID = int((Power_network_inform.plants.hydro.x(hydro_iter) - Power_network_inform.points.x.minCoeff()) / Power_network_inform.points.grid_length + .5);
			y_ID = int((Power_network_inform.plants.hydro.y(hydro_iter) - Power_network_inform.points.y.minCoeff()) / Power_network_inform.points.grid_length + .5);
			point_ID = Power_network_inform.points.coordinate_grid(x_ID, y_ID);
			if(point_ID == -1){
				continue;
			}
			node_ID = Power_network_inform.points.node(point_ID);
			DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);
			bid_vec = International_Market.merit_order_curve.col(Power_network_inform.points.bidding_zone(point_ID))
				* Power_network_inform.plants.hydro.cap(hydro_iter) / International_Market.merit_order_curve.col(Power_network_inform.points.bidding_zone(point_ID)).sum();
			
			DSO_Markets.markets[DSO_ID].submitted_supply.col(Power_network_inform.points.in_cluster_ID(point_ID)) += bid_vec;
		}
		TSO_Market.submitted_supply.col(node_ID) += bid_vec;
		International_Market.submitted_supply.col(bz_ID) += bid_vec;		
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