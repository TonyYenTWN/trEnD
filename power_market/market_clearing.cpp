// Main Source File
#include <iostream>
//#include <chrono>
#include <omp.h>
#include "../basic/rw_csv.cpp"

struct network_graph{
	// Input parameters
	int num_vertice;
	int num_edges;
	// Compact Incidence Matrix
	Eigen::MatrixXi incidence_matrix;		// 0th col: start; 1st col: end
	// Power flow constraint
	Eigen::MatrixXd power_constraint;  		// 0th col: from start to end; 1st col: from end to start

	// Output variables
	Eigen::MatrixXd confirmed_power;		// Power flow across each edge; positive indicates flowing from start to end
	Eigen::MatrixXd confirmed_voltage;		// Voltage at each vertice
};

struct market_inform{
	// Input parameters
	int num_zone;							// Can be the actual bidding zones, or just a node / spatial element
	int cross_border_zone_start;			// Indicate the index, after whose corresponding bidding zones are on the boundary and only act by cross-border flow
	int time_intervals;
	int price_intervals;
	std::vector<std::string> zone_names;
	Eigen::Vector2d price_range_inflex;
	Eigen::Vector2d price_range_flex;
	Eigen::VectorXd bidded_price;
	Eigen::MatrixXd merit_order_curve;
	//Eigen::MatrixXd cross_border_price;		// Default price of cross border zones without cross-border flow with bidding zones inside the boundary
	
	// Output variables
	Eigen::MatrixXd confirmed_supply;		// Confirmed supply quantity from MO or TSO 
	Eigen::MatrixXd confirmed_demand;		// Confirmed demand quantity from MO or TSO 
	Eigen::MatrixXd confirmed_price;		// Confirmed market clearing price from MO or TSO 
	Eigen::MatrixXd confirmed_bid_supply; 	// Confirmed supply bids from DSOs forwarding to wholesale electricity market
	Eigen::MatrixXd confirmed_bid_demand;	// Confirmed demand bids from DSOs forwarding to wholesale electricity market
	
	// Mixed Substructure
	network_graph network;
};

market_inform International_Market_Set(int Time, std::string fin_name_moc){
	market_inform International_Market;
	
	// Input Parameters of international market
	International_Market.num_zone = 13;
	International_Market.cross_border_zone_start = 5;
	International_Market.time_intervals = Time;
	International_Market.price_intervals = 600;
	International_Market.zone_names = {"NO1","NO2","NO3","NO4","NO5","DE-LU","DK1","FI","GB","NL","SE1","SE2","SE3"};
	International_Market.price_range_inflex << -500, 3000;
	International_Market.price_range_flex << -100, 500;
	International_Market.network.num_vertice = International_Market.num_zone;
	International_Market.network.num_edges = 15;
	International_Market.network.incidence_matrix = Eigen::MatrixXi(International_Market.network.num_edges, 2);
	International_Market.network.incidence_matrix.col(0) << 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3;
	International_Market.network.incidence_matrix.col(1) << 1, 2, 4, 12, 4, 5, 6, 8, 9, 3, 4, 11, 7, 10, 11;
	International_Market.network.power_constraint = Eigen::MatrixXd(International_Market.network.num_edges, 2);
	International_Market.network.power_constraint.col(0) << 1900, 100, 500, 2130, 300, 1400, 1680, 720, 720, 300, 500, 600, 0, 650, 200;
	International_Market.network.power_constraint.col(1) << 3400, 350, 3900, 2095, 500, 1400, 1150, 720, 720, 1100, 450, 1000, 0, 600, 250;
	International_Market.bidded_price = Eigen::VectorXd(International_Market.price_intervals + 2);
	International_Market.bidded_price(0) = International_Market.price_range_inflex(0);
	International_Market.bidded_price.array().tail(1) = International_Market.price_range_inflex(1);
	International_Market.bidded_price.array().segment(1, International_Market.price_intervals) = Eigen::VectorXd::LinSpaced(International_Market.price_intervals, International_Market.price_range_flex(0) + .5, International_Market.price_range_flex(1) - .5);
	
	// Quantity density at each price
	// Read inferred merit order curve data
	int num_row = International_Market.price_intervals + 2; 
	int num_col = International_Market.num_zone;
	Eigen::MatrixXd merit_order_curve_q = read_file(num_row, num_col, fin_name_moc);
	Eigen::MatrixXd diff_merit_order_curve_q = merit_order_curve_q.bottomRows(num_row - 1) - merit_order_curve_q.topRows(num_row - 1);
	International_Market.merit_order_curve = merit_order_curve_q;
	International_Market.merit_order_curve.bottomRows(num_row - 1) = diff_merit_order_curve_q;
	International_Market.merit_order_curve = Eigen::MatrixXd::Ones(num_row, num_col).array() * International_Market.merit_order_curve.array().max(0);
	
	// Initialization of output variables
	International_Market.confirmed_supply = Eigen::MatrixXd(Time, International_Market.num_zone);
	International_Market.confirmed_demand = Eigen::MatrixXd(Time, International_Market.num_zone);
	International_Market.confirmed_price = Eigen::MatrixXd(Time, International_Market.num_zone);
	International_Market.network.confirmed_power = Eigen::MatrixXd(Time, International_Market.network.num_edges);

	return(International_Market);
}

void International_Market_Optimization(int tick, market_inform &International_Market){
	// Initialization of process variables
	int type_capacity_exchange;
	double exchange_quantity;
	double maximum_price_diff;
	Eigen::Vector2i maximum_price_diff_ID;
	Eigen::VectorXi default_price_ID;
	Eigen::VectorXi price_demand_ID;
	Eigen::VectorXi price_supply_ID;
	Eigen::VectorXd demand_inflex(International_Market.num_zone);
	demand_inflex << 4500, 4000, 3500, 2500, 2500, 15000, 500, 2000, 10000, 5000, 5000, 5000, 10000;
	Eigen::MatrixXd bidded_supply_default = Eigen::MatrixXd::Zero(International_Market.price_intervals + 2, International_Market.num_zone);
	Eigen::MatrixXd bidded_demand_default = Eigen::MatrixXd::Zero(International_Market.price_intervals + 2, International_Market.num_zone);	
	Eigen::MatrixXd bidded_supply = Eigen::MatrixXd::Zero(International_Market.price_intervals + 2, International_Market.num_zone);
	Eigen::MatrixXd bidded_demand = Eigen::MatrixXd::Zero(International_Market.price_intervals + 2, International_Market.num_zone);	
	Eigen::MatrixXd bidded_supply_export = Eigen::MatrixXd::Zero(International_Market.price_intervals + 2, International_Market.num_zone);
	Eigen::MatrixXd bidded_demand_export = Eigen::MatrixXd::Zero(International_Market.price_intervals + 2, International_Market.num_zone);
	Eigen::MatrixXd bidded_supply_import = Eigen::MatrixXd::Zero(International_Market.price_intervals + 2, International_Market.num_zone);
	Eigen::MatrixXd bidded_demand_import = Eigen::MatrixXd::Zero(International_Market.price_intervals + 2, International_Market.num_zone);
	Eigen::MatrixXd maximum_capacity_exchange = Eigen::MatrixXd::Zero(International_Market.num_zone, International_Market.num_zone);
	for(int edge_ID = 0; edge_ID < International_Market.network.num_edges; ++ edge_ID){
		maximum_capacity_exchange(International_Market.network.incidence_matrix(edge_ID, 0), International_Market.network.incidence_matrix(edge_ID, 1)) = International_Market.network.power_constraint(edge_ID, 0);
		maximum_capacity_exchange(International_Market.network.incidence_matrix(edge_ID, 1), International_Market.network.incidence_matrix(edge_ID, 0)) = International_Market.network.power_constraint(edge_ID, 1);
	}
	Eigen::MatrixXi available_capacity_exchange;
	Eigen::MatrixXd remaining_capacity_exchange;
	Eigen::MatrixXd surplus_exchange;
	Eigen::MatrixXd actual_capacity_exchange;
	
	// Execution of optimization process
	bidded_supply_default = International_Market.merit_order_curve;
	bidded_demand_default.bottomRows(1) = demand_inflex.transpose();
	bidded_supply = bidded_supply_default;
	bidded_demand = bidded_demand_default;
	
	// Initial optimization within each bidding zones
	default_price_ID = Eigen::VectorXi(International_Market.num_zone);
	price_demand_ID = (International_Market.price_intervals + 1) * Eigen::VectorXi::Ones(International_Market.num_zone);
	price_supply_ID = Eigen::VectorXi::Zero(International_Market.num_zone);
	
	double trade_quantity;
	for(int zone_ID = 0; zone_ID < International_Market.num_zone; ++ zone_ID){
		
		while(price_demand_ID(zone_ID) > price_supply_ID(zone_ID)){
			// Check if there are demand bids at current price interval
			while(bidded_demand(price_demand_ID(zone_ID), zone_ID) == 0){
				if(price_demand_ID(zone_ID) > 0){
					price_demand_ID(zone_ID) -= 1;
				}
				else{
					// No available buyer left to buy electricity
					default_price_ID(zone_ID) = price_supply_ID(zone_ID);
					break;
				}
			}
			
			// Check if there are supply bids at current price interval
			while(bidded_supply(price_supply_ID(zone_ID), zone_ID) == 0){
				if(price_supply_ID(zone_ID) < International_Market.bidded_price.size() - 1){
					price_supply_ID(zone_ID) += 1;
				}
				else{
					// No available seller left to sell electricity
					default_price_ID(zone_ID) = price_demand_ID(zone_ID);
					break;				
				}			
			}
			
			if(price_demand_ID(zone_ID) > price_supply_ID(zone_ID)){
				trade_quantity = std::min(bidded_supply(price_supply_ID(zone_ID), zone_ID), bidded_demand(price_demand_ID(zone_ID), zone_ID));
				International_Market.confirmed_supply(tick, zone_ID) += trade_quantity;
				International_Market.confirmed_demand(tick, zone_ID) += trade_quantity;
				bidded_supply(price_supply_ID(zone_ID), zone_ID) -= trade_quantity;
				bidded_demand(price_demand_ID(zone_ID), zone_ID) -= trade_quantity;
			}
			else{
				if(bidded_supply(price_supply_ID(zone_ID), zone_ID) >= 0){
					default_price_ID(zone_ID) = price_supply_ID(zone_ID);
				}
				else{
					default_price_ID(zone_ID) = price_demand_ID(zone_ID);
				}
			}
		}
		
		// Record default market clearing prices
		International_Market.confirmed_price(tick, zone_ID) = International_Market.bidded_price(default_price_ID(zone_ID));		
	}
	
	std::cout << "  Default Price: " << default_price_ID.transpose() << std::endl;
	std::cout << "  Sell Quantity: " << International_Market.confirmed_supply << std::endl;
	std::cout << "   Buy Quantity: " << International_Market.confirmed_demand << std::endl;
	std::cout << "Residual Demand: " << bidded_demand.bottomRows(1) << "\n" << std::endl;
	
	// Optimization of cross border exchange
	for(int zone_ID = 0; zone_ID < International_Market.num_zone; ++ zone_ID){
		// Supply curves
		bidded_supply_export.col(zone_ID).array().tail(International_Market.price_intervals + 1 - default_price_ID(zone_ID)) = bidded_supply_default.col(zone_ID).array().tail(International_Market.price_intervals + 1 - default_price_ID(zone_ID));
		bidded_supply_export(default_price_ID(zone_ID), zone_ID) = bidded_supply(default_price_ID(zone_ID), zone_ID);
		bidded_supply_import.col(zone_ID).array().head(default_price_ID(zone_ID)) = bidded_supply_default.col(zone_ID).array().head(default_price_ID(zone_ID));
		bidded_supply_import(default_price_ID(zone_ID), zone_ID) = bidded_supply_default(default_price_ID(zone_ID), zone_ID) - bidded_supply(default_price_ID(zone_ID), zone_ID);
		
		// Demand curves: inflexible case
		bidded_demand_export.col(zone_ID).array().tail(International_Market.price_intervals + 1 - default_price_ID(zone_ID)) = bidded_demand_default.col(zone_ID).array().tail(International_Market.price_intervals + 1 - default_price_ID(zone_ID));
		bidded_demand_export(default_price_ID(zone_ID), zone_ID) = bidded_demand_default(default_price_ID(zone_ID), zone_ID) - bidded_demand(default_price_ID(zone_ID), zone_ID);
		bidded_demand_import.col(zone_ID).array().head(default_price_ID(zone_ID)) = bidded_demand_default.col(zone_ID).array().head(default_price_ID(zone_ID));
		bidded_demand_import(default_price_ID(zone_ID), zone_ID) = bidded_demand(default_price_ID(zone_ID), zone_ID);		
	}
	
	price_demand_ID = default_price_ID;
	price_supply_ID = default_price_ID;
	remaining_capacity_exchange = maximum_capacity_exchange;
	surplus_exchange = (International_Market.bidded_price(International_Market.bidded_price.size() - 1) - International_Market.bidded_price(0)) * Eigen::MatrixXd::Ones(International_Market.num_zone, International_Market.num_zone);
	available_capacity_exchange = Eigen::MatrixXi::Ones(International_Market.num_zone, International_Market.num_zone);
	available_capacity_exchange.bottomRightCorner(International_Market.num_zone - International_Market.cross_border_zone_start, International_Market.num_zone - International_Market.cross_border_zone_start) = Eigen::MatrixXi::Identity(International_Market.num_zone - International_Market.cross_border_zone_start, International_Market.num_zone - International_Market.cross_border_zone_start);
	
	// Main loop for optimization
	while(available_capacity_exchange.sum() > 0){
		// Update price of demand and supply at each bidding zone
		for(int zone_ID = 0; zone_ID < International_Market.num_zone; ++ zone_ID){
			// Check if there are demand bids at current price interval
			while(bidded_demand_import(price_demand_ID(zone_ID), zone_ID) == 0 && bidded_supply_import(price_demand_ID(zone_ID), zone_ID) == 0){
				if(price_demand_ID(zone_ID) > 0){
					price_demand_ID(zone_ID) -= 1;
				}
				else{
					// No available capacity left to import electricity
					available_capacity_exchange.col(zone_ID) = Eigen::VectorXi::Zero(International_Market.num_zone);
					break;
				}
			}
			
			// Check if there are supply bids at current price interval
			while(bidded_supply_export(price_supply_ID(zone_ID), zone_ID) == 0 && bidded_demand_export(price_supply_ID(zone_ID), zone_ID) == 0){
				if(price_supply_ID(zone_ID) < International_Market.bidded_price.size() - 1){
					price_supply_ID(zone_ID) += 1;
				}
				else{
					// No available capacity left to export electricity
					available_capacity_exchange.row(zone_ID) = Eigen::VectorXi::Zero(International_Market.num_zone);
					break;				
				}			
			}
		}
		
		// Find the exchange with the greatest surplus
		maximum_price_diff = 0;
		for(int row_ID = 0; row_ID < International_Market.num_zone; ++ row_ID){
			for(int col_ID = 0; col_ID < International_Market.num_zone; ++ col_ID){
				if(available_capacity_exchange(row_ID, col_ID)){
					// Check if surplus according to updated price is still positive
					if(price_demand_ID(col_ID) >= price_supply_ID(row_ID)){
						// Update surplus for each possible exchange configuration
						surplus_exchange(row_ID, col_ID) = International_Market.bidded_price(price_demand_ID(col_ID)) - International_Market.bidded_price(price_supply_ID(row_ID));
						
						// Check if the surplus is the current maximum
						if(surplus_exchange(row_ID, col_ID) >= maximum_price_diff){
							// Find the type of exchange
							if(bidded_supply_export(price_supply_ID(row_ID), row_ID) > 0){
								if(bidded_demand_import(price_supply_ID(col_ID), col_ID) > 0){
									type_capacity_exchange = 0;			// Increase total trade quantity
								}
								else{
									type_capacity_exchange = 1;			// Replace supply in importing zone with export
								}
							}
							else{
								type_capacity_exchange = 2;				// Replace demand in exporting zone with import
								std::cout << type_capacity_exchange << " " << row_ID << " " << col_ID << std::endl;
								std::cout << price_demand_ID(col_ID) << " " << price_supply_ID(row_ID) << std::endl;
							}
							
							// Encode maximum surplus and occuring zones
							maximum_price_diff = surplus_exchange(row_ID, col_ID);
							maximum_price_diff_ID << row_ID, col_ID;
						}					
					}
					else{	
						available_capacity_exchange(row_ID, col_ID) = 0;
					}
					
					// Check if limit of edge capacity is reached
					if(remaining_capacity_exchange(row_ID, col_ID) == 0){
						available_capacity_exchange(row_ID, col_ID) = 0;
					}
				}
			}
		}
		
		// Exchange between bidding zones
		if(maximum_price_diff_ID(0) != maximum_price_diff_ID(1)){
			switch(type_capacity_exchange){
				case 0:
					// Update traded quantity
					exchange_quantity = std::min(bidded_demand_import(price_demand_ID(maximum_price_diff_ID(1)), maximum_price_diff_ID(1)), bidded_supply_export(price_supply_ID(maximum_price_diff_ID(0)), maximum_price_diff_ID(0)));
					exchange_quantity = std::min(exchange_quantity, remaining_capacity_exchange(maximum_price_diff_ID(0), maximum_price_diff_ID(1)));
					International_Market.confirmed_supply(tick, maximum_price_diff_ID(0)) += exchange_quantity;
					International_Market.confirmed_demand(tick, maximum_price_diff_ID(1)) += exchange_quantity;
					bidded_demand_import(price_demand_ID(maximum_price_diff_ID(1)), maximum_price_diff_ID(1)) -= exchange_quantity;
					bidded_supply_export(price_supply_ID(maximum_price_diff_ID(0)), maximum_price_diff_ID(0)) -= exchange_quantity;
					
					// Update market clearing price
					if(exchange_quantity == remaining_capacity_exchange(maximum_price_diff_ID(0), maximum_price_diff_ID(1))){
						International_Market.confirmed_price(tick, maximum_price_diff_ID(0)) == International_Market.bidded_price(price_supply_ID(maximum_price_diff_ID(0)));
						International_Market.confirmed_price(tick, maximum_price_diff_ID(1)) == International_Market.bidded_price(price_demand_ID(maximum_price_diff_ID(1)));
					}
					else{
						if(bidded_demand_import(price_demand_ID(maximum_price_diff_ID(1)), maximum_price_diff_ID(1)) <= bidded_supply_export(price_supply_ID(maximum_price_diff_ID(0)), maximum_price_diff_ID(0))){
							International_Market.confirmed_price(tick, maximum_price_diff_ID(0)) = International_Market.bidded_price(price_supply_ID(maximum_price_diff_ID(0)));
							International_Market.confirmed_price(tick, maximum_price_diff_ID(1)) = International_Market.bidded_price(price_supply_ID(maximum_price_diff_ID(0)));	
						}
						else{
							International_Market.confirmed_price(tick, maximum_price_diff_ID(0)) = International_Market.bidded_price(price_demand_ID(maximum_price_diff_ID(1)));
							International_Market.confirmed_price(tick, maximum_price_diff_ID(1)) = International_Market.bidded_price(price_demand_ID(maximum_price_diff_ID(1)));						
						}
					}
					break;
				case 1:
					// Update traded quantity
					exchange_quantity = std::min(bidded_supply_import(price_demand_ID(maximum_price_diff_ID(1)), maximum_price_diff_ID(1)), bidded_supply_export(price_supply_ID(maximum_price_diff_ID(0)), maximum_price_diff_ID(0)));
					exchange_quantity = std::min(exchange_quantity, remaining_capacity_exchange(maximum_price_diff_ID(0), maximum_price_diff_ID(1)));
					International_Market.confirmed_supply(tick, maximum_price_diff_ID(0)) += exchange_quantity;
					International_Market.confirmed_supply(tick, maximum_price_diff_ID(1)) -= exchange_quantity;
					bidded_supply_import(price_demand_ID(maximum_price_diff_ID(1)), maximum_price_diff_ID(1)) -= exchange_quantity;
					bidded_supply_export(price_supply_ID(maximum_price_diff_ID(0)), maximum_price_diff_ID(0)) -= exchange_quantity;		
					
					// Update market clearing price
					if(exchange_quantity == remaining_capacity_exchange(maximum_price_diff_ID(0), maximum_price_diff_ID(1))){
						International_Market.confirmed_price(tick, maximum_price_diff_ID(0)) == International_Market.bidded_price(price_supply_ID(maximum_price_diff_ID(0)));
						International_Market.confirmed_price(tick, maximum_price_diff_ID(1)) == International_Market.bidded_price(price_demand_ID(maximum_price_diff_ID(1)));
					}
					else{
						if(bidded_supply_import(price_demand_ID(maximum_price_diff_ID(1)), maximum_price_diff_ID(1)) <= bidded_supply_export(price_supply_ID(maximum_price_diff_ID(0)), maximum_price_diff_ID(0))){
							International_Market.confirmed_price(tick, maximum_price_diff_ID(0)) = International_Market.bidded_price(price_supply_ID(maximum_price_diff_ID(0)));
							International_Market.confirmed_price(tick, maximum_price_diff_ID(1)) = International_Market.bidded_price(price_supply_ID(maximum_price_diff_ID(0)));	
						}
						else{
							International_Market.confirmed_price(tick, maximum_price_diff_ID(0)) = International_Market.bidded_price(price_demand_ID(maximum_price_diff_ID(1)));
							International_Market.confirmed_price(tick, maximum_price_diff_ID(1)) = International_Market.bidded_price(price_demand_ID(maximum_price_diff_ID(1)));						
						}						
					}					
					break;
				case 2:
					exchange_quantity = std::min(bidded_demand_import(price_demand_ID(maximum_price_diff_ID(1)), maximum_price_diff_ID(1)), bidded_demand_export(price_supply_ID(maximum_price_diff_ID(0)), maximum_price_diff_ID(0)));
					exchange_quantity = std::min(exchange_quantity, remaining_capacity_exchange(maximum_price_diff_ID(0), maximum_price_diff_ID(1)));
					International_Market.confirmed_demand(tick, maximum_price_diff_ID(0)) -= exchange_quantity;
					International_Market.confirmed_demand(tick, maximum_price_diff_ID(1)) += exchange_quantity;
					bidded_demand_import(price_demand_ID(maximum_price_diff_ID(1)), maximum_price_diff_ID(1)) -= exchange_quantity;
					bidded_demand_export(price_supply_ID(maximum_price_diff_ID(0)), maximum_price_diff_ID(0)) -= exchange_quantity;					
					
					// Update market clearing price
					if(exchange_quantity == remaining_capacity_exchange(maximum_price_diff_ID(0), maximum_price_diff_ID(1))){
						International_Market.confirmed_price(tick, maximum_price_diff_ID(0)) == International_Market.bidded_price(price_supply_ID(maximum_price_diff_ID(0)));
						International_Market.confirmed_price(tick, maximum_price_diff_ID(1)) == International_Market.bidded_price(price_demand_ID(maximum_price_diff_ID(1)));
					}
					else{
						if(bidded_demand_import(price_demand_ID(maximum_price_diff_ID(1)), maximum_price_diff_ID(1)) <= bidded_demand_export(price_supply_ID(maximum_price_diff_ID(0)), maximum_price_diff_ID(0))){
							International_Market.confirmed_price(tick, maximum_price_diff_ID(0)) = International_Market.bidded_price(price_supply_ID(maximum_price_diff_ID(0)));
							International_Market.confirmed_price(tick, maximum_price_diff_ID(1)) = International_Market.bidded_price(price_supply_ID(maximum_price_diff_ID(0)));	
						}
						else{
							International_Market.confirmed_price(tick, maximum_price_diff_ID(0)) = International_Market.bidded_price(price_demand_ID(maximum_price_diff_ID(1)));
							International_Market.confirmed_price(tick, maximum_price_diff_ID(1)) = International_Market.bidded_price(price_demand_ID(maximum_price_diff_ID(1)));						
						}
					}					
					break;	
			}
			remaining_capacity_exchange(maximum_price_diff_ID(0), maximum_price_diff_ID(1)) -= exchange_quantity;			
		}
	}

	// Calculate exchange flow and market clearing prices of the bidding zones
	actual_capacity_exchange = maximum_capacity_exchange - remaining_capacity_exchange;

	// Results Output
	std::cout << "Export Price ID: " << price_supply_ID.transpose() << std::endl;
	std::cout << "Import Price ID: " << price_demand_ID.transpose() << std::endl;
	std::cout << "   Market Price: " << International_Market.confirmed_price.row(0) << std::endl;
	std::cout << "  Sell Quantity: " << International_Market.confirmed_supply << std::endl;
	std::cout << "   Buy Quantity: " << International_Market.confirmed_demand << std::endl;
	std::cout << "Residual Demand: " << bidded_demand_default.bottomRows(1) - International_Market.confirmed_demand << "\n" << std::endl;
	std::cout << "    Trade Zones: " << maximum_price_diff_ID.transpose() << "\n" << std::endl;
	std::cout << "        Surplus: " << maximum_price_diff << "\n" << std::endl;
	std::cout << "Available Trade: " << std::endl; 
	std::cout << available_capacity_exchange << "\n" << std::endl;
	std::cout << "Used Capacity: " << std::endl; 
	std::cout << actual_capacity_exchange << "\n" << std::endl;
	std::cout << "Remaining Capacity: " << std::endl; 
	std::cout << remaining_capacity_exchange << "\n" << std::endl;
}

int main(){
	// Input Variables
	int Time = 1;
	std::string fin_name_moc = "input/merit_order_curve_q_assimilated_2021.csv";
	market_inform International_Market = International_Market_Set(Time, fin_name_moc);
	International_Market_Optimization(0, International_Market);
	
	//std::cout << "Supply" << std::endl;
	//std::cout << International_Market.bidded_supply << "\n" << std::endl;
	//std::cout << "Demand" << std::endl;
	//std::cout << International_Market.bidded_demand << std::endl;
}