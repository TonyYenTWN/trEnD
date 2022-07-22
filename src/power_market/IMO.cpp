// Source File for market clearing of the international market operator of the energy-only-market
#include <iostream>
//#include <chrono>
#include "src/basic/rw_csv.h"
#include "power_market.h"

void power_market::International_Market_Set(market_inform &International_Market, power_network::network_inform &Power_network_inform, int Time, std::string fin_name_moc, std::string fin_name_demand){
	// Input Parameters of international market
	International_Market.num_zone = 13;
	International_Market.cross_border_zone_start = 5;
	International_Market.time_intervals = Time;
	International_Market.zone_names = Power_network_inform.cbt.bz_names;
	International_Market.set_bidded_price();
	International_Market.network.num_vertice = International_Market.num_zone;
	International_Market.network.line_capacity_dense = Power_network_inform.cbt.flow_constraint;

	// Construct incidence vector matrix
	International_Market.network.incidence.reserve(International_Market.network.num_vertice * International_Market.network.num_vertice);
	for(int row_iter = 0; row_iter < International_Market.network.num_vertice - 1; ++ row_iter){
		for(int col_iter = row_iter + 1; col_iter < International_Market.network.num_vertice; ++ col_iter){
			bool add_flag = 1 - (International_Market.network.line_capacity_dense(row_iter, col_iter) == 0.) * (International_Market.network.line_capacity_dense(col_iter, row_iter) == 0.);
			if(add_flag){
				International_Market.network.incidence.push_back(Eigen::Vector2i(row_iter, col_iter));
			}
		}
	}

	// Construct power constraint matrix
	International_Market.network.num_edges = International_Market.network.incidence.size();
	International_Market.network.power_constraint = Eigen::MatrixXd(International_Market.network.num_edges, 2);
	for(int edge_iter = 0; edge_iter < International_Market.network.num_edges; ++ edge_iter){
		int row_ID = International_Market.network.incidence[edge_iter](0);
		int col_ID = International_Market.network.incidence[edge_iter](1);
		International_Market.network.power_constraint(edge_iter, 0) = International_Market.network.line_capacity_dense(row_ID, col_ID);
		International_Market.network.power_constraint(edge_iter, 1) = International_Market.network.line_capacity_dense(col_ID, row_ID);
	}

	// Quantity density at each price
	// Read inferred merit order curve data
	int num_row = International_Market.price_intervals + 2;
	int num_col = International_Market.num_zone;
	Eigen::MatrixXd merit_order_curve_q = basic::read_file(num_row, num_col, fin_name_moc);
	Eigen::MatrixXd diff_merit_order_curve_q = merit_order_curve_q.bottomRows(num_row - 1) - merit_order_curve_q.topRows(num_row - 1);
	International_Market.merit_order_curve = merit_order_curve_q;
	International_Market.merit_order_curve.bottomRows(num_row - 1) = diff_merit_order_curve_q;
	International_Market.merit_order_curve = Eigen::MatrixXd::Ones(num_row, num_col).array() * International_Market.merit_order_curve.array().max(0);

	// Default demand time series
	// Read default demand data
	num_row = International_Market.time_intervals;
	num_col = International_Market.num_zone + 1;
	International_Market.demand_default = basic::read_file(num_row, num_col, fin_name_demand).rightCols(International_Market.num_zone);

	// Initialization of process variables
	power_market::Market_Initialization(International_Market);

	// Initialization of output variables
	International_Market.confirmed_supply = Eigen::MatrixXd::Zero(Time, International_Market.num_zone);
	International_Market.confirmed_demand = Eigen::MatrixXd::Zero(Time, International_Market.num_zone);
	International_Market.confirmed_price = Eigen::MatrixXd(Time, International_Market.num_zone);
	International_Market.network.confirmed_power = Eigen::MatrixXd(Time, International_Market.network.num_edges);
}

void power_market::International_Market_Optimization(int tick, market_inform &International_Market, bool print_result){
	// Initialization of process variables
	int type_capacity_exchange;
	double exchange_quantity;
	double maximum_price_diff;
	Eigen::Vector2i maximum_price_diff_ID;
	Eigen::VectorXi default_price_ID;
	Eigen::VectorXi price_demand_ID;
	Eigen::VectorXi price_supply_ID;
	Eigen::MatrixXd bidded_supply;
	Eigen::MatrixXd bidded_demand;
	Eigen::MatrixXd bidded_supply_export = Eigen::MatrixXd::Zero(International_Market.price_intervals + 2, International_Market.num_zone);
	Eigen::MatrixXd bidded_demand_export = Eigen::MatrixXd::Zero(International_Market.price_intervals + 2, International_Market.num_zone);
	Eigen::MatrixXd bidded_supply_import = Eigen::MatrixXd::Zero(International_Market.price_intervals + 2, International_Market.num_zone);
	Eigen::MatrixXd bidded_demand_import = Eigen::MatrixXd::Zero(International_Market.price_intervals + 2, International_Market.num_zone);
	Eigen::MatrixXd maximum_capacity_exchange = International_Market.network.line_capacity_dense;
	Eigen::MatrixXi available_capacity_exchange;
	Eigen::MatrixXd remaining_capacity_exchange;
	Eigen::MatrixXd surplus_exchange;
	Eigen::MatrixXd actual_capacity_exchange;

	// Update of submitted supply and demand bids at other nations
	// Bidding zones of other nations: assume inflexible supply or demand
	International_Market.submitted_supply = International_Market.merit_order_curve;
	for(int zone_ID = International_Market.cross_border_zone_start; zone_ID < International_Market.num_zone; ++ zone_ID){
		if(International_Market.demand_default(tick, zone_ID) >= 0){
			International_Market.submitted_demand(International_Market.price_intervals, zone_ID) += International_Market.demand_default(tick, zone_ID);
		}
		else{
			International_Market.submitted_supply(0, zone_ID) += -International_Market.demand_default(tick, zone_ID);
			//std::cout << "Negative residual demand!!! \n\n";
		}
	}
	bidded_supply = International_Market.submitted_supply;
	bidded_demand = International_Market.submitted_demand;

	// Initial market clearing within each bidding zones
	default_price_ID = Eigen::VectorXi(International_Market.num_zone);
	power_market::Market_clearing_nodal(tick, International_Market, default_price_ID, bidded_supply, bidded_demand);

//	if(print_result){
//		std::cout << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n";
//		std::cout << "Tick: " << tick << "\n";
//		std::cout << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n";
//		std::cout << "  Default Price: " << default_price_ID.transpose() << "\n";
//		std::cout << "  Sell Quantity: " << International_Market.confirmed_supply.row(tick) << "\n";
//		std::cout << "   Buy Quantity: " << International_Market.confirmed_demand.row(tick) << "\n";
//		std::cout << "Residual Demand: " << bidded_demand.bottomRows(1) << "\n\n";
//	}

	// Optimization of cross border exchange
	for(int zone_ID = 0; zone_ID < International_Market.num_zone; ++ zone_ID){
		// Supply curves
		bidded_supply_export.col(zone_ID).array().tail(International_Market.price_intervals + 1 - default_price_ID(zone_ID)) = International_Market.submitted_supply.col(zone_ID).array().tail(International_Market.price_intervals + 1 - default_price_ID(zone_ID));
		bidded_supply_export(default_price_ID(zone_ID), zone_ID) = bidded_supply(default_price_ID(zone_ID), zone_ID);
		bidded_supply_import.col(zone_ID).array().head(default_price_ID(zone_ID)) = International_Market.submitted_supply.col(zone_ID).array().head(default_price_ID(zone_ID));
		bidded_supply_import(default_price_ID(zone_ID), zone_ID) = International_Market.submitted_supply(default_price_ID(zone_ID), zone_ID) - bidded_supply(default_price_ID(zone_ID), zone_ID);

		// Demand curves
		bidded_demand_export.col(zone_ID).array().tail(International_Market.price_intervals + 1 - default_price_ID(zone_ID)) = International_Market.submitted_demand.col(zone_ID).array().tail(International_Market.price_intervals + 1 - default_price_ID(zone_ID));
		bidded_demand_export(default_price_ID(zone_ID), zone_ID) = International_Market.submitted_demand(default_price_ID(zone_ID), zone_ID) - bidded_demand(default_price_ID(zone_ID), zone_ID);
		bidded_demand_import.col(zone_ID).array().head(default_price_ID(zone_ID)) = International_Market.submitted_demand.col(zone_ID).array().head(default_price_ID(zone_ID));
		bidded_demand_import(default_price_ID(zone_ID), zone_ID) = bidded_demand(default_price_ID(zone_ID), zone_ID);
	}

	price_demand_ID = default_price_ID;
	price_supply_ID = default_price_ID;
	remaining_capacity_exchange = maximum_capacity_exchange;
	surplus_exchange = (International_Market.bidded_price(International_Market.bidded_price.size() - 1) - International_Market.bidded_price(0)) * Eigen::MatrixXd::Ones(International_Market.num_zone, International_Market.num_zone);
	available_capacity_exchange = Eigen::MatrixXi::Ones(International_Market.num_zone, International_Market.num_zone);
	available_capacity_exchange.bottomRightCorner(International_Market.num_zone - International_Market.cross_border_zone_start, International_Market.num_zone - International_Market.cross_border_zone_start) = Eigen::MatrixXi::Identity(International_Market.num_zone - International_Market.cross_border_zone_start, International_Market.num_zone - International_Market.cross_border_zone_start);

	// Main loop for optimization of cross border exchange
	int loop_count = 0;
	//while(loop_count < 100){
	while(available_capacity_exchange.sum() > 0){
		loop_count += 1;
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
//		for(int row_ID = 0; row_ID < International_Market.num_zone; ++ row_ID){
//			for(int col_ID = 0; col_ID < International_Market.num_zone; ++ col_ID){
		for(int row_ID = International_Market.num_zone - 1; row_ID > -1; -- row_ID){
			for(int col_ID = International_Market.num_zone - 1; col_ID > - 1; -- col_ID){
				if(available_capacity_exchange(row_ID, col_ID)){
					// Check if limit of edge capacity is reached
					if(remaining_capacity_exchange(row_ID, col_ID) == 0){
						available_capacity_exchange(row_ID, col_ID) = 0;
						// Skip the comparison below
						continue;
					}

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
								//std::cout << type_capacity_exchange << " " << row_ID << " " << col_ID << "\n";
								//std::cout << price_demand_ID(col_ID) << " " << price_supply_ID(row_ID) << "\n";
							}

							// Encode maximum surplus and occuring zones
							maximum_price_diff = surplus_exchange(row_ID, col_ID);
							maximum_price_diff_ID << row_ID, col_ID;
						}
					}
					else{
						available_capacity_exchange(row_ID, col_ID) = 0;
					}

//					// Check if limit of edge capacity is reached
//					if(remaining_capacity_exchange(row_ID, col_ID) == 0){
//						available_capacity_exchange(row_ID, col_ID) = 0;
//					}
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
					//std::cout << bidded_supply_export(price_supply_ID(maximum_price_diff_ID(0)), maximum_price_diff_ID(0)) << "\n";

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
					//std::cout << bidded_supply_export(price_supply_ID(maximum_price_diff_ID(0)), maximum_price_diff_ID(0)) << "\n";
					//std::cout << bidded_supply_import(price_demand_ID(maximum_price_diff_ID(1)), maximum_price_diff_ID(1)) << "\n";
					exchange_quantity = std::min(bidded_supply_import(price_demand_ID(maximum_price_diff_ID(1)), maximum_price_diff_ID(1)), bidded_supply_export(price_supply_ID(maximum_price_diff_ID(0)), maximum_price_diff_ID(0)));
					//std::cout << exchange_quantity << "\n";
					exchange_quantity = std::min(exchange_quantity, remaining_capacity_exchange(maximum_price_diff_ID(0), maximum_price_diff_ID(1)));
					//std::cout << exchange_quantity << "\n";
					International_Market.confirmed_supply(tick, maximum_price_diff_ID(0)) += exchange_quantity;
					International_Market.confirmed_supply(tick, maximum_price_diff_ID(1)) -= exchange_quantity;
					bidded_supply_import(price_demand_ID(maximum_price_diff_ID(1)), maximum_price_diff_ID(1)) -= exchange_quantity;
					bidded_supply_export(price_supply_ID(maximum_price_diff_ID(0)), maximum_price_diff_ID(0)) -= exchange_quantity;
					//std::cout << bidded_supply_export(price_supply_ID(maximum_price_diff_ID(0)), maximum_price_diff_ID(0)) << "\n";
					//std::cout << bidded_supply_import(price_demand_ID(maximum_price_diff_ID(1)), maximum_price_diff_ID(1)) << "\n";

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
			//std::cout << remaining_capacity_exchange(maximum_price_diff_ID(0), maximum_price_diff_ID(1)) << "\n";
			remaining_capacity_exchange(maximum_price_diff_ID(0), maximum_price_diff_ID(1)) -= exchange_quantity;
			//std::cout << remaining_capacity_exchange(maximum_price_diff_ID(0), maximum_price_diff_ID(1)) << "\n";
		}

//		std::cout << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n";
//		std::cout << "Tick: " << tick << "\n";
//		std::cout << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n";
//		std::cout << "  Export Price ID: " << price_supply_ID.transpose() << "\n";
//		std::cout << "  Import Price ID: " << price_demand_ID.transpose() << "\n";
//		std::cout << "    Sell Quantity: " << International_Market.confirmed_supply.row(tick) << "\n";
//		std::cout << "     Buy Quantity: " << International_Market.confirmed_demand.row(tick) << "\n";
//		std::cout << type_capacity_exchange << ": " << maximum_price_diff_ID.transpose() << "\n\n";
	}

	// Calculate exchange flow and market clearing prices of the bidding zones
	actual_capacity_exchange = maximum_capacity_exchange - remaining_capacity_exchange;

	// Record the cross-border flow
	for(int edge_ID = 0; edge_ID < International_Market.network.num_edges; ++ edge_ID){
		if(actual_capacity_exchange(International_Market.network.incidence[edge_ID](0), International_Market.network.incidence[edge_ID](1)) > 0){
			International_Market.network.confirmed_power(tick, edge_ID) = actual_capacity_exchange(International_Market.network.incidence[edge_ID](0), International_Market.network.incidence[edge_ID](1));
		}
		else{
			International_Market.network.confirmed_power(tick, edge_ID) = -actual_capacity_exchange(International_Market.network.incidence[edge_ID](1), International_Market.network.incidence[edge_ID](0));
		}
	}

	// Results Output
	if(print_result){
		std::cout << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n";
		std::cout << "Tick: " << tick << "\n";
		std::cout << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n";
		std::cout << "  Export Price ID: " << price_supply_ID.transpose() << "\n";
		std::cout << "  Import Price ID: " << price_demand_ID.transpose() << "\n";
		std::cout << "     Market Price: " << International_Market.confirmed_price.row(tick) << "\n";
		std::cout << "    Sell Quantity: " << International_Market.confirmed_supply.row(tick) << "\n";
		std::cout << "     Buy Quantity: " << International_Market.confirmed_demand.row(tick) << "\n";
		std::cout << "Cross border flow: " << International_Market.network.confirmed_power.row(tick) << "\n\n";
	}
}

void power_market::International_Market_Output(market_inform &International_Market){
	std::string fout_name = "output/IMO_confirmed_price.csv";
	basic::write_file(International_Market.confirmed_price, fout_name, International_Market.zone_names);
}

//int main(){
//	// Input variables
//	int Time = 8760;
//	std::string fin_name_moc = "input/merit_order_curve_q_assimilated_2021.csv";
//	std::string fin_name_demand = "input/residual_load_default_forecast_2021.csv";
//	market_inform International_Market;
//	International_Market_Set(International_Market, Time, fin_name_moc, fin_name_demand);
//
//	// Naive market clearing
//	for(int tick = 0; tick < 10; ++ tick){
//		Market_Initialization(International_Market);
//		International_Market_Optimization(tick, International_Market, 1);
//	}
//
//	// Write csv file
//	// International_Market_Output(International_Market);
//}
