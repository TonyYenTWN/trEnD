// Source file for printing results from simulation and settlement of the power market
#include <filesystem>
#include <iostream>
#include "src/basic/rw_csv.h"
#include "power_market.h"

namespace local{
	void Simplified_network_print(power_market::market_inform &Market, std::string name){
		int edge_num = Market.network.num_edges;
		std::string fout_name;

		// Set edge names
		std::vector <std::string> col_names(3);
		col_names[0] = "from";
		col_names[1] = "to";
		col_names[2] = "power_capacity";

		Eigen::MatrixXd output(edge_num, 3);
		for(int edge_iter = 0; edge_iter < edge_num; ++ edge_iter){
			output(edge_iter, 0) = Market.network.incidence[edge_iter](0);
			output(edge_iter, 1) = Market.network.incidence[edge_iter](1);
			output(edge_iter, 2) = Market.network.power_constraint(edge_iter, 1);
		}

		fout_name = name + "_simplified_network.csv";
		basic::write_file(output, fout_name, col_names);
	}

	void Market_results_print(power_market::market_inform &Market, std::string name, bool nodal = 1){
		std::string fout_name;

		// Set edge names
		std::vector <std::string> edge_names;
		edge_names.reserve(Market.network.num_edges);
		for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
			std::string edge_name = "Edge_" + std::to_string(edge_iter);
			edge_names.push_back(edge_name);
		}

		// Confirmed results
		fout_name = name + "_confirmed_price.csv";
		basic::write_file(Market.confirmed.price, fout_name, Market.zone_names);
		fout_name = name + "_confirmed_demand.csv";
		basic::write_file(Market.confirmed.demand, fout_name, Market.zone_names);
		fout_name = name + "_confirmed_supply.csv";
		basic::write_file(Market.confirmed.supply, fout_name, Market.zone_names);
		fout_name = name + "_confirmed_DC_power.csv";
		basic::write_file(Market.network.confirmed_power, fout_name, edge_names);

		if(!nodal){
			// Operation schedule
			fout_name = name + "_operation_EOM_cross_border.csv";
			basic::write_file(Market.operation.cross_border.EOM, fout_name, Market.zone_names);
			fout_name = name + "_operation_EOM_end_user.csv";
			basic::write_file(Market.operation.end_user.EOM, fout_name, Market.zone_names);
			fout_name = name + "_operation_EOM_hydro.csv";
			basic::write_file(Market.operation.hydro.EOM, fout_name, Market.zone_names);
			fout_name = name + "_operation_EOM_industrial.csv";
			basic::write_file(Market.operation.industrial.EOM, fout_name, Market.zone_names);
			fout_name = name + "_operation_EOM_pump_storage.csv";
			basic::write_file(Market.operation.pump_storage.EOM, fout_name, Market.zone_names);
			fout_name = name + "_operation_EOM_wind.csv";
			basic::write_file(Market.operation.wind.redispatch, fout_name, Market.zone_names);
			fout_name = name + "_operation_redispatch_cross_border.csv";
			basic::write_file(Market.operation.cross_border.redispatch, fout_name, Market.zone_names);
			fout_name = name + "_operation_redispatch_end_user.csv";
			basic::write_file(Market.operation.end_user.redispatch, fout_name, Market.zone_names);
			fout_name = name + "_operation_redispatch_hydro.csv";
			basic::write_file(Market.operation.hydro.redispatch, fout_name, Market.zone_names);
			fout_name = name + "_operation_redispatch_industrial.csv";
			basic::write_file(Market.operation.industrial.redispatch, fout_name, Market.zone_names);
			fout_name = name + "_operation_redispatch_pump_storage.csv";
			basic::write_file(Market.operation.pump_storage.redispatch, fout_name, Market.zone_names);
			fout_name = name + "_operation_redispatch_wind.csv";
			basic::write_file(Market.operation.wind.redispatch, fout_name, Market.zone_names);
			fout_name = name + "_operation_balancing_cross_border.csv";
			basic::write_file(Market.operation.cross_border.balancing, fout_name, Market.zone_names);
			fout_name = name + "_operation_balancing_end_user.csv";
			basic::write_file(Market.operation.end_user.balancing, fout_name, Market.zone_names);
			fout_name = name + "_operation_balancing_hydro.csv";
			basic::write_file(Market.operation.hydro.balancing, fout_name, Market.zone_names);
			fout_name = name + "_operation_balancing_industrial.csv";
			basic::write_file(Market.operation.industrial.balancing, fout_name, Market.zone_names);
			fout_name = name + "_operation_balancing_pump_storage.csv";
			basic::write_file(Market.operation.pump_storage.balancing, fout_name, Market.zone_names);
			fout_name = name + "_operation_balancing_wind.csv";
			basic::write_file(Market.operation.wind.balancing, fout_name, Market.zone_names);

			return;
		}

		// EOM settlement
		fout_name = name + "_EOM_cost.csv";
		basic::write_file(Market.EOM.cost, fout_name, Market.zone_names);
		fout_name = name + "_EOM_utility.csv";
		basic::write_file(Market.EOM.utility, fout_name, Market.zone_names);

		// Redispatch settlement
		fout_name = name + "_redispatch_cost.csv";
		basic::write_file(Market.redispatch.cost, fout_name, Market.zone_names);
		fout_name = name + "_redispatch_utility.csv";
		basic::write_file(Market.redispatch.utility, fout_name, Market.zone_names);
		fout_name = name + "_redispatch_demand.csv";
		basic::write_file(Market.redispatch.demand_up - Market.redispatch.demand_down, fout_name, Market.zone_names);
		fout_name = name + "_redispatch_supply.csv";
		basic::write_file(Market.redispatch.supply_up - Market.redispatch.supply_down, fout_name, Market.zone_names);

		// Imbalance
		fout_name = name + "_imbalance_demand.csv";
		basic::write_file(Market.imbalance.demand_up - Market.redispatch.demand_down, fout_name, Market.zone_names);
		fout_name = name + "_imbalance_supply.csv";
		basic::write_file(Market.imbalance.supply_up - Market.redispatch.supply_down, fout_name, Market.zone_names);

		// Actual results
		fout_name = name + "_actual_price.csv";
		basic::write_file(Market.actual.price, fout_name, Market.zone_names);
		fout_name = name + "_actual_demand.csv";
		basic::write_file(Market.actual.demand, fout_name, Market.zone_names);
		fout_name = name + "_actual_supply.csv";
		basic::write_file(Market.actual.supply, fout_name, Market.zone_names);
		fout_name = name + "_actual_DC_power.csv";
		basic::write_file(Market.network.actual_power, fout_name, edge_names);

		// Balancing settlement
		fout_name = name + "_balancing_cost.csv";
		basic::write_file(Market.balancing.cost, fout_name, Market.zone_names);
		fout_name = name + "_balancing_utility.csv";
		basic::write_file(Market.balancing.utility, fout_name, Market.zone_names);
		fout_name = name + "_balancing_demand.csv";
		basic::write_file(Market.balancing.demand_up - Market.balancing.demand_down, fout_name, Market.zone_names);
		fout_name = name + "_balancing_supply.csv";
		basic::write_file(Market.balancing.supply_up - Market.balancing.supply_down, fout_name, Market.zone_names);
	}
}

void power_market::Markets_results_print(market_whole_inform &Power_market_inform){
	// Create a folder to store the file
	std::string dir_name = "csv/output/power_market";
	std::filesystem::create_directories(dir_name);
	dir_name += "/";

	local::Market_results_print(Power_market_inform.International_Market, dir_name + "IMO", 0);
	local::Market_results_print(Power_market_inform.TSO_Market, dir_name + "TSO");
}

void power_market::Simplified_network_print(market_whole_inform &Power_market_inform){
	// Create a folder to store the file
	std::string dir_name = "csv/processed/power_market";
	std::filesystem::create_directories(dir_name);
	dir_name += "/";

	local::Simplified_network_print(Power_market_inform.TSO_Market, dir_name + "TSO");
}
