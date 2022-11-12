// Source file for printing results from simulation and settlement of the results of power flow analysis
#include <filesystem>
#include <iostream>
#include "src/basic/rw_csv.h"
#include "power_flow_analysis.h"

namespace local{
	void power_flow_results_print(power_market::market_inform &Market, std::string name){
		std::string fout_name;

		// Set edge names
		int edge_num = Market.power_flow.current_abs.cols();
		std::vector <std::string> edge_names;
		edge_names.reserve(edge_num);
		for(int edge_iter = 0; edge_iter < edge_num; ++ edge_iter){
			std::string edge_name = "Edge_" + std::to_string(edge_iter);
			edge_names.push_back(edge_name);
		}

		// Output results
		fout_name = name + "_actual_AC_voltage_magnitude.csv";
		basic::write_file(Market.power_flow.voltage_abs, fout_name, Market.zone_names);
		fout_name = name + "_actual_AC_voltage_angle.csv";
		basic::write_file(Market.power_flow.voltage_arg, fout_name, Market.zone_names);
		fout_name = name + "_actual_AC_P_node.csv";
		basic::write_file(Market.power_flow.P_node, fout_name, Market.zone_names);
		fout_name = name + "_actual_AC_Q_node.csv";
		basic::write_file(Market.power_flow.Q_node, fout_name, Market.zone_names);
		fout_name = name + "_actual_AC_current_magnitude.csv";
		basic::write_file(Market.power_flow.current_abs, fout_name, edge_names);
		fout_name = name + "_actual_AC_current_angle.csv";
		basic::write_file(Market.power_flow.current_arg, fout_name, edge_names);
	}
}

void power_network::power_flow_results_print(power_market::market_whole_inform &Power_market_inform){
	// Create a folder to store the file
	std::string dir_name = "csv/output/power_network";
	std::filesystem::create_directories(dir_name);
	dir_name += "/";

	local::power_flow_results_print(Power_market_inform.TSO_Market, dir_name + "TSO");
}
