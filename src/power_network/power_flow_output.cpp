// Source file for printing results from simulation and settlement of the results of power flow analysis
#include <filesystem>
#include <iostream>
#include "src/basic/rw_csv.h"
#include "power_flow_analysis.h"

namespace local{
	void power_flow_transmission_results_print(power_market::market_inform &Market, std::string name, configuration::process_config &process_par){
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
		basic::write_file(Market.power_flow.voltage_abs.middleRows(process_par.time_boundary[0], process_par.time_boundary[1]), fout_name, Market.zone_names);
		fout_name = name + "_actual_AC_voltage_angle.csv";
		basic::write_file(Market.power_flow.voltage_arg.middleRows(process_par.time_boundary[0], process_par.time_boundary[1]), fout_name, Market.zone_names);
		fout_name = name + "_actual_AC_P_node.csv";
		basic::write_file(Market.power_flow.P_node.middleRows(process_par.time_boundary[0], process_par.time_boundary[1]), fout_name, Market.zone_names);
		fout_name = name + "_actual_AC_Q_node.csv";
		basic::write_file(Market.power_flow.Q_node.middleRows(process_par.time_boundary[0], process_par.time_boundary[1]), fout_name, Market.zone_names);
		fout_name = name + "_actual_AC_current_magnitude.csv";
		basic::write_file(Market.power_flow.current_abs.middleRows(process_par.time_boundary[0], process_par.time_boundary[1]), fout_name, edge_names);
		fout_name = name + "_actual_AC_current_angle.csv";
		basic::write_file(Market.power_flow.current_arg.middleRows(process_par.time_boundary[0], process_par.time_boundary[1]), fout_name, edge_names);
	}

	void power_flow_whole_results_print(power_network::network_inform &Power_network_inform, std::string name, configuration::process_config &process_par){
		std::string fout_name;

		// Set bus names
		int node_num = Power_network_inform.nodes.bidding_zone.size();
		int point_num = Power_network_inform.points.bidding_zone.size();
		int bus_num = node_num + point_num;
		std::vector <std::string> bus_names;
		bus_names.reserve(bus_num);
		for(int bus_iter = 0; bus_iter < bus_num; ++ bus_iter){
			std::string bus_name = "Bus_" + std::to_string(bus_iter);
			bus_names.push_back(bus_name);
		}

		fout_name = name + "_actual_AC_voltage_magnitude.csv";
		basic::write_file(Power_network_inform.power_flow.voltage_abs.middleRows(process_par.time_boundary[0], process_par.time_boundary[1]), fout_name, bus_names);
		fout_name = name + "_actual_AC_voltage_angle.csv";
		basic::write_file(Power_network_inform.power_flow.voltage_arg.middleRows(process_par.time_boundary[0], process_par.time_boundary[1]), fout_name, bus_names);
	}
}

void power_network::power_flow_results_print(power_market::market_whole_inform &Power_market_inform, network_inform &Power_network_inform, configuration::process_config &process_par){
	// Create a folder to store the file
	std::string dir_name = "csv/output/power_network";
	std::filesystem::create_directories(dir_name);
	dir_name += "/";

	local::power_flow_transmission_results_print(Power_market_inform.TSO_Market, dir_name + "TSO", process_par);
	local::power_flow_whole_results_print(Power_network_inform, dir_name + "whole", process_par);
}
