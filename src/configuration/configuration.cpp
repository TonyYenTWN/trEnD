// Read and store configuration data for the program and model
#include "configuration.h"

void configuration::process_config_input(process_config &process_par, std::string dir_name){
    // Read process configuration data
    auto fin_process = dir_name + "process.csv";
    auto process_inform = basic::read_config_file(fin_process);

    // Create a map for the data
    process_par.estimation_flag = (bool) process_inform["estimation"];
    process_par.estimation_demand_flag = (bool) process_inform["est_demand"];
    process_par.estimation_wind_flag = (bool) process_inform["est_wind"];
    process_par.estimation_solar_flag = (bool) process_inform["est_solar"];
    process_par.simulation_flag = (bool) process_inform["simulation"];
    process_par.DSO_filter_flag = (bool) process_inform["DSO_filter"];
    process_par.control_reserve_flag = (bool) process_inform["control_reserve"];
    process_par.encourage_redispatch= (bool) process_inform["encourage_redispatch"];
    process_par.total_time = (int) process_inform["total_time"];
    process_par.time_boundary.push_back((int) process_inform["start_time"]);
    process_par.time_boundary.push_back((int) process_inform["duration"]);
}