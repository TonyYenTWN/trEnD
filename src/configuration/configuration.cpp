// Read and store configuration data for the program and model
#include "configuration.h"

void configuration::process_config_input(process_config &process_par, std::string dir_name){
    // Read process configuration data
    auto fin_process = dir_name + "process.csv";
    auto process_inform = basic::read_config_file(fin_process);

    // Create a map for the data
    process_par.estimation_flag = (bool) stod(process_inform["estimation"][0]);
    process_par.estimation_demand_flag = (bool) stod(process_inform["est_demand"][0]);
    process_par.estimation_wind_flag = (bool) stod(process_inform["est_wind"][0]);
    process_par.estimation_solar_flag = (bool) stod(process_inform["est_solar"][0]);
    process_par.simulation_flag = (bool) stod(process_inform["simulation"][0]);
    process_par.DSO_filter_flag = (bool) stod(process_inform["DSO_filter"][0]);
    process_par.control_reserve_flag = (bool) stod(process_inform["control_reserve"][0]);
    process_par.encourage_redispatch= (bool) stod(process_inform["encourage_redispatch"][0]);
    process_par.total_time = (int) stod(process_inform["total_time"][0]);
    process_par.time_boundary.push_back((int) stod(process_inform["start_time"][0]));
    process_par.time_boundary.push_back((int) stod(process_inform["duration"][0]));
    process_par.power_flow = (bool) stod(process_inform["power_flow"][0]);
    process_par.contingency_flag = (bool) stod(process_inform["contingency"][0]);
    process_par.contingency_sampling = (bool) stod(process_inform["contingency_sampling"][0]);
    process_par.folder_name = process_inform["folder_name"][0];
}
