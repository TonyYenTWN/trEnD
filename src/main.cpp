// Source file for the main procedure of the power market clearing
#include <filesystem>
#include "src/agent/agent_func.h"
#include "src/power_market/power_market.h"
#include "src/power_network/contingency_analysis.h"
#include "src/power_network/power_flow_analysis.h"
#include "src/power_network/power_network.h"
#include "src/spatial_field/spatial_field.h"

int main(int argc, char** argv){
	// Set booleans for the process
	configuration::process_config process_par;
	if(argc == 1){
        std::cout << "Folder name?    | ";
        std::cin >> process_par.folder_name;
        std::cout << "\n";
	}
	else{
        process_par.folder_name = argv[1];

        if(argc == 3){
            process_par.contingency_sample_number = argv[2];
        }
        else{
            process_par.contingency_sample_number = 10;
        }
	}

    process_config_input(process_par, "csv/" + process_par.folder_name + "/configuration/");

	if(process_par.estimation_flag + process_par.simulation_flag + process_par.contingency_flag == 0){
		std::cout << "No process selected. Exit program...";
		return 0;
	}

	// Set folder and file name for log messages
	std::string output_dir_name = "csv/" + process_par.folder_name + "/output";
	std::string output_log_name = output_dir_name + "/log.txt";
	std::filesystem::create_directories(output_dir_name);
	std::freopen(output_log_name.c_str() , "w", stdout);
	process_par.process_bool_output();

	// Initialization of power network information
	power_network::network_inform Power_network_inform;
	power_network::power_network_input_process(Power_network_inform, process_par.hydro_factor, "csv/" + process_par.folder_name + "/input/power_network/");

	// Set bidding prices and default (residual) demand time series
	power_market::market_whole_inform Power_market_inform;
	power_market::parameters::bidded_price(Power_market_inform.price_map);
	power_market::default_demand_set(Power_network_inform, Power_market_inform, process_par);

	// Spatial fields estimation
	if(process_par.estimation_flag){
        // Create a folder to store the file
		std::filesystem::create_directories("csv/" + process_par.folder_name + "/processed/spatial_field");

		if(process_par.estimation_demand_flag){
            std::filesystem::create_directories("csv/" + process_par.folder_name + "/processed/spatial_field/demand");
            std::filesystem::create_directories("csv/" + process_par.folder_name + "/processed/spatial_field/imbalance");
            spatial_field::demand_imbalance_estimation(Power_network_inform, Power_market_inform.International_Market, process_par);
		}

		if(process_par.estimation_wind_flag){
            std::filesystem::create_directories("csv/" + process_par.folder_name + "/processed/spatial_field/wind");
            spatial_field::wind_on_cf_estimation(Power_network_inform, process_par);
		}

		if(process_par.estimation_solar_flag){
            std::filesystem::create_directories("csv/" + process_par.folder_name + "/processed/spatial_field/solar");
            spatial_field::solar_radiation_estimation(Power_network_inform, process_par);
		}
	}

	// Power market processes
	if(process_par.simulation_flag){
        power_market::power_market_process_set(Power_network_inform, Power_market_inform, process_par);
		power_market::power_market_process_update(Power_network_inform, Power_market_inform, process_par);

		// Output results
		power_market::Markets_results_print(Power_market_inform, process_par);
		power_network::power_flow_results_print(Power_market_inform, Power_network_inform, process_par);
		agent::agents_results_print(Power_market_inform, Power_network_inform, process_par);
	}

	// Contigency analysis
	if(process_par.contingency_flag){
        //  Read and store flex_stat data to TSO
        power_network::flex_stat_input(Power_market_inform, Power_network_inform, process_par);

        // Initialization of contingency analysis object
        power_network::contingency_analysis_struct contingency_analysis;

        // Sampling of contingencies
        power_network::contingency_analysis_set(contingency_analysis, Power_market_inform, process_par);
        power_network::contigency_sampling(contingency_analysis, 0, process_par); // default samples = 1E5
        power_network::contingency_analysis_solve(contingency_analysis, Power_market_inform, Power_network_inform, process_par);
        power_network::contingency_analysis_print(contingency_analysis, Power_market_inform, process_par);
	}

	// Close log file
	std::fclose(stdout);
}
