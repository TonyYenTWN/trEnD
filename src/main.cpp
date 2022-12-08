// Source file for the main procedure of the power market clearing
#include <filesystem>
#include "src/agent/agent_func.h"
//#include "src/configuration/configuration.h"
#include "src/power_market/power_market.h"
#include "src/power_network/power_flow_analysis.h"
#include "src/power_network/power_network.h"
#include "src/spatial_field/spatial_field.h"

int main(){
	// Set booleans for the process
	configuration::process_config process_par;
	process_par.process_default_get();
	if(!process_par.default_flag){
		process_par.process_bool_input();
	}
	else{
		process_par.process_bool_set();
	}

	if(!process_par.estimation_flag && !process_par.simulation_flag){
		std::cout << "No process selected. Exit program...";
		return 0;
	}

	// Initialization of power network information
	power_network::network_inform Power_network_inform;
	power_network::power_network_input_process(Power_network_inform, "csv/input/power_network/");

	// Set bidding prices and default (residual) demand time series
	power_market::market_whole_inform Power_market_inform;
	power_market::parameters::bidded_price(Power_market_inform.price_map);
	power_market::default_demand_set(Power_network_inform, Power_market_inform);

	// Spatial fields estimation
	if(process_par.estimation_flag){
		// Create a folder to store the file
		std::filesystem::create_directories("csv/processed/spatial_field");

		spatial_field::demand_imbalance_estimation(Power_network_inform, Power_market_inform.International_Market, process_par);
		spatial_field::wind_on_cf_estimation(Power_network_inform, process_par);
		//spatial_field::solar_radiation_estimation(Power_network_inform, process_par);
	}

	// Power market processes
	if(process_par.simulation_flag){
		power_market::power_market_process_set(Power_network_inform, Power_market_inform, process_par);
		power_market::power_market_process_update(Power_network_inform, Power_market_inform, process_par);

		// Output results
//		power_market::Markets_results_print(Power_market_inform, process_par);
//		power_network::power_flow_results_print(Power_market_inform, Power_network_inform, process_par);
//		agent::agents_results_print(Power_market_inform, Power_network_inform);
	}
}
//	std::cin.get();
