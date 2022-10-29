// Source file for the main procedure of the power market clearing
#include "power_network/power_network.h"
#include "power_market/power_market.h"
#include "src/spatial_field/spatial_field.h"

namespace{
	struct process_bool{
		bool default_flag;
		bool estimation_flag;
		bool simulation_flag;
		bool DSO_filter_flag;
		bool control_reserve_flag;

		void process_default_get(){
			std::cout << "Default procedure?        Yes: 1 / No: 0 | ";
			std::cin >> this->default_flag;
			std::cout << "\n";
		}

		void process_bool_set(){
			this->estimation_flag = 0;
			this->simulation_flag = 1;
			this->DSO_filter_flag = 0;
			this->control_reserve_flag = 1;
		}

		void process_bool_input(){
			std::cout << "Estimate spatial fields?  Yes: 1 / No: 0 | ";
			std::cin >> this->estimation_flag;
			std::cout << "\n";

			std::cout << "Simulate operation?       Yes: 1 / No: 0 | ";
			std::cin >> this->simulation_flag;
			std::cout << "\n";

			if(this->simulation_flag == 1){
				std::cout << "DSOs filter bids?         Yes: 1 / No: 0 | ";
				std::cin >> this->DSO_filter_flag;
				std::cout << "\n";

				std::cout << "Control reserve?          Yes: 1 / No: 0 | ";
				std::cin >> this->control_reserve_flag;
				std::cout << "\n";
			}
		}
	};
}

int main(){
	// Set booleans for the process
	process_bool process_par;
	process_par.process_default_get();
	if(!process_par.default_flag){
		process_par.process_bool_input();
	}
	else{
		process_par.process_bool_set();
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
		spatial_field::demand_imbalance_estimation(Power_network_inform, Power_market_inform.International_Market);
		spatial_field::wind_on_cf_estimation(Power_network_inform);
		spatial_field::solar_radiation_estimation(Power_network_inform);
	}

	// Power market processes
	if(process_par.simulation_flag){
		power_market::power_market_process_set(Power_network_inform, Power_market_inform, process_par.DSO_filter_flag, process_par.control_reserve_flag);
		power_market::power_market_process_update(Power_network_inform, Power_market_inform, process_par.DSO_filter_flag, process_par.control_reserve_flag);

		// Output results
		power_market::Markets_results_print(Power_market_inform);
	}
}
//	std::cin.get();
