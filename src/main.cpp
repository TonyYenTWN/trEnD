// Source file for the main procedure of the power market clearing
#include "power_network/power_network.h"
#include "power_market/power_market.h"
#include "src/spatial_field/spatial_field.h"

namespace{
	struct process_bool{
		bool inference_flag;
		bool DSO_filter_flag;

		void process_bool_set(){
			this->inference_flag = 0;
			this->DSO_filter_flag = 1;
		}

		void process_bool_input(){
			std::cout << "Inference spatial fields? Yes: 1 / No: 0 | ";
			std::cin >> this->inference_flag;
			std::cout << "\n";

			std::cout << "DSOs filter bids? Yes: 1 / No: 0 | ";
			std::cin >> this->DSO_filter_flag;
			std::cout << "\n";
		}
	};
}

int main(){
	// Initialization of power network information
	power_network::network_inform Power_network_inform;
	power_network::power_network_input_process(Power_network_inform, "csv/input/power_network/");

	// Set booleans for the process
	process_bool process_par;
	//process_par.process_bool_input();
	process_par.process_bool_set();

	// Spatial fields inference
	if(process_par.inference_flag){
		spatial_field::spatial_field_inference(Power_network_inform);
	}

	// Power market processes
	power_market::market_whole_inform Power_market_inform;
	power_market::power_market_process_set(Power_network_inform, Power_market_inform, process_par.DSO_filter_flag);
	power_market::power_market_process_update(Power_network_inform, Power_market_inform, process_par.DSO_filter_flag);
}
//	std::cin.get();
