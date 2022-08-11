// Source file for the main procedure of the power market clearing
#include "power_network/power_network.h"
#include "power_market/power_market.h"
#include "src/spatial_field/spatial_field.h"

int main(){
	// Initialization of power network information
	power_network::network_inform Power_network_inform;
	power_network::power_network_input_process(Power_network_inform, "csv/input/power_network/");

	spatial_field::nominal_demand_inference(Power_network_inform);

	power_market::market_whole_inform Power_market_inform;
	//power_market::power_market_process_set(Power_network_inform, Power_market_inform, 0);
	//power_market::power_market_process_update(Power_network_inform, Power_market_inform, 0);
}
//	std::cin.get();
