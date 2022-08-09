// Source file for the main procedure of the power market clearing
#include <iostream>
#include <omp.h>
#include <chrono>
//#include "../basic/LP_gpa.h"
#include "basic/rw_csv.h"
#include "alglib/optimization.h"
#include "power_network/power_network.h"
#include "power_market/power_market.h"

int main(){
	power_network::network_inform Power_network_inform;
	power_market::market_whole_inform Power_market_inform;

	power_market::power_market_process_set(Power_network_inform, Power_market_inform, 0);
	power_market::power_market_process_update(Power_network_inform, Power_market_inform, 0);
}
//	std::cin.get();
