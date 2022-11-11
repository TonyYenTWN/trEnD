// Header file for power flow analysis
//#include "src/configuration/configuration.h"
#include "src/power_market/power_market.h"
#include "src/spatial_field/geostat.h"
#include "power_network.h"

namespace power_network{
	// Function for power flow analysis
	void HELM_Transmission_Set(network_inform&, power_market::market_whole_inform&);
	void HELM_Node_Update(int, network_inform&, power_market::market_whole_inform&);
	void HELM_Transmission_Solve(int, network_inform&, power_market::market_whole_inform&);
	//void HELM_Set(network_inform&, power_market::market_whole_inform&);
	//void HELM_Node_Update(int, network_inform&, power_market::market_whole_inform&);
	//void HELM_Solve(int, network_inform&);
}
