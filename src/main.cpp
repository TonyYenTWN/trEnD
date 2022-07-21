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
	auto start = std::chrono::high_resolution_clock::now();
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast <std::chrono::microseconds> (stop - start);
	std::cout << "Current time: " << duration.count() << " microseconds" << "\n\n";

	start = std::chrono::high_resolution_clock::now();

	// Initialization of power network information
	power_network::network_inform Power_network_inform;
	power_network::power_network_input_process(Power_network_inform, "csv/input/power_network/");

	// Initialization of the IMO
	int Time = 8760;
	std::string fin_name_moc = "csv/input/power_market/merit_order_curve_q_assimilated_2021.csv";
	std::string fin_name_demand = "csv/input/power_market/residual_load_default_forecast_2021.csv";
	power_market::market_inform International_Market;
	power_market::International_Market_Set(International_Market, Time, fin_name_moc, fin_name_demand);

	// Initialization of the TSO
	power_market::market_inform TSO_Market;
	power_market::TSO_Market_Set(TSO_Market, Power_network_inform, Time);
	alglib::minlpstate TSO_Problem;
	power_market::Flow_Based_Market_LP_Set(TSO_Market, TSO_Problem);

	// Initialization of the DSO
	power_market::markets_inform DSO_Markets;
	power_market::DSO_Markets_Set(DSO_Markets, Power_network_inform, Time);
	std::vector <alglib::minlpstate> DSO_Problems (DSO_Markets.size());
	for(int DSO_iter = 0; DSO_iter < DSO_Markets.size(); ++ DSO_iter){
		power_market::Flow_Based_Market_LP_Set(DSO_Markets[DSO_iter], DSO_Problems[DSO_iter]);
	}

	// Re-initialization of submitted bids in DSOs and TSOs
	std::string fin_point_demand = "csv/processed/spatial_field/nominal_mean_demand_field_10km_annual_mean.csv";
	power_market::Submitted_bid_calculation(DSO_Markets, TSO_Market, International_Market, Power_network_inform, fin_point_demand);

	// Ideal market clearing in IMO
	power_market::International_Market_Optimization(0, International_Market, 0);

	stop = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast <std::chrono::microseconds> (stop - start);
	std::cout << "Set time: " << duration.count() << " microseconds" << "\n\n";

	// Bid-filtering in DSOs
	//#pragma omp parallel
	{
		//#pragma omp for schedule(dynamic)
		//for(int DSO_iter = 0; DSO_iter < 5; ++ DSO_iter){
		for(int DSO_iter = 0; DSO_iter < DSO_Markets.size(); ++ DSO_iter){
			start = std::chrono::high_resolution_clock::now();
			std::cout << DSO_iter << "-th DSO \n";

			power_market::Source_Node_Set(DSO_Markets[DSO_iter], Power_network_inform.DSO_cluster[DSO_iter]);
			power_market::Flow_Based_Market_Optimization(DSO_Markets[DSO_iter], DSO_Problems[DSO_iter]);
			power_market::DSO_Market_Results_Get(DSO_Markets[DSO_iter], DSO_Problems[DSO_iter], Power_network_inform.DSO_cluster[DSO_iter], 0);

			power_market::Sink_Node_Set(DSO_Markets[DSO_iter], Power_network_inform.DSO_cluster[DSO_iter]);
			power_market::Flow_Based_Market_Optimization(DSO_Markets[DSO_iter], DSO_Problems[DSO_iter]);
			power_market::DSO_Market_Results_Get(DSO_Markets[DSO_iter], DSO_Problems[DSO_iter], Power_network_inform.DSO_cluster[DSO_iter], 1);

			stop = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast <std::chrono::microseconds> (stop - start);
			std::cout << duration.count() << " microseconds" << "\n";
			std::cout << "--------------------------------------------------------------------------------------------------------------\n\n";
		}
	}

	// Re-dispatch + tertiary control reserve in TSO
	start = std::chrono::high_resolution_clock::now();

	std::cout << "TSO: \n";
	power_market::Flow_Based_Market_Optimization(TSO_Market, TSO_Problem);

	stop = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast <std::chrono::microseconds> (stop - start);
	std::cout << "TSO optimization time: " << duration.count() << " microseconds" << "\n\n";

	power_market::TSO_Market_Results_Get(0, TSO_Market, TSO_Problem);
//	//std::cin.get();
}
