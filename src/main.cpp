// Source file for the main procedure of the power market clearing
#include <iostream>
#include <omp.h>
//#include <chrono>
//#include "../basic/LP_gpa.h"
#include "basic/rw_csv.h"
#include "alglib/optimization.h"
#include "power_network/power_network.h"
#include "power_market/power_market.h"

int main(){
	// Initialization of power network information
	power_network::network_inform Power_network_inform;
	power_network::power_network_input_process(Power_network_inform, "csv/input/power_network/");
	Power_network_inform.set_line_density();

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
	//power_market::Flow_Based_Market_LP_Set(TSO_Market, TSO_Problem);

	// Initialization of the DSO
	power_market::markets_inform DSO_Markets;
	power_market::DSO_Markets_Set(DSO_Markets, Power_network_inform, Time);
	std::vector <alglib::minlpstate> DSO_Problems (DSO_Markets.size());
	for(int DSO_iter = 0; DSO_iter < DSO_Markets.size(); ++ DSO_iter){
		if(DSO_Markets[DSO_iter].num_zone > 1){
			std::cout << DSO_iter << "\t";
			power_market::Flow_Based_Market_LP_Set(DSO_Markets[DSO_iter], DSO_Problems[DSO_iter]);
		}
		else{
			//std::cout << DSO_iter << "\n\n";
		}
	}

//	alglib::minlpstate* DSO_Problem = new alglib::minlpstate;
//	power_market::Flow_Based_Market_LP_Set(DSO_Markets[0], *DSO_Problem);
//	delete DSO_Problem;
//	std::cout << DSO_Markets[0].num_zone << "\n";
//	DSO_Problem = new alglib::minlpstate;
//	power_market::Flow_Based_Market_LP_Set(DSO_Markets[1], *DSO_Problem);
//	std::cout << DSO_Markets[1].num_zone << "\n";
//	delete DSO_Problem;
//	DSO_Problem = new alglib::minlpstate;
//	power_market::Flow_Based_Market_LP_Set(DSO_Markets[2], *DSO_Problem);
//	std::cout << DSO_Markets[2].num_zone << "\n";
//	delete DSO_Problem;
//	DSO_Problem = new alglib::minlpstate;
//	power_market::Flow_Based_Market_LP_Set(DSO_Markets[3], *DSO_Problem);
//	std::cout << DSO_Markets[3].num_zone << "\n";
//	delete DSO_Problem;
//	DSO_Problem = new alglib::minlpstate;
//	power_market::Flow_Based_Market_LP_Set(DSO_Markets[4], *DSO_Problem);
//	std::cout << DSO_Markets[4].num_zone << "\n";
//	delete DSO_Problem;
//	DSO_Problem = new alglib::minlpstate;
//	power_market::Flow_Based_Market_LP_Set(DSO_Markets[5], *DSO_Problem);
//	std::cout << DSO_Markets[5].num_zone << "\n";
//	delete DSO_Problem;
//	alglib::minlpstate* DSO_Problem_2 = new alglib::minlpstate;
//	power_market::Flow_Based_Market_LP_Set(DSO_Markets[6], *DSO_Problem_2);
//	std::cout << DSO_Markets[6].num_zone << "\n";
//	delete DSO_Problem_2;

//	for(int DSO_iter = 0; DSO_iter < 5; ++ DSO_iter){
//		alglib::minlpstate DSO_Problem;
//		power_market::Flow_Based_Market_LP_Set(DSO_Markets[0], DSO_Problem);
//	}
//	auto DSO_Market_temp = DSO_Markets[0];
//	power_market::Flow_Based_Market_LP_Set(DSO_Market_temp);
//	std::thread th0(power_market::Flow_Based_Market_LP_Set, DSO_Market_temp);
//	DSO_Markets[1].Problem = DSO_Markets[0].Problem;
	//for(int DSO_iter = 0; DSO_iter < DSO_Markets.size(); ++ DSO_iter){
//	for(int DSO_iter = 0; DSO_iter < 12; ++ DSO_iter){
//		auto DSO_Market_temp = DSO_Markets[DSO_iter];
//		if(DSO_Markets[DSO_iter].num_zone > 1){
//			//power_market::Flow_Based_Market_LP_Set(DSO_Markets[DSO_iter]);
//			//power_market::Flow_Based_Market_LP_Set(*DSO_Market_temp);
//		}
//	}
	//std::cout << DSO_Markets[6].num_zone << "\n\n";

	// Re-initialization of submitted bids
	std::string fin_point_demand = "csv/processed/spatial_field/nominal_mean_demand_field_10km_annual_mean.csv";
	power_market::Submitted_bid_calculation(0, DSO_Markets, TSO_Market, International_Market, Power_network_inform, fin_point_demand);

	// Ideal market clearing in IMO

	// Bid-filtering in DSOs

	// Re-dispatch + tertiary control reserve in TSO
	//power_market::Flow_Based_Market_Optimization(0, TSO_Market, TSO_Problem);
}
