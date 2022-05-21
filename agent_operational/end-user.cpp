// Main Source File
#include <iostream>
#include <chrono>
#include <omp.h>
//#include <boost/math/distributions/normal.hpp>
#include "../basic/Basic_Definitions.h"

struct end_user_operation{
	// Input parameters
	int point_ID;
	bool dynamic_tariff;
	bool PV_BESS;
	bool EV_self_charging;
	bool smart_appliance;
	bool reverse_flow;										// Whether the end-user can inject power flow back to grid; false when active_flow is false
	bool active_flex;										// Whether the end-user can provide flexibility to the aggregator; false when dynamic_tariff is false, or when end-user does not have PV + BESS, EV, or smart appliance
	double PV_scale;										// Unless mentioned otherwise, all the scaling factor = capacity / normalized default base value
	double BESS_energy_scale;								// kWh per person
	double BESS_capacity_scale;
	double BESS_efficiency;
	double EV_energy_scale;
	double EV_capacity_scale;
	double EV_efficiency;
	double smart_appliance_scale;							// Indicates how much ratio of the total demand in the time interval can be shifted around flexibly; assuming the default is constant profile before shifting
	double smart_appliance_flexibility_factor;				// Indicates how flexible the smart appliances are; e.g. 1 / 2 = can concentrate the demand within half of the time interval
	Eigen::VectorXd normalized_default_demand_profile;		// Normalized to nominal value (kWh per hour per person)
	Eigen::VectorXd normalized_default_PV_profile;			// Normalized with the same base as demand
	Eigen::VectorXi EV_usage_default_period;				// The time intervals when EV is actually used
	Eigen::VectorXi EV_house_default_period;				// The time intervals when EV is parked in the house
	
	// Output variables
	Eigen::VectorXd normalized_scheduled_residual_demand_inflex_profile;
	Eigen::VectorXd normalized_scheduled_BESS_capacity_profile;
	Eigen::VectorXd normalized_scheduled_BESS_energy_profile;
	Eigen::VectorXd normalized_scheduled_EV_capcity_profile;
	Eigen::VectorXd normalized_scheduled_EV_energy_profile;
	Eigen::VectorXd normalized_scheduled_smart_appliance_profile;
	Eigen::VectorXd normalized_scheduled_pos_cr_profile;
	Eigen::VectorXd normalized_scheduled_neg_cr_profile;
};

struct sorted_vector{
	Eigen::VectorXi id;
	Eigen::VectorXd value;
};

sorted_vector sort(Eigen::VectorXd original){
	// Sort of vector
 	std::vector<int> item_seq(original.size());
 	std::iota(item_seq.begin(), item_seq.end(), 0); 
 	sort(item_seq.begin(), item_seq.end(), [&](int i,int j){return original(i) < original(j);});
 	
 	// Output of vector
 	sorted_vector result;
 	result.id = Eigen::VectorXi(original.size());
 	result.value = Eigen::VectorXd(original.size());
 	for(int item_ID = 0; item_ID < original.size(); ++ item_ID){
 		result.id(item_ID) = item_seq[item_ID];
 		result.value(item_ID) = original(item_seq[item_ID]);
	}
	
	return(result);
}

int main(){
	// Test case
	Eigen::VectorXd subscept_tariff(10);
	subscept_tariff << 5, 6, 9, 2, 3, 1, 6, 9, 10, 2;
	sorted_vector sorted_tariff = sort(subscept_tariff);
	
	end_user_operation test_user;
	test_user.normalized_default_demand_profile = Eigen::VectorXd(subscept_tariff.size());
	test_user.normalized_default_demand_profile << 1, 1.2, 1.5, .8, .6, 2, 3, .1, .5, 1.2;
	test_user.smart_appliance_scale = .2;
	test_user.smart_appliance_flexibility_factor = .5;
	
	// Smart appliance test
	double test_total_flex_demand_energy = test_user.smart_appliance_scale * test_user.normalized_default_demand_profile.sum();
	double test_max_flex_demand_capacity = test_total_flex_demand_energy / test_user.normalized_default_demand_profile.size() / test_user.smart_appliance_flexibility_factor;
	int flex_demand_duration = int (test_user.smart_appliance_flexibility_factor * test_user.normalized_default_demand_profile.size());
	if(double (test_user.smart_appliance_flexibility_factor * test_user.normalized_default_demand_profile.size() - flex_demand_duration) != 0){
		flex_demand_duration += 1;
	}
	test_user.normalized_scheduled_smart_appliance_profile = Eigen::VectorXd::Zero(test_user.normalized_default_demand_profile.size());
	for(int tock = 0; tock < flex_demand_duration - 1; ++ tock){
		test_user.normalized_scheduled_smart_appliance_profile(sorted_tariff.id(tock)) = test_max_flex_demand_capacity;
	}
	test_user.normalized_scheduled_smart_appliance_profile(sorted_tariff.id(flex_demand_duration - 1)) = test_total_flex_demand_energy - (flex_demand_duration - 1) * test_max_flex_demand_capacity;
	
	std::cout << test_user.normalized_scheduled_smart_appliance_profile.transpose() << std::endl;
	std::cout << test_user.normalized_scheduled_smart_appliance_profile.transpose() + (1 - test_user.smart_appliance_scale) * test_user.normalized_default_demand_profile.transpose() << std::endl;
	
	// BESS test, naive
	test_user.BESS_energy_scale = 4;
	test_user.BESS_capacity_scale = 1;
	test_user.BESS_efficiency = .95;
	test_user.normalized_scheduled_BESS_energy_profile = Eigen::VectorXd::Zero(subscept_tariff.size());
	// Need to determine best end storage level (given start storage level) in the future
	test_user.normalized_scheduled_BESS_energy_profile(0) = test_user.BESS_energy_scale;
	test_user.normalized_scheduled_BESS_energy_profile(subscept_tariff.size() - 1) = test_user.BESS_energy_scale;
	int BESS_ch_duration = int ((test_user.BESS_energy_scale - (std::max(test_user.normalized_scheduled_BESS_energy_profile(0) - test_user.normalized_scheduled_BESS_energy_profile(subscept_tariff.size() - 1), 0.))) / test_user.BESS_capacity_scale);
	int BESS_dc_duration = int ((test_user.BESS_energy_scale - (std::max(test_user.normalized_scheduled_BESS_energy_profile(subscept_tariff.size() - 1) - test_user.normalized_scheduled_BESS_energy_profile(0), 0.))) / test_user.BESS_capacity_scale);
	//std::cout << BESS_ch_duration << " " << BESS_dc_duration << std::endl;
	
	for(int tock = 0; tock < BESS_ch_duration; ++ tock){
		test_user.normalized_scheduled_BESS_energy_profile(sorted_tariff.id(tock)) = -test_user.BESS_capacity_scale / test_user.BESS_efficiency;
	}
	for(int tock = 0; tock < BESS_dc_duration; ++ tock){
		test_user.normalized_scheduled_BESS_energy_profile(sorted_tariff.id(sorted_tariff.id.size() - tock)) = test_user.BESS_capacity_scale * test_user.BESS_efficiency;
	}
	std::cout << test_user.normalized_scheduled_BESS_energy_profile.transpose() << std::endl;
	
	// EV test
	int EV_activity_duration = 0;
	test_user.EV_usage_default_period = Eigen::VectorXi::Zero(subscept_tariff.size());
	test_user.EV_usage_default_period(3) = 1;
	test_user.EV_usage_default_period(6) = 1;
	test_user.EV_house_default_period = Eigen::VectorXi::Zero(subscept_tariff.size());
	test_user.EV_house_default_period.head(3) << 1, 1, 1;
	test_user.EV_house_default_period.tail(3) << 1, 1, 1;
	if(test_user.EV_house_default_period(0) == 1){
		int tock = 0;
		while(test_user.EV_house_default_period(tock) == 1){
			tock += 1;
		}
		
		// Then use the BESS charge / discharge function
		// In fact, all flexibility resources mentioned above should be able to capture in the same function
	}
}