// Main Source File
#include <iostream>
#include <chrono>
#include <omp.h>
#include <boost/math/distributions/normal.hpp>
#include "../basic/Basic_Definitions.h"

struct end_user_operation{
	// Input parameters
	int point_ID;
	bool dynamic_tariff;
	bool PV_BESS;
	bool EV_self_charging;
	bool smart_appliance;
	bool active_flex;
	double PV_scale;										// All the scaling factor = capacity / normalized default base value
	double BESS_scale;
	double EV_scale;
	int EV_full_charge_duration;
	double smart_appliance_scale;							// Indicates how much ratio of the total demand in the time interval can be shifted around flexibly; assuming the default is constant profile before shifting
	double smart_appliance_flexibility_factor;				// Indicates how flexible the smart appliances are; e.g. 1 / 2 = can concentrate the demand within half of the time interval
	Eigen::VectorXd normalized_default_demand_profile;		// Normalized to nominal value (kWh per hour per person)
	Eigen::VectorXd normalized_default_PV_profile;			// Normalized with the same base as demand
	Eigen::VectorXi EV_usage_default_period;				// The time intervals when EV is actually used
	Eigen::VectorXi EV_house_default_period;				// The time intervals when EV is parked in the house
	
	// Output variables
	Eigen::VectorXd normalized_scheduled_BESS_profile;
	Eigen::VectorXd normalized_scheduled_EV_profile;
	Eigen::VectorXd normalized_scheduled_smart_appliance_profile;
	Eigen::VectorXd normalized_scheduled_pos_flex_profile;
	Eigen::VectorXd normalized_scheduled_neg_flex_profile;
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
	Eigen::VectorXd subscept_tariff(5);
	subscept_tariff << 5, 6, 9, 2, 3;
	sorted_vector sorted_tariff = sort(subscept_tariff);
	
	end_user_operation test_user;
	test_user.normalized_default_demand_profile = Eigen::VectorXd(5);
	test_user.normalized_default_demand_profile << 1, 1.2, 1.5, .8, .6;
	test_user.smart_appliance_scale = .2;
	test_user.smart_appliance_flexibility_factor = .5;
	
	// Flexible demand test
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
	
//	std::cout << test_total_flex_demand_energy << " " << test_max_flex_demand_capacity << " " << flex_demand_duration << std::endl;
//	std::cout << subscept_tariff.transpose() << std::endl;
	std::cout << test_user.normalized_scheduled_smart_appliance_profile.transpose() << std::endl;
	std::cout << test_user.normalized_scheduled_smart_appliance_profile.transpose() + (1 - test_user.smart_appliance_scale) * test_user.normalized_default_demand_profile.transpose() << std::endl;
}