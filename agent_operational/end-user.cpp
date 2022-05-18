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
	double smart_appliance_scale;							// Indicates how much of the total demand in the time interval can be shifted around flexibly; assuming the default is constant profile before shifting
	double smart_appliance_flexibility_factor;				// Indicates how flexible the smart appliances are; e.g. 1 / 2 = can concentrate the demand within half of the time interval
	Eigen::VectorXd normalized_default_demand_profile;		// Normalized s.t. base value corresponds to meand demand at the spatial point
	Eigen::VectorXd normalized_PV_profile;					// Normalized with the same base as demand
	Eigen::VectorXi EV_usage_period;						// The time intervals when EV is actually used
	Eigen::VectorXi EV_house_period;						// The time intervals when EV is parked in the house
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
	subscept_tariff << 5, 6, 9, 2, 3;
	sorted_vector sorted_tariff = sort(subscept_tariff);
	
	end_user_operation test_user;
	test_user.normalized_default_demand_profile = Eigen::VectorXd(5);
	test_user.normalized_default_demand_profile << 1, 1.2, 1.5, .8, .6;
}