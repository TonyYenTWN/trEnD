// Source file for end-user operational strategy
#include <iostream>
//#include <chrono>
#include <omp.h>
#include "../basic/Basic_Definitions.h"

struct end_user_decision{
	bool dynamic_tariff;
	bool smart_appliance;
	bool PV_BESS;
	bool EV_self_charging;
	bool reverse_flow;										// Whether the end-user can inject power flow back to grid; false when active_flex is false
	bool active_flex;										// Whether the end-user can provide flexibility to the aggregator; false when dynamic_tariff is false, or when end-user does not have PV + BESS, EV, or smart appliance
};

struct smart_appliance_inform{
	// Input parameters
	double scale;											// Indicates how much ratio of the total demand in the time interval can be shifted around flexibly; assuming the default is constant profile before shifting
	double flexibility_factor;								// Indicates how flexible the smart appliances are; e.g. 1 / 2 = can concentrate the demand within half of the time interval

	// Output variables
	Eigen::VectorXd normalized_scheduled_profile;
};

struct storage_inform{
	// Input parameters
	double energy_scale;									// kWh per person
	double capacity_scale;									// kW per person
	double efficiency;
	double soc_ini;
	double soc_final;
	
	// Output variables
	Eigen::VectorXd normalized_scheduled_capacity_profile;
	Eigen::VectorXd normalized_scheduled_soc_profile;
};

struct EV_inform{
	// Input parameters
	double energy_demand;									// kWh per person per hour of usage
	Eigen::VectorXi usage_default_period;					// The time intervals when EV is actually used
	Eigen::VectorXi house_default_period;					// The time intervals when EV is parked in the house
	
	// Mixed Substructure
	storage_inform BESS;
};

struct end_user_operation{
	// Input parameters
	int point_ID;
	end_user_decision decision;
	double PV_scale;										// Unless mentioned otherwise, all the scaling factor = capacity / normalized default base value
	Eigen::VectorXd normalized_default_demand_profile;		// Normalized to nominal value (kWh per hour per person)
	Eigen::VectorXd normalized_default_PV_profile;			// Normalized with the same base as demand
	
	// Output variables
	Eigen::VectorXd normalized_scheduled_residual_demand_inflex_profile;
	Eigen::VectorXd normalized_scheduled_residual_demand_flex_profile;
	Eigen::VectorXd normalized_scheduled_pos_cr_profile;
	Eigen::VectorXd normalized_scheduled_neg_cr_profile;
	
	// Mixed Substructure
	smart_appliance_inform smart_appliance;
	storage_inform BESS;
	EV_inform EV;
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

void smart_appliance_schedule(sorted_vector sorted_tariff, Eigen::VectorXd normalized_default_demand_profile, double residual_unfulfilled_demand, smart_appliance_inform &result){
	// Initialization
	double total_flex_demand_energy = result.scale * normalized_default_demand_profile.sum() + residual_unfulfilled_demand;
	double flex_demand_capacity_max = total_flex_demand_energy / normalized_default_demand_profile.size() / result.flexibility_factor;
	int flex_demand_duration = int(result.flexibility_factor * normalized_default_demand_profile.size());
	if(double(result.flexibility_factor * normalized_default_demand_profile.size() - flex_demand_duration) != 0){
		flex_demand_duration += 1;
	}
	
	// Schedule the demand to low price periods
	result.normalized_scheduled_profile = Eigen::VectorXd::Zero(normalized_default_demand_profile.size());
	for(int tick = 0; tick < flex_demand_duration - 1; ++ tick){
		result.normalized_scheduled_profile(sorted_tariff.id(tick)) = flex_demand_capacity_max;
	}	
	result.normalized_scheduled_profile(sorted_tariff.id(flex_demand_duration - 1)) = total_flex_demand_energy - (flex_demand_duration - 1) * flex_demand_capacity_max;
}

void storage_schedule(sorted_vector sorted_tariff, storage_inform &result){
	// Initialization
	int total_duration = sorted_tariff.value.size();
	int ch_surplus_duration = int(std::max(result.soc_final - result.soc_ini, 0.) / result.capacity_scale);
	int dc_surplus_duration = int(std::max(result.soc_ini - result.soc_final, 0.) / result.capacity_scale);
	result.normalized_scheduled_capacity_profile = Eigen::VectorXd::Zero(total_duration);
	result.normalized_scheduled_soc_profile = result.soc_ini * Eigen::VectorXd::Ones(total_duration);
	// Allowed exchange matrix: row = charge period; col = discharge period
	Eigen::MatrixXi Allowed_Exchange = Eigen::MatrixXi::Ones(total_duration, total_duration);
	
	// Schedule of surplus dc / ch so that the soc change is close to as expected
	if(ch_surplus_duration > 0){
		for(int tick = 0; tick < ch_surplus_duration; ++ tick){
			result.normalized_scheduled_capacity_profile(sorted_tariff.id(tick)) = -result.capacity_scale / result.efficiency;
			for(int tock = sorted_tariff.id(tick); tock < total_duration; ++ tock){
				result.normalized_scheduled_soc_profile(tock) += result.capacity_scale;
			}
			Allowed_Exchange.row(sorted_tariff.id(tick)) = Eigen::VectorXi::Zero(total_duration);
			Allowed_Exchange.col(sorted_tariff.id(tick)) = Eigen::VectorXi::Zero(total_duration);
		}
	}
	else{
		for(int tick = 0; tick < dc_surplus_duration; ++ tick){
			result.normalized_scheduled_capacity_profile(sorted_tariff.id(total_duration - tick - 1)) = result.capacity_scale * result.efficiency;
			for(int tock = sorted_tariff.id(total_duration - tick - 1); tock < total_duration; ++ tock){
				result.normalized_scheduled_soc_profile(tock) -= result.capacity_scale;
			}
			Allowed_Exchange.row(sorted_tariff.id(total_duration - tick - 1)) = Eigen::VectorXi::Zero(total_duration);
			Allowed_Exchange.col(sorted_tariff.id(total_duration - tick - 1)) = Eigen::VectorXi::Zero(total_duration);
		}
	}
	
	// Check profitability of exchange
	for(int tick = 0; tick < total_duration; ++ tick){
		for(int tock = 0; tock < total_duration; ++ tock){
			if(Allowed_Exchange(sorted_tariff.id(tick), sorted_tariff.id(tock)) == 1){
				if(sorted_tariff.value(tick) >= pow(result.efficiency, 2) * sorted_tariff.value(tock)){
					Allowed_Exchange(sorted_tariff.id(tick), sorted_tariff.id(tock)) = 0;
				}
			}
		}
	}
	
	// Pair 2 time intervals for maximum possible price arbitrage
	double maximum_price_diff;
	Eigen::Vector2i maximum_price_diff_ID;
	Eigen::MatrixXi Allowed_Exchange_temp = Allowed_Exchange;
	Eigen::MatrixXd soc_profile_temp;
	int count = 0;
	while(Allowed_Exchange_temp.sum() > 0){
	//while(count < 10){
		count += 1;
		maximum_price_diff = 0;
		Allowed_Exchange_temp = Allowed_Exchange;
		
		// Check if the excahnge pair is feasible and if yes, most profitable up to point
		for(int tick = 0; tick < total_duration; ++ tick){
			for(int tock = 0; tock < total_duration; ++ tock){
				if(Allowed_Exchange(sorted_tariff.id(tick), sorted_tariff.id(tock)) == 1){
					soc_profile_temp = result.normalized_scheduled_soc_profile;
					if(sorted_tariff.id(tick) < sorted_tariff.id(tock)){
						for(int tuck = sorted_tariff.id(tick); tuck < sorted_tariff.id(tock); ++ tuck){
							soc_profile_temp(tuck) += result.capacity_scale;
							if(soc_profile_temp(tuck) > result.energy_scale){
								Allowed_Exchange_temp(sorted_tariff.id(tick), sorted_tariff.id(tock)) = 0;
								break;
							}
						}						
					}
					else{
						for(int tuck = sorted_tariff.id(tock); tuck < sorted_tariff.id(tick); ++ tuck){
							soc_profile_temp(tuck) -= result.capacity_scale;
							if(soc_profile_temp(tuck) < 0){
								Allowed_Exchange_temp(sorted_tariff.id(tick), sorted_tariff.id(tock)) = 0;
								break;
							}							
						}						
					}
					
					if(Allowed_Exchange_temp(sorted_tariff.id(tick), sorted_tariff.id(tock)) == 1){
						if(sorted_tariff.value(tock) * result.efficiency - sorted_tariff.value(tick) / result.efficiency > maximum_price_diff){
							maximum_price_diff = sorted_tariff.value(tock) * result.efficiency - sorted_tariff.value(tick) / result.efficiency;
							maximum_price_diff_ID << sorted_tariff.id(tick), sorted_tariff.id(tock);
						}
					}
				}
			}
		}
		
		// Update if feasible price arbitrage between pairs exists
		if(Allowed_Exchange_temp.sum() > 0){
			result.normalized_scheduled_capacity_profile(maximum_price_diff_ID(0)) = -result.capacity_scale / result.efficiency;
			result.normalized_scheduled_capacity_profile(maximum_price_diff_ID(1)) = result.capacity_scale * result.efficiency;
			
			if(maximum_price_diff_ID(0) < maximum_price_diff_ID(1)){
				for(int tuck = maximum_price_diff_ID(0); tuck < maximum_price_diff_ID(1); ++ tuck){
					result.normalized_scheduled_soc_profile(tuck) += result.capacity_scale;
				}			
			}
			else{
				for(int tuck = maximum_price_diff_ID(1); tuck < maximum_price_diff_ID(0); ++ tuck){
					result.normalized_scheduled_soc_profile(tuck) -= result.capacity_scale;
				}				
			}
			
			Allowed_Exchange.row(maximum_price_diff_ID(0)) = Eigen::VectorXi::Zero(total_duration);
			Allowed_Exchange.col(maximum_price_diff_ID(0)) = Eigen::VectorXi::Zero(total_duration);
			Allowed_Exchange.row(maximum_price_diff_ID(1)) = Eigen::VectorXi::Zero(total_duration);
			Allowed_Exchange.col(maximum_price_diff_ID(1)) = Eigen::VectorXi::Zero(total_duration);
		}
	}
}

void EV_schedule(Eigen::VectorXd subscept_tariff, EV_inform &result){
	// Check if EV currently undergoes standby period at house
	if(result.house_default_period(0) == 1){
		int tick = 0;
		while(result.house_default_period(tick) == 1){
			tick += 1;
		}
		
		result.BESS.soc_final = result.BESS.energy_scale;
		sorted_vector sorted_tariff = sort(subscept_tariff.head(tick));
		storage_schedule(sorted_tariff, result.BESS);
	}
}

int main(){
	// Test case Initialization
	Eigen::VectorXd subscept_tariff(10);
	subscept_tariff << 6, 6, 9, 2, 6, 1, 6, 9, 10, 2;
	sorted_vector sorted_tariff = sort(subscept_tariff);
	end_user_operation test_user;
	test_user.normalized_default_demand_profile = Eigen::VectorXd(subscept_tariff.size());
	test_user.normalized_default_demand_profile << 1, 1.2, 1.5, .8, .6, 2, 3, .1, .5, 1.2;
	
	// Smart appliance test
	test_user.smart_appliance.scale = .2;
	test_user.smart_appliance.flexibility_factor = .5;
	smart_appliance_schedule(sorted_tariff, test_user.normalized_default_demand_profile, 0, test_user.smart_appliance);
	std::cout << test_user.smart_appliance.normalized_scheduled_profile.transpose() << "\n" << std::endl;
	
	// BESS test, naive
	test_user.BESS.energy_scale = 4;
	test_user.BESS.capacity_scale = 1;
	test_user.BESS.efficiency = .95;
	test_user.BESS.soc_ini = 2;
	test_user.BESS.soc_final = 2;
	storage_schedule(sorted_tariff, test_user.BESS);
	std::cout << test_user.BESS.normalized_scheduled_capacity_profile.transpose() << "\n" << std::endl;
	
	// EV test
	test_user.EV.energy_demand = 1;
	test_user.EV.BESS = test_user.BESS;
	test_user.EV.usage_default_period = Eigen::VectorXi::Zero(subscept_tariff.size());
	test_user.EV.usage_default_period(3) = 1;
	test_user.EV.usage_default_period(6) = 1;
	test_user.EV.house_default_period = Eigen::VectorXi::Zero(subscept_tariff.size());
	test_user.EV.house_default_period.head(3) << 1, 1, 1;
	test_user.EV.house_default_period.tail(3) << 1, 1, 1;
	EV_schedule(subscept_tariff, test_user.EV);
	std::cout << test_user.EV.BESS.normalized_scheduled_capacity_profile.transpose() << std::endl;
}