// Source file for end-user operational strategy
#include <iostream>
//#include <chrono>
#include "../basic/Basic_Definitions.h"
#include "../basic/LP_GPA.cpp"

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

void storage_schedule_LP(Eigen::VectorXd subscept_tariff, storage_inform &result){
	LP_object Problem;
	
	// Set dimension of the problem
	Problem.Constraints_eq_num = subscept_tariff.size();
	Problem.Constraints_ie_num = 0;
	Problem.Variables_num = 3 * subscept_tariff.size();
	
	// Set objective vector
	// The variables are ordered as {q_dc(0), q_ch(0), s(0), q_dc(1), q_ch(1), s(1), ...}
	Problem.Objective.orig_vector = Eigen::VectorXd::Zero(Problem.Variables_num);
	for(int tick = 0; tick < subscept_tariff.size(); ++ tick){
		Problem.Objective.orig_vector(3 * tick) = subscept_tariff(tick) * result.efficiency;
		Problem.Objective.orig_vector(3 * tick + 1) = -subscept_tariff(tick) / result.efficiency;
	}
	// Can add terminal value for final storage level in the future
	//std::cout << std::fixed << std::setprecision(3) << Problem.Objective.orig_vector.transpose() << "\n" << std::endl;
	
	// Set boudary values for equality and inequality constraints
	Problem.Boundary.eq_vector = Eigen::VectorXd::Zero(Problem.Constraints_eq_num);
	Problem.Boundary.eq_vector(0) = result.soc_ini;		// Initial storage level
	Problem.Boundary.ie_orig_matrix = Eigen::MatrixXd::Zero(Problem.Variables_num + Problem.Constraints_ie_num, 2);
	for(int tick = 0; tick < subscept_tariff.size(); ++ tick){
		Problem.Boundary.ie_orig_matrix(3 * tick, 1) = result.capacity_scale;
		Problem.Boundary.ie_orig_matrix(3 * tick + 1, 1) = result.capacity_scale;
		Problem.Boundary.ie_orig_matrix(3 * tick + 2, 1) = result.energy_scale;
	}
	//std::cout << Problem.Boundary.eq_vector.transpose() << "\n" << std::endl;
	//std::cout << Problem.Boundary.ie_matrix.transpose() << "\n" << std::endl;
	
	// Set sparse matrix for equality constraints
	std::vector<Trip> Constraint_eq_trip;
	Constraint_eq_trip.reserve(Problem.Constraints_eq_num * Problem.Variables_num);
	for(int tick = 0; tick < subscept_tariff.size(); ++ tick){
		if(tick > 0){
			Constraint_eq_trip.push_back(Trip(tick, 3 * tick - 1, -1));
		}
		Constraint_eq_trip.push_back(Trip(tick, 3 * tick, 1));
		Constraint_eq_trip.push_back(Trip(tick, 3 * tick + 1, -1));
		Constraint_eq_trip.push_back(Trip(tick, 3 * tick + 2, 1));
	}
	Problem.Constraint.eq_orig_matrix = Eigen::SparseMatrix <double> (Problem.Constraints_eq_num, Problem.Variables_num);
	Problem.Constraint.eq_orig_matrix.setFromTriplets(Constraint_eq_trip.begin(), Constraint_eq_trip.end());
	//std::cout << Problem.Constraint.eq_matrix << "\n" << std::endl;
	
	// Set sparse matrix for original inequality constraints
	Problem.Constraint.ie_orig_matrix = Eigen::SparseMatrix <double> (Problem.Variables_num + Problem.Constraints_ie_num, Problem.Variables_num);
	std::vector<Trip> Constraint_ie_trip;
	Constraint_ie_trip.reserve((Problem.Constraints_ie_num + 1) * Problem.Variables_num);
	for(int var_id = 0; var_id < Problem.Variables_num; ++ var_id){
		Constraint_ie_trip.push_back(Trip(var_id, var_id, 1));
	}
	Problem.Constraint.ie_orig_matrix.setFromTriplets(Constraint_ie_trip.begin(), Constraint_ie_trip.end());
	//std::cout << Problem.Constraint.ie_matrix << std::endl;
	
	// Set initial feasible solutions and also reference point for contraction
	// The first feasible solution is where soc goes towards maximum asap and stays there; then q_dc(t) = q_ch(t) = 0
	// The second feasible solution is where soc goes towards minimum asap and stays there; then q_dc(t) / eff = q_ch(t) * eff = maximum possible value
	int transient_time_start = int (std::min(result.energy_scale - result.soc_ini, result.soc_ini)  / result.capacity_scale);
	transient_time_start = std::min(transient_time_start, int (subscept_tariff.size() / 2));
	Eigen::VectorXd feasible_point_1 = Eigen::VectorXd::Zero(Problem.Variables_num);
	Eigen::VectorXd feasible_point_2 = Eigen::VectorXd::Zero(Problem.Variables_num);
	if(transient_time_start > 0){
		feasible_point_1(1) = result.capacity_scale;
		feasible_point_2(0) = result.capacity_scale;
		feasible_point_1(2) = result.soc_ini + result.capacity_scale;
		feasible_point_2(2) = result.soc_ini - result.capacity_scale;
	}
	else{
		feasible_point_1(2) = result.soc_ini;
		feasible_point_2(2) = result.soc_ini;
	}
	for(int tick = 1; tick < subscept_tariff.size(); ++ tick){
		if(tick < transient_time_start){
			feasible_point_1(3 * tick + 1) = result.capacity_scale;
			feasible_point_2(3 * tick) = result.capacity_scale;
			feasible_point_1(3 * tick + 2) = feasible_point_1(3 * tick - 1) + result.capacity_scale;
			feasible_point_2(3 * tick + 2) = feasible_point_2(3 * tick - 1) - result.capacity_scale;			
		}
		else{
			feasible_point_2(3 * tick) = result.capacity_scale;
			feasible_point_2(3 * tick + 1) = result.capacity_scale;
			feasible_point_1(3 * tick + 2) = feasible_point_1(3 * tick - 1);
			feasible_point_2(3 * tick + 2) = feasible_point_2(3 * tick - 1);
		}
	}
	Problem.Solution.orig_vector = .5 * (feasible_point_1 + feasible_point_2);
	
	// Solve the LP and store the output schedule
	// Solution is degenerate; need to minimize energy loss due to efficiency by minimizing total amount of dc or ch
	double tol = pow(10, -12);
	LP_process(Problem, "BESS Schedule", 0);
	//std::cout << std::fixed << std::setprecision(3) << Problem.Solution.orig_vector.transpose() << std::endl;
	result.normalized_scheduled_capacity_profile = Eigen::VectorXd::Zero(subscept_tariff.size());
	result.normalized_scheduled_soc_profile = Eigen::VectorXd::Zero(subscept_tariff.size());
	for(int tick = 0; tick < subscept_tariff.size(); ++ tick){
		if(Problem.Solution.orig_vector(3 * tick) > tol){
			result.normalized_scheduled_capacity_profile(tick) = Problem.Solution.orig_vector(3 * tick) * result.efficiency;
		}
		else{
			result.normalized_scheduled_capacity_profile(tick) = -Problem.Solution.orig_vector(3 * tick + 1) / result.efficiency;
		}
		result.normalized_scheduled_soc_profile(tick) = Problem.Solution.orig_vector(3 * tick + 2);
	}
	//std::cout << result.normalized_scheduled_capacity_profile.transpose() << std::endl;
	//std::cout << result.normalized_scheduled_soc_profile.transpose() << std::endl;
}

void storage_schedule_naive(sorted_vector sorted_tariff, storage_inform &result){
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
		storage_schedule_naive(sorted_tariff, result.BESS);
	}
}

int main(){
	// Test case Initialization
	Eigen::VectorXd subscept_tariff(9);
	subscept_tariff << 1, 3, 1, -1, -1, -.5, 2, 2.5, 1.5;
	sorted_vector sorted_tariff = sort(subscept_tariff);
	end_user_operation test_user;
	test_user.normalized_default_demand_profile = Eigen::VectorXd(subscept_tariff.size());
	test_user.normalized_default_demand_profile << 1, 1.2;
	
	// Smart appliance test
	test_user.smart_appliance.scale = .2;
	test_user.smart_appliance.flexibility_factor = .5;
	//smart_appliance_schedule(sorted_tariff, test_user.normalized_default_demand_profile, 0, test_user.smart_appliance);
	//std::cout << test_user.smart_appliance.normalized_scheduled_profile.transpose() << "\n" << std::endl;
	
	// BESS test, naive
	test_user.BESS.energy_scale = 4;
	test_user.BESS.capacity_scale = 1;
	test_user.BESS.efficiency = .95;
	test_user.BESS.soc_ini = 2;
	test_user.BESS.soc_final = 2;
	storage_schedule_naive(sorted_tariff, test_user.BESS);
	std::cout << test_user.BESS.normalized_scheduled_capacity_profile.transpose() << "\n" << std::endl;
	storage_schedule_LP(subscept_tariff, test_user.BESS);
	std::cout << test_user.BESS.normalized_scheduled_capacity_profile.transpose() << "\n" << std::endl;
	
	// EV test
//	test_user.EV.energy_demand = 1;
//	test_user.EV.BESS = test_user.BESS;
//	test_user.EV.usage_default_period = Eigen::VectorXi::Zero(subscept_tariff.size());
//	test_user.EV.usage_default_period(3) = 1;
//	test_user.EV.usage_default_period(6) = 1;
//	test_user.EV.house_default_period = Eigen::VectorXi::Zero(subscept_tariff.size());
//	test_user.EV.house_default_period.head(3) << 1, 1, 1;
//	test_user.EV.house_default_period.tail(3) << 1, 1, 1;
//	EV_schedule(subscept_tariff, test_user.EV);
	//std::cout << test_user.EV.BESS.normalized_scheduled_capacity_profile.transpose() << std::endl;
}