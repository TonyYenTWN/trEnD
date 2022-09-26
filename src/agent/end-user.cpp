// Source file for end-user operational strategy
//#include <iostream>
//#include <chrono>
//#include "../basic/Basic_Definitions.h"
//#include "../basic/LP_gpa.cpp"
//#include "../basic/LP_gpa_fast.cpp"
#include "end-user.h"

void agent::end_user::end_user_LP_set(profile &profile){
	// -------------------------------------------------------------------------------
	// LP Solver initialization for operation strategy of end-users
	// Warm-up once and reuse for the rest of time slices
	// Variables are sorted as the following order (repeat for every t)
	// U^s(t), U^d(t): aggregated scheduled supply and demand
	// U^b(t): electricity demand from BESS
	// U^ev(t): electricity demand from EV
	// U^sa(t): electricity demand from smart appliance
	// RL_0(t): default residual load
	// ch^b(t): electricity charge from BESS
	// dc^b(t): electricity discharge from BESS
	// d^b(t): self consumption from BESS
	// soc^b(t): state of charge of BESS
	// ch^ev(t): electricity charge from EV
	// dc^ev(t): electricity discharge from EV
	// d^ev(t): electricity demand from EV
	// soc^ev(t): state of charge of EV
	// d^sa(..., t): scheduled smart appliance load at time t
	// d^sa(t): default smart appliance load at time t
	// -------------------------------------------------------------------------------
	int foresight_time = profile.operation.foresight_time;
	int load_shift_time = profile.operation.smart_appliance.shift_time;
//	load_shift_time = std::min(load_shift_time, foresight_time / 2);

	// -------------------------------------------------------------------------------
	// Set matrix for general constraints
	// -------------------------------------------------------------------------------
	// Generate sparse matrix for general equality constraints of dynamic equation
	int constrant_num = 7 * foresight_time + load_shift_time;
	int variable_per_time = 14 + foresight_time + load_shift_time;
	int variable_num = variable_per_time * foresight_time + foresight_time + load_shift_time;
	Eigen::VectorXpd non_zero_num(constrant_num);
	non_zero_num.head(5 * foresight_time) << Eigen::VectorXpd::Constant(foresight_time, 6), Eigen::VectorXpd::Constant(foresight_time, 3), 4, Eigen::VectorXpd::Constant(foresight_time - 1, 5), Eigen::VectorXpd::Constant(foresight_time, 3), 4, Eigen::VectorXpd::Constant(foresight_time - 1, 5);
	non_zero_num.segment(5 * foresight_time, foresight_time - load_shift_time) =  Eigen::VectorXpd::Constant(foresight_time - load_shift_time, 2 * load_shift_time + 1);
	if(load_shift_time > 0){
		for(int tick = 0; tick < load_shift_time; ++ tick){
			int row_ID = 6 * foresight_time - load_shift_time + tick;
			non_zero_num(row_ID) = 2 * load_shift_time - tick;
		}

		for(int tick = 0; tick < foresight_time + load_shift_time; ++ tick){
			int row_ID = 6 * foresight_time + tick;
			int tock_start = std::max(tick - 2 * load_shift_time, 0);
			int tock_end = std::min(tick, foresight_time - 1);
			non_zero_num(row_ID) = tock_end - tock_start + 1;
		}
	}
	non_zero_num.segment(5 * foresight_time, foresight_time) += Eigen::VectorXpd::Ones(foresight_time);
	non_zero_num.tail(foresight_time + load_shift_time) += Eigen::VectorXpd::Ones(foresight_time + load_shift_time);
	alglib::integer_1d_array row_sizes_general;
	row_sizes_general.setcontent(non_zero_num.size(), non_zero_num.data());
	alglib::sparsematrix constraint_general;
	alglib::sparsecreatecrs(constrant_num, variable_num, row_sizes_general, constraint_general);

	// Fill in the coefficients for the sparse matrix
	// U^d(t) - U^s(t) - U^b(t) - U^ev(t) - U^sa(t) - RL_0(t) = 0
	for(int row_iter = 0; row_iter < foresight_time; ++ row_iter){
		int U_d_ID = row_iter * variable_per_time;
		int U_s_ID = row_iter * variable_per_time + 1;
		int U_b_ID = row_iter * variable_per_time + 2;
		int U_ev_ID = row_iter * variable_per_time + 3;
		int U_sa_ID = row_iter * variable_per_time + 4;
		int RL_ID = row_iter * variable_per_time + 5;

		alglib::sparseset(constraint_general, row_iter, U_d_ID, 1.);
		alglib::sparseset(constraint_general, row_iter, U_s_ID, -1.);
		alglib::sparseset(constraint_general, row_iter, U_b_ID, -1.);
		alglib::sparseset(constraint_general, row_iter, U_ev_ID, -1.);
		alglib::sparseset(constraint_general, row_iter, U_sa_ID, -1.);
		alglib::sparseset(constraint_general, row_iter, RL_ID, -1.);
	}

	// U^b(t) - 1 / eta * ch^b(t) + eta * dc^b(t) = 0
	for(int row_iter = 0; row_iter < foresight_time; ++ row_iter){
		int row_ID = foresight_time + row_iter;
		int U_b_ID = row_iter * variable_per_time + 2;
		int ch_b_ID = row_iter * variable_per_time + 6;
		int dc_b_ID = row_iter * variable_per_time + 7;

		alglib::sparseset(constraint_general, row_ID, U_b_ID, 1.);
		alglib::sparseset(constraint_general, row_ID, ch_b_ID, -1. / profile.operation.BESS.efficiency);
		alglib::sparseset(constraint_general, row_ID, dc_b_ID, profile.operation.BESS.efficiency);
	}

	// soc^b(t) - soc^b(t - 1) - ch^b(t) + dc^b(t) + d^b(t) = 0
	// when t = 0, soc(0) also stored in d^b(t)
	{
		int row_iter = 0;
		int row_ID = 2 * foresight_time + row_iter;
		int ch_b_ID = row_iter * variable_per_time + 6;
		int dc_b_ID = row_iter * variable_per_time + 7;
		int d_b_ID = row_iter * variable_per_time + 8;
		int s_b_now_ID = row_iter * variable_per_time + 9;
		alglib::sparseset(constraint_general, row_ID, ch_b_ID, -1);
		alglib::sparseset(constraint_general, row_ID, dc_b_ID, 1.);
		alglib::sparseset(constraint_general, row_ID, d_b_ID, 1.);
		alglib::sparseset(constraint_general, row_ID, s_b_now_ID, 1.);
	}
	for(int row_iter = 1; row_iter < foresight_time; ++ row_iter){
		int row_ID = 2 * foresight_time + row_iter;
		int s_b_prev_ID = (row_iter - 1) * variable_per_time + 9;
		int ch_b_ID = row_iter * variable_per_time + 6;
		int dc_b_ID = row_iter * variable_per_time + 7;
		int d_b_ID = row_iter * variable_per_time + 8;
		int s_b_now_ID = row_iter * variable_per_time + 9;

		alglib::sparseset(constraint_general, row_ID, s_b_prev_ID, -1.);
		alglib::sparseset(constraint_general, row_ID, ch_b_ID, -1);
		alglib::sparseset(constraint_general, row_ID, dc_b_ID, 1.);
		alglib::sparseset(constraint_general, row_ID, d_b_ID, 1.);
		alglib::sparseset(constraint_general, row_ID, s_b_now_ID, 1.);
	}

	// U^ev(t) - 1 / eta * ch^ev(t) + eta * dc^ev(t) = 0
	for(int row_iter = 0; row_iter < foresight_time; ++ row_iter){
		int row_ID = 3 * foresight_time + row_iter;
		int U_ev_ID = row_iter * variable_per_time + 3;
		int ch_ev_ID = row_iter * variable_per_time + 10;
		int dc_ev_ID = row_iter * variable_per_time + 11;

		alglib::sparseset(constraint_general, row_ID, U_ev_ID, 1.);
		alglib::sparseset(constraint_general, row_ID, ch_ev_ID, -1. / profile.operation.EV.BESS.efficiency);
		alglib::sparseset(constraint_general, row_ID, dc_ev_ID, profile.operation.EV.BESS.efficiency);
	}

	// soc^ev(t) - soc^ev(t - 1) - ch^ev(t) + dc^ev(t) + d^ev(t) = 0
	// when t = 0, soc(0) also stored in d^ev(t)
	{
		int row_iter = 0;
		int row_ID = 4 * foresight_time + row_iter;
		int ch_ev_ID = row_iter * variable_per_time + 10;
		int dc_ev_ID = row_iter * variable_per_time + 11;
		int d_ev_ID = row_iter * variable_per_time + 12;
		int s_ev_now_ID = row_iter * variable_per_time + 13;

		alglib::sparseset(constraint_general, row_ID, ch_ev_ID, -1);
		alglib::sparseset(constraint_general, row_ID, dc_ev_ID, 1.);
		alglib::sparseset(constraint_general, row_ID, d_ev_ID, 1.);
		alglib::sparseset(constraint_general, row_ID, s_ev_now_ID, 1.);
	}
	for(int row_iter = 1; row_iter < foresight_time; ++ row_iter){
		int row_ID = 4 * foresight_time + row_iter;
		int s_ev_prev_ID = (row_iter - 1) * variable_per_time + 13;
		int ch_ev_ID = row_iter * variable_per_time + 10;
		int dc_ev_ID = row_iter * variable_per_time + 11;
		int d_ev_ID = row_iter * variable_per_time + 12;
		int s_ev_now_ID = row_iter * variable_per_time + 13;

		alglib::sparseset(constraint_general, row_ID, s_ev_prev_ID, -1.);
		alglib::sparseset(constraint_general, row_ID, ch_ev_ID, -1);
		alglib::sparseset(constraint_general, row_ID, dc_ev_ID, 1.);
		alglib::sparseset(constraint_general, row_ID, d_ev_ID, 1.);
		alglib::sparseset(constraint_general, row_ID, s_ev_now_ID, 1.);
	}

	// U^sa(t) - \sum_{tau} d^sa(tau, t) = 0
	for(int row_iter = 0; row_iter < foresight_time; ++ row_iter){
		int row_ID = 5 * foresight_time + row_iter;
		int U_sa_ID = row_iter * variable_per_time + 4;
		alglib::sparseset(constraint_general, row_ID, U_sa_ID, 1.);

		for(int tick = 0; tick < 2 * load_shift_time + 1; ++ tick){
			int tick_ID = row_iter - load_shift_time + tick;
			if(tick_ID >= -load_shift_time && tick_ID < foresight_time){
				int d_sa_ID = row_iter * variable_per_time + 14 + tick_ID + load_shift_time;
				alglib::sparseset(constraint_general, row_ID, d_sa_ID, -1.);
			}
		}
	}

	// d^sa(t) - \sum_{tau} d^sa(t, tau) = 0
	for(int row_iter = 0; row_iter < foresight_time + load_shift_time; ++ row_iter){
		int row_ID = 6 * foresight_time + row_iter;

		for(int tick = 0; tick < 2 * load_shift_time + 1; ++ tick){
			int tick_ID = row_iter - 2 * load_shift_time + tick;
			if(tick_ID >= 0 && tick_ID < foresight_time){
				int d_sa_ID = tick_ID * variable_per_time + 14 + row_iter;
				alglib::sparseset(constraint_general, row_ID, d_sa_ID, -1.);
			}
		}

		int d_sa_total_ID = variable_per_time * foresight_time + row_iter;
		alglib::sparseset(constraint_general, row_ID, d_sa_total_ID, 1.);
	}

	// -------------------------------------------------------------------------------
	// Set bounds for general constraints
	// -------------------------------------------------------------------------------
	Eigen::MatrixXd bound_general = Eigen::MatrixXd::Zero(constrant_num, 2);
	alglib::real_1d_array lb_general;
	alglib::real_1d_array ub_general;
	lb_general.setcontent(bound_general.rows(), bound_general.col(0).data());
	ub_general.setcontent(bound_general.rows(), bound_general.col(1).data());

	// -------------------------------------------------------------------------------
	// Set the LP problem object
	// -------------------------------------------------------------------------------
	alglib::minlpcreate(variable_num, profile.operation.Problem);
	alglib::minlpsetlc2(profile.operation.Problem, constraint_general, lb_general, ub_general, constrant_num);
	alglib::minlpsetalgodss(profile.operation.Problem, 0.);
}

void agent::end_user::end_user_LP_optimize(int tick, profile &profile){

}

agent::sorted_vector agent::sort(Eigen::VectorXd original){
//	// Sort of vector
// 	std::vector<int> item_seq(original.size());
// 	std::iota(item_seq.begin(), item_seq.end(), 0);
// 	std::sort(item_seq.begin(), item_seq.end(), [&](int i,int j){return original(i) < original(j);});
//
// 	// Output of vector
// 	agent::sorted_vector result;
// 	result.id = Eigen::VectorXi(original.size());
// 	result.value = Eigen::VectorXd(original.size());
// 	for(int item_ID = 0; item_ID < original.size(); ++ item_ID){
// 		result.id(item_ID) = item_seq[item_ID];
// 		result.value(item_ID) = original(item_seq[item_ID]);
//	}
//
//	return(result);
}

void agent::end_user::smart_appliance_schedule(agent::sorted_vector sorted_tariff, Eigen::VectorXd normalized_default_demand_profile, agent::end_user::smart_appliance_inform &result){
//	// Initialization
//	double total_flex_demand_energy = result.scale * normalized_default_demand_profile.sum() + result.unfulfilled_demand;
//	double flex_demand_capacity_max = total_flex_demand_energy / normalized_default_demand_profile.size() / result.flexibility_factor;
//	int flex_demand_duration = int(result.flexibility_factor * normalized_default_demand_profile.size());
//	if(double(result.flexibility_factor * normalized_default_demand_profile.size() - flex_demand_duration) != 0){
//		flex_demand_duration += 1;
//	}
//
//	// Schedule the demand to low price periods
//	result.normalized_scheduled_profile = Eigen::VectorXd::Zero(normalized_default_demand_profile.size());
//	for(int tick = 0; tick < flex_demand_duration - 1; ++ tick){
//		result.normalized_scheduled_profile(sorted_tariff.id(tick)) = flex_demand_capacity_max;
//	}
//	result.normalized_scheduled_profile(sorted_tariff.id(flex_demand_duration - 1)) = total_flex_demand_energy - (flex_demand_duration - 1) * flex_demand_capacity_max;
}

alglib::minlpstate agent::end_user::storage_schedule_LP_mold(int foresight_time){
//	alglib::minlpstate Problem;
//
//	// -------------------------------------------------------------------------------
//	// LP Solver initialization for BESS schedule
//	// Warm-up once and reuse for the rest of time slices
//	// Variables are sorted as {s_i, q_self(0), q_sell(0), q_dc(0), q_ch(0), s(0), ...}
//	// -------------------------------------------------------------------------------
//
//	// -------------------------------------------------------------------------------
//	// Set matrix for general constraints
//	// -------------------------------------------------------------------------------
//	// Generate sparse matrix for general equality constraints of dynamic equation
//	int constrant_num = 2 * foresight_time;
//	int variable_num = 5 * foresight_time + 1;
//	Eigen::VectorXpd non_zero_num(constrant_num);
//	non_zero_num << Eigen::VectorXpd::Constant(foresight_time, 4), Eigen::VectorXpd::Constant(foresight_time, 3);
//	alglib::integer_1d_array row_sizes_general;
//	row_sizes_general.setcontent(non_zero_num.size(), non_zero_num.data());
//	alglib::sparsematrix constraint_general;
//	alglib::sparsecreatecrs(constrant_num, variable_num, row_sizes_general, constraint_general);
//
//	// Fill in the coefficients for the sparse matrix
//	// -s(t - 1) + q_dc(t) - q_ch(t) + s(t) = 0
//	for(int row_iter = 0; row_iter < foresight_time; ++ row_iter){
//		int s_prev_ID = 5 * row_iter;
//		int q_dc_ID = 5 * row_iter + 3;
//		int q_ch_ID = 5 * row_iter + 4;
//		int s_ID = 5 * row_iter + 5;
//
//		alglib::sparseset(constraint_general, row_iter, s_prev_ID, -1.);
//		alglib::sparseset(constraint_general, row_iter, q_dc_ID, 1.);
//		alglib::sparseset(constraint_general, row_iter, q_ch_ID, -1.);
//		alglib::sparseset(constraint_general, row_iter, s_ID, 1.);
//	}
//	// q_self(t) + q_sell(t) - q_dc(t) = 0
//	for(int row_iter = 0; row_iter < foresight_time; ++ row_iter){
//		int row_ID = foresight_time + row_iter;
//		int q_self_ID = 5 * row_iter + 1;
//		int q_sell_ID = 5 * row_iter + 2;
//		int q_dc_ID = 5 * row_iter + 3;
//
//		alglib::sparseset(constraint_general, row_ID, q_self_ID, 1.);
//		alglib::sparseset(constraint_general, row_ID, q_sell_ID, 1.);
//		alglib::sparseset(constraint_general, row_ID, q_dc_ID, -1.);
//	}
//
//	// -------------------------------------------------------------------------------
//	// Set bounds for general constraints
//	// -------------------------------------------------------------------------------
//	Eigen::MatrixXd bound_general = Eigen::MatrixXd::Zero(constrant_num, 2);
//	alglib::real_1d_array lb_general;
//	alglib::real_1d_array ub_general;
//	lb_general.setcontent(bound_general.rows(), bound_general.col(0).data());
//	ub_general.setcontent(bound_general.rows(), bound_general.col(1).data());
//
//	// -------------------------------------------------------------------------------
//	// Set the LP problem object
//	// -------------------------------------------------------------------------------
//	alglib::minlpcreate(variable_num, Problem);
//	alglib::minlpsetlc2(Problem, constraint_general, lb_general, ub_general, constrant_num);
//	//alglib::minlpsetalgoipm(Problem);
//	alglib::minlpsetalgodss(Problem, 0.);
//	return Problem;
}

void agent::end_user::storage_schedule_LP_optimize(int foresight_time, sorted_vector expected_price_sorted, storage_inform &result, bool fixed_end){
//	int variable_num = 5 * foresight_time + 1;
//
//	// -------------------------------------------------------------------------------
//	// Set bounds for box constraints
//	// -------------------------------------------------------------------------------
//	Eigen::MatrixXd bound_box(variable_num, 2);
//	bound_box.row(0) << result.soc_ini, result.soc_ini;
//	for(int tick = 0; tick < foresight_time; ++ tick){
//		int q_self_ID = 5 * tick + 1;
//		int q_sell_ID = 5 * tick + 2;
//		int q_dc_ID = 5 * tick + 3;
//		int q_ch_ID = 5 * tick + 4;
//		int s_ID = 5 * tick + 5;
//
//		bound_box.row(q_self_ID) << 0, result.capacity_scale;
//		bound_box.row(q_sell_ID) << 0, result.capacity_scale;
//		bound_box.row(q_dc_ID) << 0, result.capacity_scale;
//		bound_box.row(q_ch_ID) << 0, result.capacity_scale;
//		bound_box.row(s_ID) << 0, result.energy_scale;
//	}
//	bound_box.bottomRows(1) << fixed_end * result.soc_final, fixed_end * result.soc_final + (1 - fixed_end) * result.energy_scale;
//
//	// Bounds of box constraints
//	alglib::real_1d_array lb_box;
//	alglib::real_1d_array ub_box;
//	lb_box.setcontent(bound_box.rows(), bound_box.col(0).data());
//	ub_box.setcontent(bound_box.rows(), bound_box.col(1).data());
//	alglib::minlpsetbc(result.Problem, lb_box, ub_box);
//
//	// -------------------------------------------------------------------------------
//	// Set objective coefficients of variables
//	// -------------------------------------------------------------------------------
//	Eigen::VectorXd obj_vec = Eigen::VectorXd::Zero(variable_num);
//	for(int tick = 0; tick < foresight_time; ++ tick){
//		int tick_ID = expected_price_sorted.id[tick];
//		int q_self_ID = 5 * tick_ID + 1;
//		int q_sell_ID = 5 * tick_ID + 2;
//		int q_ch_ID = 5 * tick_ID + 4;
//
//		obj_vec(q_self_ID) = -expected_price_sorted.value[tick] * result.efficiency;
//		obj_vec(q_sell_ID) = -expected_price_sorted.value[tick] * result.efficiency;
//		obj_vec(q_ch_ID) = expected_price_sorted.value[tick] / result.efficiency;
//	}
//	alglib::real_1d_array obj_coeff;
//	obj_coeff.setcontent(obj_vec.size(), obj_vec.data());
//	alglib::minlpsetcost(result.Problem, obj_coeff);
//
//	// -------------------------------------------------------------------------------
//	// Solve the problem and store the results
//	// -------------------------------------------------------------------------------
//	alglib::minlpoptimize(result.Problem);
//	alglib::real_1d_array sol;
//	alglib::minlpreport rep;
//	alglib::minlpresults(result.Problem, sol, rep);
//	Eigen::VectorXd sol_vec = Eigen::Map <Eigen::VectorXd> (sol.getcontent(), sol.length());
//
//	result.normalized_scheduled_capacity_profile = Eigen::VectorXd::Zero(foresight_time);
//	result.normalized_scheduled_soc_profile = Eigen::VectorXd::Zero(foresight_time);
//	for(int tick = 0; tick < foresight_time; ++ tick){
//		int q_dc_ID = 5 * tick + 3;
//		int q_ch_ID = 5 * tick + 4;
//		int s_ID = 5 * tick + 5;
//
//		result.normalized_scheduled_capacity_profile(tick) = sol_vec(q_dc_ID) - sol_vec(q_ch_ID);
//		result.normalized_scheduled_soc_profile(tick) = sol_vec(s_ID);
//		//std::cout << result.normalized_scheduled_capacity_profile.transpose() << "\n\n";
//	}
}

void agent::end_user::EV_schedule(int foresight_time, sorted_vector expected_price_sorted, EV_inform &result){
//	// Check if EV currently undergoes standby period at house
//	if(result.house_default_period(0) == 1){
//		// If currently at standby mode, find remaining duration of the mode
//		int tick = 0;
//		while(result.house_default_period(tick) == 1){
//			tick += 1;
//		}
//
//		// Set fixed end and optimize BESS schedule of the EV
//		result.BESS.soc_final = result.BESS.energy_scale;
//		result.BESS.Problem = storage_schedule_LP_mold(tick);
//		storage_schedule_LP_optimize(tick, expected_price_sorted, result.BESS, 1);
//	}
}
