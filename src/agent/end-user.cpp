// Source file for end-user operational strategy
#include "src/power_market/power_market.h"
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
		alglib::sparseset(constraint_general, row_ID, ch_b_ID, -1.);
		alglib::sparseset(constraint_general, row_ID, dc_b_ID, 1.);
		alglib::sparseset(constraint_general, row_ID, d_b_ID, 1.);
		alglib::sparseset(constraint_general, row_ID, s_b_now_ID, 1.);
//		std::cout << row_ID << ":\t";
//		std::cout << ch_b_ID << "\t";
//		std::cout << dc_b_ID << "\t";
//		std::cout << d_b_ID << "\t";
//		std::cout << s_b_now_ID  << "\n";
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
	int foresight_time = profile.operation.foresight_time;
	int load_shift_time = profile.operation.smart_appliance.shift_time;
	int price_interval = power_market::parameters::price_interval();
	int variable_per_time = 14 + foresight_time + load_shift_time;
	int variable_num = variable_per_time * foresight_time + foresight_time + load_shift_time;
	power_market::parameters::price_ID_bimap bidded_price_map;
	power_market::parameters::bidded_price(bidded_price_map);

	// -------------------------------------------------------------------------------
	// Set bounds for box constraints
	// -------------------------------------------------------------------------------
	Eigen::MatrixXd bound_box = Eigen::MatrixXd::Zero(variable_num, 2);
	for(int tock = 0; tock < foresight_time; ++ tock){
		int row_start = tock * variable_per_time;
		int U_d_ID = row_start;
		int U_s_ID = row_start + 1;
		int U_b_ID = row_start + 2;
		int U_ev_ID = row_start + 3;
		int U_sa_ID = row_start + 4;
		int RL_ID = row_start + 5;
		int ch_b_ID = row_start + 6;
		int dc_b_ID = row_start + 7;
		int d_b_ID = row_start + 8;
		int s_b_ID = row_start + 9;
		int ch_ev_ID = row_start + 10;
		int dc_ev_ID = row_start + 11;
		int d_ev_ID = row_start + 12;
		int s_ev_ID = row_start + 13;

//		bound_box.row(U_d_ID) << -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();
//		bound_box.row(U_s_ID) << -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();
		bound_box.row(U_d_ID) << 0, std::numeric_limits<double>::infinity();
		bound_box.row(U_s_ID) << 0, std::numeric_limits<double>::infinity();
		bound_box.row(U_b_ID) << -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();
		bound_box.row(U_ev_ID) << -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();
		bound_box.row(U_sa_ID) << -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();
		bound_box.row(RL_ID) << Eigen::RowVector2d::Constant(profile.operation.default_demand_profile(tock) - profile.operation.default_PV_profile(tock));
		bound_box.row(ch_b_ID) << 0., profile.operation.BESS.capacity_scale;
		bound_box.row(ch_b_ID) *= profile.investment.decision.BESS;
		bound_box.row(dc_b_ID) << 0., profile.operation.BESS.capacity_scale;
		bound_box.row(dc_b_ID) *= profile.investment.decision.BESS;
		bound_box.row(s_b_ID) << 0., profile.operation.BESS.energy_scale;
		bound_box.row(d_b_ID) << Eigen::RowVector2d::Constant(profile.operation.BESS.self_consumption);
		bound_box.row(ch_ev_ID) << 0., profile.operation.EV.BESS.capacity_scale;
		bound_box.row(ch_ev_ID) *= profile.investment.decision.EV_self_charging;
		bound_box.row(dc_ev_ID) << 0., profile.operation.EV.BESS.capacity_scale;
		bound_box.row(dc_ev_ID) *= profile.investment.decision.EV_self_charging;
		bound_box.row(s_ev_ID) << 0., profile.operation.EV.BESS.energy_scale;
		bound_box.row(d_ev_ID) << Eigen::RowVector2d::Constant(profile.operation.EV.default_demand_profile(tock) + profile.operation.EV.BESS.self_consumption);
		if(tock == 0){
			bound_box.row(d_b_ID) -= Eigen::RowVector2d::Constant(profile.operation.BESS.soc);
			bound_box.row(d_ev_ID) -= Eigen::RowVector2d::Constant(profile.operation.EV.BESS.soc);
		}

		for(int tuck = 0; tuck < 2 * load_shift_time + 1; ++ tuck){
			int tuck_ID = tock - load_shift_time + tuck;
			if(tuck_ID >= -load_shift_time && tuck_ID < foresight_time){
				int d_sa_ID = tock * variable_per_time + 14 + tuck_ID + load_shift_time;
				bound_box.row(d_sa_ID) << 0., profile.operation.smart_appliance.unfulfilled_demand(tuck_ID + load_shift_time);
				//std::cout << bound_box.row(d_sa_ID) << "\n";
			}
		}
		//std::cout << bound_box.row(s_b_ID) << "\n";
		//std::cout << "\n";
	}
	//std::cout << "\n";

	for(int tock = 0; tock < foresight_time + load_shift_time; ++ tock){
		int d_sa_total_ID = variable_per_time * foresight_time + tock;
		bound_box.row(d_sa_total_ID) = Eigen::RowVector2d::Constant(profile.operation.smart_appliance.unfulfilled_demand(tock));
	}

	// Bounds of box constraints
	alglib::real_1d_array lb_box;
	alglib::real_1d_array ub_box;
	lb_box.setcontent(bound_box.rows(), bound_box.col(0).data());
	ub_box.setcontent(bound_box.rows(), bound_box.col(1).data());
	alglib::minlpsetbc(profile.operation.Problem, lb_box, ub_box);

	// -------------------------------------------------------------------------------
	// Set objective coefficients of variables
	// -------------------------------------------------------------------------------
	Eigen::VectorXd obj_vec = Eigen::VectorXd::Zero(variable_num);
	for(int tock = 0; tock < foresight_time; ++ tock){
		int row_start = tock * variable_per_time;
		int U_d_ID = row_start;
		int U_s_ID = row_start + 1;

		obj_vec(U_d_ID) = profile.operation.price_demand_profile(tock);
		obj_vec(U_s_ID) = -profile.operation.price_supply_profile(tock);
	}
	alglib::real_1d_array obj_coeff;
	obj_coeff.setcontent(obj_vec.size(), obj_vec.data());
	alglib::minlpsetcost(profile.operation.Problem, obj_coeff);

	// -------------------------------------------------------------------------------
	// Solve the problem and store the results
	// -------------------------------------------------------------------------------
	alglib::minlpoptimize(profile.operation.Problem);
	alglib::real_1d_array sol;
	alglib::minlpreport rep;
	alglib::minlpresults(profile.operation.Problem, sol, rep);
	profile.operation.BESS.scheduled_capacity = sol[6] - sol[7];
	profile.operation.EV.BESS.scheduled_capacity = sol[10] - sol[11];
	profile.operation.smart_appliance.scheduled_demand = Eigen::VectorXd::Zero(2 * load_shift_time + 1);
	for(int tock = 0; tock < 2 * load_shift_time + 1; ++ tock){
		profile.operation.smart_appliance.scheduled_demand(tock) = sol[14 + tock];
	}
//	std::cout << rep.terminationtype << "\n";
//	std::cout << sol[0] << "\t";
//	std::cout << sol[1] << "\t";
//	std::cout << sol[2] << "\t";
//	std::cout << sol[3] << "\t";
//	std::cout << sol[4] << "\t";
//	std::cout << sol[5] << "\t";
//	std::cout << sol[6] << "\t";
//	std::cout << sol[7] << "\t";
//	std::cout << sol[8] << "\t";
//	std::cout << sol[9] << "\n";
//	std::cout << sol[variable_per_time] << "\t";
//	std::cout << sol[variable_per_time + 1] << "\t";
//	std::cout << sol[variable_per_time + 2] << "\t";
//	std::cout << sol[variable_per_time + 3] << "\t";
//	std::cout << sol[variable_per_time + 4] << "\t";
//	std::cout << sol[variable_per_time + 5] << "\t";
//	std::cout << sol[variable_per_time + 6] << "\t";
//	std::cout << sol[variable_per_time + 7] << "\t";
//	std::cout << sol[variable_per_time + 8] << "\t";
//	std::cout << sol[variable_per_time + 9] << "\n\n";
//	std::cout << -sol[6] + sol[7] + sol[8] + sol[9] << "\n";
//	std::cout << bidded_price_map.price_ID[profile.operation.price_demand_profile(0)] << "\n\n";

	int price_demand_inflex_ID = price_interval + 1;
	int price_demand_flex_ID = bidded_price_map.price_ID[profile.operation.price_demand_profile(0)];
	int price_supply_inflex_ID = 0;
	int price_supply_flex_ID = bidded_price_map.price_ID[profile.operation.price_supply_profile(0)];
	profile.operation.bids.submitted_demand_inflex(price_demand_inflex_ID) += std::max(sol[5], 0.);
	profile.operation.bids.submitted_supply_inflex(price_supply_inflex_ID) += -std::min(sol[5], 0.);
	profile.operation.bids.submitted_demand_inflex(price_demand_inflex_ID) += sol[14];
	//profile.operation.direct_demand = sol[0] - sol[6] / profile.operation.BESS.efficiency - sol[10] / profile.operation.EV.BESS.efficiency;
	if(profile.investment.decision.redispatch){
		profile.operation.bids.submitted_demand_flex(price_demand_flex_ID) += sol[4] - sol[14];
		profile.operation.bids.submitted_demand_flex(price_demand_flex_ID) += std::max(sol[2], 0.);
		profile.operation.bids.submitted_supply_flex(price_supply_flex_ID) += -std::min(sol[2], 0.);
		profile.operation.bids.submitted_demand_flex(price_demand_flex_ID) += std::max(sol[3], 0.);
		profile.operation.bids.submitted_supply_flex(price_supply_flex_ID) += -std::min(sol[3], 0.);
	}
	else{
		profile.operation.bids.submitted_demand_inflex(price_demand_flex_ID) += sol[4] - sol[14];
		profile.operation.bids.submitted_demand_inflex(price_demand_flex_ID) += std::max(sol[2], 0.);
		profile.operation.bids.submitted_supply_inflex(price_supply_flex_ID) += -std::min(sol[2], 0.);
		profile.operation.bids.submitted_demand_inflex(price_demand_flex_ID) += std::max(sol[3], 0.);
		profile.operation.bids.submitted_supply_inflex(price_supply_flex_ID) += -std::min(sol[3], 0.);
	}

//	std::cout << profile.operation.BESS.scheduled_capacity << "\t" << profile.operation.smart_appliance.scheduled_demand << "\n";
//	std::cout << rep.lagbc[7] << "\t" << rep.lagbc[14] << "\n\n";
//	if(abs(sol[2]) > 1E-12){
//		std::cout << sol[2] << "\t" << sol[41] << "\t" << sol[80] <<  "\t" << sol[119] << "\t" << sol[158] << "\t" << sol[197] << "\n";
//		std::cout << rep.lagbc[7] << "\t" << rep.lagbc[46] << "\t" << rep.lagbc[85] <<  "\t" << rep.lagbc[124] << "\t" << rep.lagbc[163] << "\t" << rep.lagbc[202] << "\n\n";
//	}
	//std::cout << sol[4] << "\t" << sol[43] << "\t" << sol[82] <<  "\t" << sol[121] << "\t" << sol[160] << "\n\n";
}
