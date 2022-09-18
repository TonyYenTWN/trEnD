// Source file for re-dispatch and tertiary control reserve market clearing of TSO in Norway
#include <iostream>
//#include <chrono>
//#include "../basic/LP_gpa.h"
#include "src/basic/rw_csv.h"
#include "power_market.h"

void power_market::TSO_Market_Set(market_inform &TSO_Market, power_network::network_inform &Power_network_inform, int Time){
//	double pi = boost::math::constants::pi<double>();

	// Input parameters of TSO market
	TSO_Market.num_zone = Power_network_inform.nodes.bidding_zone.size();
	TSO_Market.time_intervals = Time;
	TSO_Market.set_bidded_price();

	// Set node admittance matrix and line capacity matrix
	TSO_Market.network.num_vertice = TSO_Market.num_zone;
	Eigen::MatrixXd admittance = Eigen::MatrixXd::Zero(TSO_Market.network.num_vertice, TSO_Market.network.num_vertice);
	Eigen::MatrixXd capacity = Eigen::MatrixXd::Zero(TSO_Market.network.num_vertice, TSO_Market.network.num_vertice);
	for(int edge_iter = 0; edge_iter < Power_network_inform.edges.distance.size(); ++ edge_iter){
		int from_ID = Power_network_inform.edges.from(edge_iter );
		int to_ID = Power_network_inform.edges.to(edge_iter);
		int voltage = Power_network_inform.edges.voltage_base(edge_iter);
		double y_ij = 1. / Power_network_inform.edges.distance(edge_iter);
		y_ij /= Power_network_inform.tech_parameters.z_trans_series.imag();
		y_ij *= Power_network_inform.tech_parameters.impedenace_base_levels[voltage];
		admittance(from_ID, to_ID) += y_ij ;
		admittance(to_ID, from_ID) += y_ij ;
		capacity(from_ID, to_ID) += Power_network_inform.tech_parameters.power_limit[voltage];
	}

	// Set compact incidence matrix and edge admittance matrix
	double tol = 1E-6;
	TSO_Market.network.incidence.reserve(TSO_Market.network.num_vertice * TSO_Market.network.num_vertice);
	TSO_Market.network.admittance.reserve(TSO_Market.network.num_vertice * TSO_Market.network.num_vertice);
	std::vector <double> power_limit;
	power_limit.reserve(TSO_Market.network.num_vertice * TSO_Market.network.num_vertice);
	for(int row_iter = 0; row_iter < TSO_Market.network.num_vertice - 1; ++ row_iter){
		for(int col_iter = row_iter + 1; col_iter < TSO_Market.network.num_vertice; ++ col_iter){
			if(abs(admittance(row_iter , col_iter)) > tol){
				TSO_Market.network.incidence.push_back(Eigen::Vector2i(row_iter, col_iter));
				TSO_Market.network.admittance.push_back(admittance(row_iter , col_iter));
				power_limit.push_back(capacity(row_iter , col_iter));
			}
		}
	}
	TSO_Market.network.num_edges = TSO_Market.network.incidence.size();

	// Set voltage and power constraints at each edge
	TSO_Market.network.voltage_constraint = Eigen::MatrixXd::Ones(TSO_Market.network.num_vertice, 2);
	TSO_Market.network.voltage_constraint.col(0) *= -Power_network_inform.tech_parameters.theta_limit;
	TSO_Market.network.voltage_constraint.col(1) *= Power_network_inform.tech_parameters.theta_limit;
	TSO_Market.network.power_constraint = Eigen::MatrixXd (TSO_Market.network.num_edges, 2);
	TSO_Market.network.power_constraint.col(1) = Eigen::Map <Eigen::VectorXd> (power_limit.data(), power_limit.size());
	TSO_Market.network.power_constraint.col(1) /= Power_network_inform.tech_parameters.s_base;
	TSO_Market.network.power_constraint.col(0) = -TSO_Market.network.power_constraint.col(1);

	// Initialization of process variables
	power_market::Market_Initialization(TSO_Market);

	// Initialization of output variables
	TSO_Market.confirmed_supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_price = Eigen::MatrixXd(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_price_ratio = Eigen::MatrixXd(Time, TSO_Market.num_zone);
	TSO_Market.actual_supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.actual_demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.actual_price = Eigen::MatrixXd(Time, TSO_Market.num_zone);
	TSO_Market.actual_price_ratio = Eigen::MatrixXd(Time, TSO_Market.num_zone);
	TSO_Market.control_reserve.activated_positive = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.control_reserve.activated_negative = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.network.confirmed_voltage = Eigen::MatrixXd(Time, TSO_Market.network.num_vertice);
	TSO_Market.network.confirmed_power = Eigen::MatrixXd(Time, TSO_Market.network.num_edges);
//	TSO_Market.control_reserve.confirmed_positive = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
//	TSO_Market.control_reserve.confirmed_negative = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
}

void power_market::TSO_Market_Results_Get(int tick, market_inform &Market, alglib::minlpstate &Problem){
	alglib::real_1d_array sol;
	alglib::minlpreport rep;
	alglib::minlpresults(Problem, sol, rep);
	Eigen::VectorXd sol_vec = Eigen::Map <Eigen::VectorXd> (sol.getcontent(), sol.length());

	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		// Store power source / sink
		int row_start = 2 * Market.network.num_vertice + node_iter * (Market.price_intervals + 2);
		Market.confirmed_supply(tick, node_iter) = (sol_vec.segment(row_start, Market.price_intervals + 2).array().max(0)).sum();
		Market.confirmed_demand(tick, node_iter) = -(sol_vec.segment(row_start, Market.price_intervals + 2).array().min(0)).sum();

		// Store nodal prices
		Market.confirmed_price(tick, node_iter) = Market.bidded_price(0) + rep.lagbc[row_start];
		Market.confirmed_price(tick, node_iter) = std::min(Market.confirmed_price(tick, node_iter), Market.price_range_inflex(1));
		Market.confirmed_price(tick, node_iter) = std::max(Market.confirmed_price(tick, node_iter), Market.price_range_inflex(0));

		// Store ratio at nodes
		for(int price_iter = 0; price_iter < Market.price_intervals + 2; ++ price_iter){
			if(Market.bidded_price(price_iter) >= Market.confirmed_price(tick, node_iter) || price_iter == Market.price_intervals + 1){
				Market.confirmed_price_ratio(tick, node_iter) = sol[row_start + price_iter];
				if(sol[row_start + price_iter] >= 0.){
					Market.confirmed_price_ratio(tick, node_iter) /= Market.submitted_supply(price_iter, node_iter);
				}
				else{
					Market.confirmed_price_ratio(tick, node_iter) /= Market.submitted_demand(price_iter, node_iter);
				}
				break;
			}
		}
	}

	// Store voltage and power flow
	Market.network.confirmed_voltage.row(tick) = sol_vec.head(Market.network.num_vertice);
	Market.network.confirmed_power.row(tick) = sol_vec.tail(Market.network.num_edges);

	std::cout << Market.confirmed_supply.row(tick).sum() << "\t" << Market.confirmed_demand.row(tick).sum() << "\n";
	std::cout << sol_vec.segment(Market.network.num_vertice, Market.network.num_vertice).minCoeff() << " " << sol_vec.segment(Market.network.num_vertice, Market.network.num_vertice).maxCoeff() << " " << .5 * sol_vec.segment(Market.network.num_vertice, Market.network.num_vertice).array().abs().sum() << "\n";
	std::cout << sol_vec.head(Market.network.num_vertice).minCoeff() << " " << sol_vec.head(Market.network.num_vertice).maxCoeff()  << "\n";
	std::cout << sol_vec.tail(Market.network.num_edges).minCoeff() << " " << sol_vec.tail(Market.network.num_edges).maxCoeff() << "\n\n";
}

void power_market::TSO_Market_control_reserve(int tick, market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, bool DSO_filter_flag){
	int bz_num = Power_network_inform.points.bidding_zone.maxCoeff() + 1;
	int point_num = Power_network_inform.points.bidding_zone.size();

	// Temporary setting!!
	Power_market_inform.TSO_Market.control_reserve.available_ratio_supply = Eigen::MatrixXd::Ones(Power_market_inform.TSO_Market.price_intervals + 2, Power_market_inform.TSO_Market.num_zone);
	Power_market_inform.TSO_Market.control_reserve.available_ratio_demand = Eigen::MatrixXd::Ones(Power_market_inform.TSO_Market.price_intervals + 2, Power_market_inform.TSO_Market.num_zone);
	// Temporary setting!!

	Eigen::MatrixXd imbalance_demand = Eigen::MatrixXd::Zero(Power_market_inform.TSO_Market.price_intervals + 2, Power_market_inform.TSO_Market.num_zone);
	for(int point_iter = 0; point_iter < point_num; ++ point_iter){
		int point_ID = Power_network_inform.points.in_cluster_ID(point_iter);
		int node_ID = Power_network_inform.points.node(point_iter);
		int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
		double node_price = Power_market_inform.TSO_Market.confirmed_price(tick, node_ID);

		Eigen::VectorXd imbalance_demand_temp = Eigen::VectorXd::Zero(Power_market_inform.TSO_Market.price_intervals + 2);
		for(int price_iter = 0; price_iter < Power_market_inform.TSO_Market.price_intervals + 2; ++ price_iter){
			// Imbalance due to residential demand
			if(DSO_filter_flag){
				imbalance_demand_temp(price_iter) += Power_market_inform.DSO_Markets[DSO_ID].filtered_demand(price_iter, point_ID);
			}
			else{
				imbalance_demand_temp(price_iter) += Power_market_inform.DSO_Markets[DSO_ID].submitted_demand(price_iter, point_ID);
			}

			// Imbalance due to industrial demand
			imbalance_demand_temp(price_iter) += Power_market_inform.industrial_submitted_demand(price_iter, point_iter);

			if(Power_market_inform.TSO_Market.bidded_price(price_iter) >= Power_market_inform.TSO_Market.confirmed_price(tick, node_ID) || price_iter == Power_market_inform.TSO_Market.price_intervals + 1){
				if(Power_market_inform.TSO_Market.confirmed_price_ratio(tick, node_ID) < 0.){
					imbalance_demand_temp(price_iter) *= -Power_market_inform.TSO_Market.confirmed_price_ratio(tick, node_ID);
				}
				else{
					imbalance_demand_temp(price_iter) = 0.;
				}
				break;
			}
		}
		imbalance_demand_temp *= Power_network_inform.points.imbalance_field(point_iter, tick);
		imbalance_demand.col(node_ID) += imbalance_demand_temp;
	}

	// Exit the function if imbalance is too small
	if(imbalance_demand.array().abs().sum() < 1.){
		std::cout << "No imbalance in the time slice.\n\n";
		exit(0);
	}

	// -------------------------------------------------------------------------------
	// Update bounds for box constraints
	// -------------------------------------------------------------------------------
	int variable_num = Power_market_inform.TSO_Market.network.num_vertice * (Power_market_inform.TSO_Market.price_intervals + 4) + Power_market_inform.TSO_Market.network.num_edges;
	Eigen::MatrixXd bound_box(variable_num, 2);
	bound_box.topRows(Power_market_inform.TSO_Market.network.num_vertice) = Power_market_inform.TSO_Market.network.voltage_constraint;
	bound_box.middleRows(Power_market_inform.TSO_Market.network.num_vertice, Power_market_inform.TSO_Market.network.num_vertice).col(0) = Eigen::VectorXd::Constant(Power_market_inform.TSO_Market.network.num_vertice, -std::numeric_limits<double>::infinity());
	bound_box.middleRows(Power_market_inform.TSO_Market.network.num_vertice, Power_market_inform.TSO_Market.network.num_vertice).col(1) = Eigen::VectorXd::Constant(Power_market_inform.TSO_Market.network.num_vertice, std::numeric_limits<double>::infinity());
	for(int node_iter = 0; node_iter < Power_market_inform.TSO_Market.network.num_vertice; ++ node_iter){
		int row_start = 2 * Power_market_inform.TSO_Market.network.num_vertice + node_iter * (Power_market_inform.TSO_Market.price_intervals + 2);
		Eigen::VectorXd origin_point(Power_market_inform.TSO_Market.price_intervals + 2);
		Eigen::VectorXd flex_up(Power_market_inform.TSO_Market.price_intervals + 2);
		Eigen::VectorXd flex_down(Power_market_inform.TSO_Market.price_intervals + 2);
		bool full_supply = 1;

		for(int price_iter = 0; price_iter < Power_market_inform.TSO_Market.price_intervals + 2; ++ price_iter){
			double flex_quan;
			if(Power_market_inform.TSO_Market.bidded_price(price_iter) == Power_market_inform.TSO_Market.confirmed_price(tick, node_iter)){
				if(Power_market_inform.TSO_Market.confirmed_price_ratio(tick, node_iter) < 0.){
					origin_point(price_iter) = -Power_market_inform.TSO_Market.confirmed_demand(price_iter, node_iter) - imbalance_demand(price_iter, node_iter);
					flex_quan = Power_market_inform.TSO_Market.control_reserve.available_ratio_demand(price_iter, node_iter) * (Power_market_inform.TSO_Market.confirmed_demand(price_iter, node_iter) + imbalance_demand(price_iter, node_iter));
					flex_quan += Power_market_inform.TSO_Market.control_reserve.available_ratio_supply(price_iter, node_iter) * Power_market_inform.TSO_Market.submitted_supply(price_iter, node_iter);
					flex_up(price_iter) = flex_quan;

					flex_quan = Power_market_inform.TSO_Market.control_reserve.available_ratio_demand(price_iter, node_iter);
					flex_quan *= Power_market_inform.TSO_Market.submitted_demand(price_iter, node_iter) - Power_market_inform.TSO_Market.confirmed_demand(price_iter, node_iter) - imbalance_demand(price_iter, node_iter);
					flex_down(price_iter) = flex_quan;
				}
				else{
					origin_point(price_iter) = Power_market_inform.TSO_Market.confirmed_supply(price_iter, node_iter);
					flex_quan = Power_market_inform.TSO_Market.control_reserve.available_ratio_supply(price_iter, node_iter) * Power_market_inform.TSO_Market.confirmed_supply(price_iter, node_iter);
					flex_quan += Power_market_inform.TSO_Market.control_reserve.available_ratio_demand(price_iter, node_iter) * Power_market_inform.TSO_Market.submitted_demand(price_iter, node_iter);
					flex_down(price_iter) = flex_quan;

					flex_quan = Power_market_inform.TSO_Market.control_reserve.available_ratio_supply(price_iter, node_iter);
					flex_quan *= Power_market_inform.TSO_Market.submitted_supply(price_iter, node_iter) - Power_market_inform.TSO_Market.confirmed_supply(price_iter, node_iter);
					flex_up(price_iter) = flex_quan;
				}

				full_supply = 0;
				continue;
			}

			flex_quan = Power_market_inform.TSO_Market.control_reserve.available_ratio_demand(price_iter, node_iter) * (Power_market_inform.TSO_Market.submitted_demand(price_iter, node_iter) + imbalance_demand(price_iter, node_iter));
			flex_quan += Power_market_inform.TSO_Market.control_reserve.available_ratio_supply(price_iter, node_iter) * Power_market_inform.TSO_Market.submitted_supply(price_iter, node_iter);
			if(full_supply){
				origin_point(price_iter) = Power_market_inform.TSO_Market.submitted_supply(price_iter, node_iter);
				flex_up(price_iter) = 0.;
				flex_down(price_iter) = flex_quan;
			}
			else{
				origin_point(price_iter) = -Power_market_inform.TSO_Market.submitted_demand(price_iter, node_iter) - imbalance_demand(price_iter, node_iter);
				flex_up(price_iter) = flex_quan;
				flex_down(price_iter) = 0.;
			}
		}
		bound_box.middleRows(row_start, Power_market_inform.TSO_Market.price_intervals + 2).col(0) = origin_point - flex_down;
		bound_box.middleRows(row_start, Power_market_inform.TSO_Market.price_intervals + 2).col(1) = origin_point + flex_up;
	}
	bound_box.bottomRows(Power_market_inform.TSO_Market.network.num_edges) = Power_market_inform.TSO_Market.network.power_constraint;

	// Bounds of box constraints
	alglib::real_1d_array lb_box;
	alglib::real_1d_array ub_box;
	lb_box.setcontent(bound_box.rows(), bound_box.col(0).data());
	ub_box.setcontent(bound_box.rows(), bound_box.col(1).data());
	alglib::minlpsetbc(Power_market_inform.TSO_Problem, lb_box, ub_box);

	// -------------------------------------------------------------------------------
	// Solve the problem
	// -------------------------------------------------------------------------------
	alglib::minlpoptimize(Power_market_inform.TSO_Problem);

	alglib::real_1d_array sol;
	alglib::minlpreport rep;
	alglib::minlpresults(Power_market_inform.TSO_Problem, sol, rep);
	Eigen::VectorXd sol_vec = Eigen::Map <Eigen::VectorXd> (sol.getcontent(), sol.length());

	for(int node_iter = 0; node_iter < Power_market_inform.TSO_Market.network.num_vertice; ++ node_iter){
		// Store power source / sink
		int row_start = 2 * Power_market_inform.TSO_Market.network.num_vertice + node_iter * (Power_market_inform.TSO_Market.price_intervals + 2);
		Power_market_inform.TSO_Market.actual_supply(tick, node_iter) = (sol_vec.segment(row_start, Power_market_inform.TSO_Market.price_intervals + 2).array().max(0)).sum();
		Power_market_inform.TSO_Market.actual_demand(tick, node_iter) = -(sol_vec.segment(row_start, Power_market_inform.TSO_Market.price_intervals + 2).array().min(0)).sum();

		// Store nodal prices
		Power_market_inform.TSO_Market.actual_price(tick, node_iter) = Power_market_inform.TSO_Market.bidded_price(0) + rep.lagbc[row_start];
		Power_market_inform.TSO_Market.actual_price(tick, node_iter) = std::min(Power_market_inform.TSO_Market.actual_price(tick, node_iter), Power_market_inform.TSO_Market.price_range_inflex(1));
		Power_market_inform.TSO_Market.actual_price(tick, node_iter) = std::max(Power_market_inform.TSO_Market.actual_price(tick, node_iter), Power_market_inform.TSO_Market.price_range_inflex(0));

		// Store ratio at nodes
		for(int price_iter = 0; price_iter < Power_market_inform.TSO_Market.price_intervals + 2; ++ price_iter){
			if(Power_market_inform.TSO_Market.bidded_price(price_iter) >= Power_market_inform.TSO_Market.actual_price(tick, node_iter) || price_iter == Power_market_inform.TSO_Market.price_intervals + 1){
				Power_market_inform.TSO_Market.actual_price_ratio(tick, node_iter) =  sol[row_start + price_iter];
				if(sol[row_start + price_iter] >= 0.){
					Power_market_inform.TSO_Market.actual_price_ratio(tick, node_iter) /= Power_market_inform.TSO_Market.submitted_supply(price_iter, node_iter);
				}
				else{
					Power_market_inform.TSO_Market.actual_price_ratio(tick, node_iter) /= Power_market_inform.TSO_Market.submitted_demand(price_iter, node_iter) - imbalance_demand(price_iter, node_iter);
				}
				break;
			}
		}
	}

	// Store voltage and power flow
	Power_market_inform.TSO_Market.network.confirmed_voltage.row(tick) = sol_vec.head(Power_market_inform.TSO_Market.network.num_vertice);
	Power_market_inform.TSO_Market.network.confirmed_power.row(tick) = sol_vec.tail(Power_market_inform.TSO_Market.network.num_edges);

	std::cout << Power_market_inform.TSO_Market.actual_supply.row(tick).sum() << "\t" << Power_market_inform.TSO_Market.actual_demand.row(tick).sum() << "\n";
	std::cout << sol_vec.segment(Power_market_inform.TSO_Market.network.num_vertice, Power_market_inform.TSO_Market.network.num_vertice).minCoeff() << " " << sol_vec.segment(Power_market_inform.TSO_Market.network.num_vertice, Power_market_inform.TSO_Market.network.num_vertice).maxCoeff() << " " << .5 * sol_vec.segment(Power_market_inform.TSO_Market.network.num_vertice, Power_market_inform.TSO_Market.network.num_vertice).array().abs().sum() << "\n";
	std::cout << sol_vec.head(Power_market_inform.TSO_Market.network.num_vertice).minCoeff() << " " << sol_vec.head(Power_market_inform.TSO_Market.network.num_vertice).maxCoeff()  << "\n";
	std::cout << sol_vec.tail(Power_market_inform.TSO_Market.network.num_edges).minCoeff() << " " << sol_vec.tail(Power_market_inform.TSO_Market.network.num_edges).maxCoeff() << "\n\n";
}
