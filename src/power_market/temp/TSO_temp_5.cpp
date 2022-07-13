// Source file for re-dispatch and tertiary control reserve market clearing of TSO in Norway
#include <iostream>
#include <iomanip>
//#include <chrono>
#include "../../basic/LP_gpa.h"
#include "../../basic/rw_csv.cpp"
#include "../../power_network/power_network_input.cpp"
#include "../power_market.h"
#include "../power_market_func.cpp"

void TSO_Market_Set_Test_1(market_inform &TSO_Market, int Time){	
	// Input parameters of TSO market
	// A trivial test case with 3 nodes, the first bus being the reference connecting to the remaining 2
	TSO_Market.num_zone = 3;
	TSO_Market.time_intervals = Time;
	TSO_Market.price_intervals = 600;
	TSO_Market.price_range_inflex << -500., 3000.;
	TSO_Market.price_range_flex << -100., 500.;
	TSO_Market.bidded_price = Eigen::VectorXd(TSO_Market.price_intervals + 2);
	TSO_Market.bidded_price(0) = TSO_Market.price_range_inflex(0);
	TSO_Market.bidded_price.array().tail(1) = TSO_Market.price_range_inflex(1);
	TSO_Market.bidded_price.array().segment(1, TSO_Market.price_intervals) = Eigen::VectorXd::LinSpaced(TSO_Market.price_intervals, TSO_Market.price_range_flex(0) + .5 * (TSO_Market.price_range_flex(1) - TSO_Market.price_range_flex(0)) / TSO_Market.price_intervals, TSO_Market.price_range_flex(1) - .5 * (TSO_Market.price_range_flex(1) - TSO_Market.price_range_flex(0)) / TSO_Market.price_intervals);
	
	// Set compact node admittance matrix Y_n
	TSO_Market.network.num_vertice = TSO_Market.num_zone;
	TSO_Market.network.num_edges = 3;
	TSO_Market.network.incidence_matrix = Eigen::MatrixXi(TSO_Market.network.num_edges, 2);
	TSO_Market.network.incidence_matrix.row(0) << 0, 1;
	TSO_Market.network.incidence_matrix.row(1) << 0, 2;
	TSO_Market.network.incidence_matrix.row(2) << 1, 2;
	TSO_Market.network.admittance_vector = Eigen::VectorXd(3);
	TSO_Market.network.admittance_vector << 100., 100., 100.;
	
	// Set voltage and power constraints at each edge
	TSO_Market.network.voltage_constraint = Eigen::MatrixXd(TSO_Market.network.num_vertice, 2);
	TSO_Market.network.voltage_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, -.2);
	TSO_Market.network.voltage_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, .2);
	TSO_Market.network.power_constraint = Eigen::MatrixXd(TSO_Market.network.num_edges, 2);
	TSO_Market.network.power_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, -50.);
	TSO_Market.network.power_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, 50.);

	// Initialization of output variables
	TSO_Market.confirmed_supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_price = Eigen::MatrixXd(Time, TSO_Market.num_zone);
	TSO_Market.network.confirmed_power = Eigen::MatrixXd(Time, TSO_Market.network.num_edges);	
	
	// For the trivial case only: initialize submitted supply and demand bids at each node
	Market_Initialization(TSO_Market);
	TSO_Market.submitted_supply(0, 0) = 10.;
	TSO_Market.submitted_supply(0, 1) = 20.;
	TSO_Market.submitted_supply(0, 2) = 20.;
	TSO_Market.submitted_supply(1, 0) = 5.;
	TSO_Market.submitted_supply(1, 1) = 15.;
	TSO_Market.submitted_supply(1, 2) = 25.;
	TSO_Market.submitted_demand(TSO_Market.price_intervals + 1, 0) = 25.;
	TSO_Market.submitted_demand(TSO_Market.price_intervals + 1, 1) = 15.;
	TSO_Market.submitted_demand(TSO_Market.price_intervals + 1, 2) = 5.;
	TSO_Market.submitted_demand(TSO_Market.price_intervals, 0) = 10.;
	TSO_Market.submitted_demand(TSO_Market.price_intervals, 1) = 7.5;
	TSO_Market.submitted_demand(TSO_Market.price_intervals, 2) = 5.;
	//S = {15, 35, 45}; D = {35, 22.5, 10}
}

void TSO_Market_Set_Test_2(market_inform &TSO_Market, int Time){
	// Input parameters of TSO market
	// A trivial test case with 20 nodes connected as a radial line
	TSO_Market.num_zone = 1000;
	TSO_Market.time_intervals = Time;
	TSO_Market.price_intervals = 600;
	TSO_Market.price_range_inflex << -500., 3000.;
	TSO_Market.price_range_flex << -100., 500.;
	TSO_Market.bidded_price = Eigen::VectorXd(TSO_Market.price_intervals + 2);
	TSO_Market.bidded_price(0) = TSO_Market.price_range_inflex(0);
	TSO_Market.bidded_price.array().tail(1) = TSO_Market.price_range_inflex(1);
	TSO_Market.bidded_price.array().segment(1, TSO_Market.price_intervals) = Eigen::VectorXd::LinSpaced(TSO_Market.price_intervals, TSO_Market.price_range_flex(0) + .5 * (TSO_Market.price_range_flex(1) - TSO_Market.price_range_flex(0)) / TSO_Market.price_intervals, TSO_Market.price_range_flex(1) - .5 * (TSO_Market.price_range_flex(1) - TSO_Market.price_range_flex(0)) / TSO_Market.price_intervals);
	
	// Set node admittance matrix Y_n
	TSO_Market.network.num_vertice = TSO_Market.num_zone;
	TSO_Market.network.num_edges = TSO_Market.num_zone - 1;
	TSO_Market.network.incidence_matrix = Eigen::MatrixXi(TSO_Market.network.num_edges, 2);
	for(int edge_iter = 0; edge_iter < TSO_Market.network.num_edges; ++ edge_iter){
		TSO_Market.network.incidence_matrix.row(edge_iter) << edge_iter, edge_iter + 1;
	}
	TSO_Market.network.admittance_vector = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, 5000.);
	
	// Set voltage and power constraints at each edge
	TSO_Market.network.voltage_constraint = Eigen::MatrixXd(TSO_Market.network.num_vertice, 2);
	TSO_Market.network.voltage_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, -.1);
	TSO_Market.network.voltage_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, .1);
	TSO_Market.network.power_constraint = Eigen::MatrixXd(TSO_Market.network.num_edges, 2);
	TSO_Market.network.power_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, -150.);
	TSO_Market.network.power_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, 150.);

	// Initialization of output variables
	TSO_Market.confirmed_supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_price = Eigen::MatrixXd(Time, TSO_Market.num_zone);
	TSO_Market.network.confirmed_power = Eigen::MatrixXd(Time, TSO_Market.network.num_edges);
	
	// For the trivial case only: initialize submitted supply and demand bids at each node
	Market_Initialization(TSO_Market);
	TSO_Market.submitted_supply.leftCols(TSO_Market.num_zone / 2 - 1) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2 - 1, 1.);
	TSO_Market.submitted_demand.rightCols(TSO_Market.num_zone / 2 - 1) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2 - 1, 1.);
}

void Flow_Based_Market_Optimization_Test(int tick, market_inform &Market, LP_object &Problem){
	Eigen::MatrixXd bidded_supply = Market.submitted_supply;
	Eigen::MatrixXd bidded_demand = Market.submitted_demand;
	
	// Initial market clearing within each nodes
	Eigen::VectorXi price_ID(Market.num_zone);
	Market_clearing_nodal(tick, Market, price_ID, bidded_supply, bidded_demand);
	//std::cout << price_ID.transpose() << "\n\n";
	
	// Initialization of process variables for the main optimization loop
	Eigen::MatrixXd bidded_supply_export = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_demand_export = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_supply_import = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_demand_import = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_supply_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);
	Eigen::MatrixXd bidded_demand_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);
	Eigen::MatrixXd bidded_total_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);
	Eigen::MatrixXd utility_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);

	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
		// Supply curves
		bidded_supply_export.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter)) = Market.submitted_supply.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter));
		bidded_supply_export(price_ID(zone_iter), zone_iter) = bidded_supply(price_ID(zone_iter), zone_iter);
		bidded_supply_import.col(zone_iter).array().head(price_ID(zone_iter)) = Market.submitted_supply.col(zone_iter).array().head(price_ID(zone_iter));
		bidded_supply_import(price_ID(zone_iter), zone_iter) = Market.submitted_supply(price_ID(zone_iter), zone_iter) - bidded_supply(price_ID(zone_iter), zone_iter);
		bidded_supply_aggregated(0, zone_iter) = -bidded_supply_import.col(zone_iter).sum();
		
		// Demand curves
		bidded_demand_export.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter)) = Market.submitted_demand.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter));
		bidded_demand_export(price_ID(zone_iter), zone_iter) = Market.submitted_demand(price_ID(zone_iter), zone_iter) - bidded_demand(price_ID(zone_iter), zone_iter);
		bidded_demand_import.col(zone_iter).array().head(price_ID(zone_iter)) = Market.submitted_demand.col(zone_iter).array().head(price_ID(zone_iter));
		bidded_demand_import(price_ID(zone_iter), zone_iter) = bidded_demand(price_ID(zone_iter), zone_iter);
		bidded_demand_aggregated(0, zone_iter) = -bidded_demand_import.col(zone_iter).sum();		
	}
	
	// Aggregated import-export curves for supply and demand
	for(int price_iter = 1; price_iter < Market.price_intervals + 3; ++ price_iter){
		bidded_supply_aggregated.row(price_iter) = bidded_supply_aggregated.row(price_iter - 1) + bidded_supply_import.row(price_iter - 1) + bidded_supply_export.row(price_iter - 1);
		bidded_demand_aggregated.row(price_iter) = bidded_demand_aggregated.row(price_iter - 1) + bidded_demand_import.row(price_iter - 1) + bidded_demand_export.row(price_iter - 1);
	}
	bidded_total_aggregated = bidded_supply_aggregated + bidded_demand_aggregated;
	
	// Update inequality constraints for import / export
	Problem.Boundary.ie_orig_matrix.col(0).segment(Market.network.num_edges, Market.network.num_vertice) = bidded_total_aggregated.row(0).transpose();
	Problem.Boundary.ie_orig_matrix.col(1).segment(Market.network.num_edges, Market.network.num_vertice) = bidded_total_aggregated.row(Market.price_intervals + 2).transpose();

	// Aggregated utility function for each bidding zone
	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){	
		// Demand
		utility_aggregated(price_ID(zone_iter), zone_iter) = (bidded_demand_import(price_ID(zone_iter), zone_iter) + bidded_supply_import(price_ID(zone_iter), zone_iter)) * Market.bidded_price(price_ID(zone_iter));
		if(price_ID(zone_iter) > 0){
			for(int price_iter = price_ID(zone_iter) - 1; price_iter >= 0; -- price_iter){
				utility_aggregated(price_iter, zone_iter) = utility_aggregated(price_iter + 1, zone_iter) + (bidded_demand_import(price_iter, zone_iter) + bidded_supply_import(price_iter, zone_iter)) * Market.bidded_price(price_iter);
			}
		}
		
		// Supply
		utility_aggregated(price_ID(zone_iter) + 1, zone_iter) = -(bidded_demand_export(price_ID(zone_iter), zone_iter) + bidded_supply_export(price_ID(zone_iter), zone_iter)) * Market.bidded_price(price_ID(zone_iter));
		if(price_ID(zone_iter) < Market.price_intervals){
			for(int price_iter = price_ID(zone_iter) + 2; price_iter < Market.price_intervals + 3; ++ price_iter){
				utility_aggregated(price_iter, zone_iter) = utility_aggregated(price_iter - 1, zone_iter) - (bidded_demand_export(price_iter - 1, zone_iter) + bidded_supply_export(price_iter - 1, zone_iter)) * Market.bidded_price(price_iter - 1);
			}			
		}
	}

	// Find nodes with no power source / sink
	std::vector <int> Non_zero;
	Non_zero.reserve(Market.num_zone);
	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
		if(abs(bidded_total_aggregated(0, zone_iter) - bidded_total_aggregated(Market.price_intervals + 2, zone_iter)) > 0.){
			Non_zero.push_back(zone_iter);
		}
	}
	Problem.Solver.qr.compute(Problem.Constraint.eq_orig_matrix * Problem.Constraint.eq_orig_matrix.transpose());

	// Declare variables for the main loop
	bool infeasbile_flag;
	bool break_flag;
	int loop_count;
	double tol = pow(10., -12.);
	double eps = pow(10., -10.);
	double mu = eps;
	double dS_max = 100.;
	double dS;
	double obj;
	double obj_temp;
	double bar;
	double bar_temp;
	double price_penalty = 10000.;
	Eigen::VectorXi price_ID_temp = price_ID;
	Eigen::VectorXd utility = Eigen::VectorXd::Zero(Market.num_zone);
	Eigen::VectorXd utility_temp = utility;
	Eigen::VectorXd sol_temp = Problem.Solution.orig_vector;
	Eigen::VectorXd grad(Problem.Variables_num);
//	std::vector <Trip> Indicator_trip;
//	Eigen::SparseMatrix <double> Indicator(Problem.Variables_num, Problem.Variables_num);
	
	// Calculate non-constrained solution
	double ratio;
	Eigen::VectorXd bidded_total_aggregated_all = bidded_total_aggregated * Eigen::VectorXd::Constant(Market.num_zone, 1.);
	for(int price_iter = 0; price_iter < Market.price_intervals + 2; ++ price_iter){
		if(bidded_total_aggregated_all(price_iter) <= 0 && bidded_total_aggregated_all(price_iter + 1) >= 0){
			if(bidded_total_aggregated_all(price_iter + 1) > bidded_total_aggregated_all(price_iter)){
				ratio = bidded_total_aggregated_all(price_iter + 1) / (bidded_total_aggregated_all(price_iter + 1) - bidded_total_aggregated_all(price_iter));				
			}
			else{
				ratio = .5;
			}
			sol_temp.segment(Market.network.num_edges, Market.network.num_vertice) = ratio * bidded_total_aggregated.row(price_iter) + (1. - ratio) * bidded_total_aggregated.row(price_iter + 1);
			
			break;
		}
	}
	//std::cout << sol_temp.segment(Market.network.num_edges, Market.network.num_vertice).transpose() << "\n\n";
	
	// Increase a liitle bit of solution so it stays in the interior
	infeasbile_flag = 1;
	ratio = .5;
	Problem.Solution.orig_vector.segment(Market.network.num_edges, Market.network.num_vertice) = 2. * sol_temp.segment(Market.network.num_edges, Market.network.num_vertice);
	sol_temp = Problem.Solution.orig_vector;
	while(infeasbile_flag){
		Problem.Solution.orig_vector.segment(Market.network.num_edges, Market.network.num_vertice) = ratio * sol_temp.segment(Market.network.num_edges, Market.network.num_vertice);
		Problem.Solution.orig_vector.tail(Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(Problem.Solution.orig_vector.segment(Market.network.num_edges + 1, Market.network.num_vertice - 1));
		Problem.Solution.orig_vector.head(Market.network.num_edges) = (Problem.Constraint.eq_orig_matrix).topRightCorner(Market.network.num_edges, Market.network.num_vertice) * Problem.Solution.orig_vector.tail(Market.network.num_vertice);
		sol_temp = Problem.Solution.orig_vector;
		
		// Calculate initial objective function
		bar = 0.;
		infeasbile_flag = 0;
		for(int constraint_iter = 0; constraint_iter < Problem.Variables_num; ++ constraint_iter){
			if(constraint_iter < Market.network.num_edges || constraint_iter >= Market.network.num_edges + Market.network.num_vertice){
				if(sol_temp(constraint_iter) <= Problem.Boundary.ie_orig_matrix(constraint_iter, 0) || sol_temp(constraint_iter) >= Problem.Boundary.ie_orig_matrix(constraint_iter, 1)){
					infeasbile_flag = 1;
					break;
				}
				bar += std::log(sol_temp(constraint_iter) - Problem.Boundary.ie_orig_matrix(constraint_iter, 0)) + std::log(Problem.Boundary.ie_orig_matrix(constraint_iter, 1) - sol_temp(constraint_iter));
				bar -= 2 * std::log(Problem.Boundary.ie_orig_matrix(constraint_iter, 1) - Problem.Boundary.ie_orig_matrix(constraint_iter, 0));				
			}
		}
		
		// Check if increment is in the interior
		if(!infeasbile_flag){
			// Update objective
			obj = bar;
			break;
		}		
	}	
	
	// Main loop
	while(mu > .5 * eps){
	//while(mu > eps * .0001){
		loop_count = 0;
		break_flag = 0;
		while(!break_flag){
			break_flag = 1;
			loop_count += 1;
			//std::cout << "---------------------------------------------------------------------------\n";
			//std::cout << loop_count << "\n";
			//std::cout << "---------------------------------------------------------------------------\n";		
			
			grad = Eigen::VectorXd::Zero(Problem.Variables_num);
			for(int constraint_iter = 0; constraint_iter < Problem.Variables_num; ++ constraint_iter){
				if(constraint_iter < Market.network.num_edges || constraint_iter >= Market.network.num_edges + Market.network.num_vertice){
					// Calculate the error and examine whether boundary is breached
					grad(constraint_iter) = 1. / (sol_temp(constraint_iter) - Problem.Boundary.ie_orig_matrix(constraint_iter, 0)) - 1. / (sol_temp(constraint_iter) - Problem.Boundary.ie_orig_matrix(constraint_iter, 1));
					grad(constraint_iter) *= mu / (Problem.Boundary.ie_orig_matrix(constraint_iter, 1) - Problem.Boundary.ie_orig_matrix(constraint_iter, 0));				
				}
				else{
					// Update objective function and gradient
					Problem.Objective.orig_vector(constraint_iter) = -Market.bidded_price(price_ID(constraint_iter - Market.network.num_edges));
					Problem.Objective.orig_vector(constraint_iter) += price_penalty * (sol_temp(constraint_iter) < Problem.Boundary.ie_orig_matrix(constraint_iter, 0));
					Problem.Objective.orig_vector(constraint_iter) -= price_penalty * (sol_temp(constraint_iter) > Problem.Boundary.ie_orig_matrix(constraint_iter, 1));
					grad(constraint_iter) += (1. - mu) * Problem.Objective.orig_vector(constraint_iter);
				}	
			}
			
			// Find the projected gradient and update the solution
			grad -= Problem.Constraint.eq_orig_matrix.transpose() * Problem.Solver.qr.solve(Problem.Constraint.eq_orig_matrix * grad);
			if(grad.norm() > tol){
				grad /= grad.norm();
			}
			else{
				break;
			}
			//std::cout << grad.segment(Market.network.num_edges, Market.network.num_vertice).transpose() << "\n\n";
			std::cout << grad.segment(494, 12).transpose() << "\n\n";
			
			// Update the solution in a loop
	//		if(loop_count > 1){
	//			// Compute psuedo 2nd derivative of objective from prior information
	//		}
	//		else{
	//			dS = dS_max;
	//		}
			dS = dS_max;
			while(dS >= pow(mu, .5) * eps){
				sol_temp = Problem.Solution.orig_vector + grad * dS;
				
				// Update temporary price and utility of each node
				for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
					// Price
					if(sol_temp(Market.network.num_edges + zone_iter) < bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)){
						while(sol_temp(Market.network.num_edges + zone_iter) < bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)){				
							if(price_ID_temp(zone_iter) == 0){
								break;
							}				
							price_ID_temp(zone_iter) -= 1;				
						}
					}
					else if(sol_temp(Market.network.num_edges + zone_iter) > bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)){
						while(sol_temp(Market.network.num_edges + zone_iter) > bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)){
							if(price_ID_temp(zone_iter) == Market.price_intervals + 1){
								break;
							}				
							price_ID_temp(zone_iter) += 1;			
						}
					}
					
					// Utility
					if(sol_temp(Market.network.num_edges + zone_iter) > Problem.Boundary.ie_orig_matrix(Market.network.num_edges + zone_iter, 1)){
						utility_temp(zone_iter) = utility_aggregated(Market.price_intervals + 2, zone_iter);
						utility_temp(zone_iter) -= (sol_temp(Market.network.num_edges + zone_iter) - bidded_total_aggregated(Market.price_intervals + 2, zone_iter)) * Market.price_range_inflex(1);
					}
					else if(sol_temp(Market.network.num_edges + zone_iter) < Problem.Boundary.ie_orig_matrix(Market.network.num_edges + zone_iter, 0)){
						utility_temp(zone_iter) = utility_aggregated(0, zone_iter);
						utility_temp(zone_iter) -= (sol_temp(Market.network.num_edges + zone_iter) - bidded_total_aggregated(0, zone_iter)) * Market.price_range_inflex(0);
					}
					else{
						utility_temp(zone_iter) = (bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - sol_temp(Market.network.num_edges + zone_iter)) * utility_aggregated(price_ID_temp(zone_iter), zone_iter)
							+ (sol_temp(Market.network.num_edges + zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)) * utility_aggregated(price_ID_temp(zone_iter) + 1, zone_iter);
						if(bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter) != 0){
							utility_temp(zone_iter) /= bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter);
						}	
					}	
				}
				
				// Examine whether boundary is breached and update log barrier function
				bar_temp = 0.;
				infeasbile_flag = 0;
				//Power flow
				for(int constraint_iter = 0; constraint_iter < Market.network.num_edges; ++ constraint_iter){
					if(sol_temp(constraint_iter) <= Problem.Boundary.ie_orig_matrix(constraint_iter, 0) || sol_temp(constraint_iter) >= Problem.Boundary.ie_orig_matrix(constraint_iter, 1)){
						infeasbile_flag = 1;
						break;
					}
					bar_temp += std::log(sol_temp(constraint_iter) - Problem.Boundary.ie_orig_matrix(constraint_iter, 0)) + std::log(Problem.Boundary.ie_orig_matrix(constraint_iter, 1) - sol_temp(constraint_iter));
					bar_temp -= 2 * std::log(Problem.Boundary.ie_orig_matrix(constraint_iter, 1) - Problem.Boundary.ie_orig_matrix(constraint_iter, 0));
				}
				//Voltage
				for(int constraint_iter = Market.network.num_edges + Market.network.num_vertice; constraint_iter < Problem.Variables_num; ++ constraint_iter){
					if(sol_temp(constraint_iter) <= Problem.Boundary.ie_orig_matrix(constraint_iter, 0) || sol_temp(constraint_iter) >= Problem.Boundary.ie_orig_matrix(constraint_iter, 1)){
						infeasbile_flag = 1;
						break;
					}
					bar_temp += std::log(sol_temp(constraint_iter) - Problem.Boundary.ie_orig_matrix(constraint_iter, 0)) + std::log(Problem.Boundary.ie_orig_matrix(constraint_iter, 1) - sol_temp(constraint_iter));
					bar_temp -= 2 * std::log(Problem.Boundary.ie_orig_matrix(constraint_iter, 1) - Problem.Boundary.ie_orig_matrix(constraint_iter, 0));
				}				
				//std::cout << infeasbile_flag << "\n";
				//std::cout << mu << " " << dS << " " << infeasbile_flag << "\n";
				
				// Check current objective function
				obj_temp = (1. - mu) * utility_temp.sum() + mu * bar_temp;
				
				// Update direction if feasible
				if(!infeasbile_flag && obj_temp > obj){
					// Update solution
					sol_temp = Problem.Solution.orig_vector + grad * dS;
					Problem.Solution.orig_vector = sol_temp;
					price_ID = price_ID_temp;
					utility = utility_temp;
					obj = obj_temp;
					bar = bar_temp;
					break_flag = 0;
					break;
				}
				else{
					dS /= 2.;
					sol_temp = Problem.Solution.orig_vector;
					price_ID_temp = price_ID;
					utility_temp = utility;
					infeasbile_flag = 0;
				}			
			}
		}		
		// Update scale factor and objective
		mu /= 2.;
		obj = (1. - mu) * utility.sum() + mu * bar;
	}
	
	//std::cout << Problem.Solution.orig_vector.segment(Market.network.num_edges, Market.network.num_vertice).transpose() << "\n\n";
	std::cout << .5 * Problem.Solution.orig_vector.segment(Market.network.num_edges, Market.network.num_vertice).array().abs().sum() << " " << utility.sum() << "\n";
	std::cout << Problem.Solution.orig_vector.segment(Market.network.num_edges, Market.network.num_vertice).minCoeff() << " " << Problem.Solution.orig_vector.segment(Market.network.num_edges, Market.network.num_vertice).maxCoeff() << "\n";
	std::cout << Problem.Solution.orig_vector.tail(Market.network.num_vertice).minCoeff() << " " << Problem.Solution.orig_vector.tail(Market.network.num_vertice).maxCoeff() << "\n";
	std::cout << Problem.Solution.orig_vector.head(Market.network.num_edges).minCoeff() << " " << Problem.Solution.orig_vector.head(Market.network.num_edges).maxCoeff() << "\n\n";	
}


int main(){
	network_inform Power_network_inform;
	power_network_input_process(Power_network_inform, "../../power_network/");	
	
	market_inform TSO_Market;
	TSO_Market_Set_Test_2(TSO_Market, 1);
	LP_object TSO_Problem;
//	Flow_Based_Market_LP_Set(TSO_Market, TSO_Problem);
//	Flow_Based_Market_Optimization(0, TSO_Market, TSO_Problem);

	Flow_Based_Market_LP_Set(TSO_Market, TSO_Problem);
	Flow_Based_Market_Optimization_Test(0, TSO_Market, TSO_Problem);
}