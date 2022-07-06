// Source file for re-dispatch and tertiary control reserve market clearing of TSO in Norway
#include <iostream>
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
	TSO_Market.submitted_supply.leftCols(TSO_Market.num_zone / 2) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 1.);
	TSO_Market.submitted_demand.rightCols(TSO_Market.num_zone / 2) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 1.);
}

void Flow_Based_Market_Optimization_Test(int tick, market_inform &Market, LP_object &Problem){
	Eigen::MatrixXd bidded_supply = Market.submitted_supply;
	Eigen::MatrixXd bidded_demand = Market.submitted_demand;
	
	// Initial market clearing within each nodes
	Eigen::VectorXi price_ID(Market.num_zone);
	Market_clearing_nodal(tick, Market, price_ID, bidded_supply, bidded_demand);
	
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

	// Declare variables for the main loop
	double tol = pow(10., -12.);
	double eps = pow(10., -10.);
	double mu = 1. - eps;
	//double mu = 1 - pow(10., -6.);
	double dS = 100.;
	double obj = 0.;
	double obj_temp;
	Eigen::VectorXi price_ID_temp = price_ID;
	Eigen::VectorXd utility = Eigen::VectorXd::Zero(Market.num_zone);
	Eigen::VectorXd sol_temp = Problem.Solution.orig_vector;
	Eigen::VectorXd grad(Problem.Variables_num);
	Eigen::VectorXd error(Problem.Variables_num);
	std::vector <Trip> Indicator_trip;
	Eigen::SparseMatrix <double> Indicator(Problem.Variables_num, Problem.Variables_num);
	
	// Find the original gradient
	int loop_count = 0;
	while(loop_count < 50){
		std::cout << "---------------------------------------------------------------------------\n";
		std::cout << loop_count << "\n";
		std::cout << "---------------------------------------------------------------------------\n";		
		loop_count += 1;		
		
		Indicator_trip.clear();
		Indicator_trip.reserve(Problem.Variables_num);
		grad = Eigen::VectorXd::Zero(Problem.Variables_num);
		error = Eigen::VectorXd::Zero(Problem.Variables_num);
		sol_temp = Problem.Solution.orig_vector;
		for(int constraint_iter = 0; constraint_iter < Problem.Variables_num; ++ constraint_iter){
			// Update objective function and gradient
			if(constraint_iter >= Market.network.num_edges && constraint_iter < Market.network.num_edges + Market.network.num_vertice){
				Problem.Objective.orig_vector(constraint_iter) = -Market.bidded_price(price_ID(constraint_iter - Market.network.num_edges));
				grad(constraint_iter) = (1. - mu) * Problem.Objective.orig_vector(constraint_iter);
			}
			
			// Calculate the error and examine whether boundary is breached
			error(constraint_iter) = std::min(sol_temp(constraint_iter) - Problem.Boundary.ie_orig_matrix(constraint_iter, 0), 0.) + std::max(sol_temp(constraint_iter) - Problem.Boundary.ie_orig_matrix(constraint_iter, 1), 0.);
			if(error(constraint_iter) != 0){
				Indicator_trip.push_back(Trip(constraint_iter, constraint_iter, 1.));
				grad(constraint_iter) -= mu * error(constraint_iter);
			}
		}
		Indicator = Eigen::SparseMatrix <double> (Problem.Variables_num, Problem.Variables_num);
		Indicator.setFromTriplets(Indicator_trip.begin(), Indicator_trip.end());
		
		// Find the projected gradient and update the solution
		grad -= Problem.Constraint.eq_orig_matrix.transpose() * Problem.Solver.qr.solve(Problem.Constraint.eq_orig_matrix * grad);
		if(grad.norm() > tol){
			grad /= grad.norm();
		}
		else{
			break;
		}
		
		// Update the solution in a loop
		sol_temp = Problem.Solution.orig_vector + grad * dS;
		std::cout << Problem.Solution.orig_vector.transpose() << "\n";
		std::cout << grad.transpose() << "\n";
		std::cout << sol_temp.transpose() << "\n\n";
		
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
			std::cout << price_ID_temp.transpose() << "\n\n";
			
			// Utility
			if(sol_temp(Market.network.num_edges + zone_iter) > Problem.Boundary.ie_orig_matrix(Market.network.num_edges + zone_iter, 1)){
				utility(zone_iter) = utility_aggregated(Market.price_intervals + 2, zone_iter);
				utility(zone_iter) -= (sol_temp(Market.network.num_edges + zone_iter) - bidded_total_aggregated(Market.price_intervals + 2, zone_iter)) * Market.price_range_inflex(1);
			}
			else if(sol_temp(Market.network.num_edges + zone_iter) < Problem.Boundary.ie_orig_matrix(Market.network.num_edges + zone_iter, 0)){
				utility(zone_iter) = utility_aggregated(0, zone_iter);
				utility(zone_iter) -= (sol_temp(Market.network.num_edges + zone_iter) - bidded_total_aggregated(0, zone_iter)) * Market.price_range_inflex(0);
			}
			else{
				utility(zone_iter) = (bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - sol_temp(Market.network.num_edges + zone_iter)) * utility_aggregated(price_ID_temp(zone_iter), zone_iter)
					+ (sol_temp(Market.network.num_edges + zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)) * utility_aggregated(price_ID_temp(zone_iter) + 1, zone_iter);
				if(bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter) != 0){
					utility(zone_iter) /= bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter);
				}	
			}	
		}
		
		// Update errors
		for(int constraint_iter = 0; constraint_iter < Problem.Variables_num; ++ constraint_iter){
			// Calculate the error and examine whether boundary is breached
			error(constraint_iter) = std::min(sol_temp(constraint_iter) - Problem.Boundary.ie_orig_matrix(constraint_iter, 0), 0.) + std::max(sol_temp(constraint_iter) - Problem.Boundary.ie_orig_matrix(constraint_iter, 1), 0.);
		}	
		
		// Check if objective actually improves
		obj_temp = (1. - mu) * utility.sum() - mu * error.squaredNorm();
		//std::cout << utility.transpose() << "\n\n";
		//std::cout << sol_temp.transpose() << "\n\n";
		//std::cout << obj_temp << "\n";
		//std::cout << obj << "\n\n";
		if(obj_temp >= obj){
			// Update solution, price, and objective if objective actually improves
			Problem.Solution.orig_vector = sol_temp;
			price_ID = price_ID_temp;
			obj = obj_temp;
		}
		else{
			// Return to original values
			dS /= 2.;
			sol_temp = Problem.Solution.orig_vector;
			price_ID_temp = price_ID;
		}	
		//std::cout << Problem.Solution.orig_vector.transpose() << "\n\n";
	}
}

void Flow_Based_Market_Optimization_Test_2(int tick, market_inform &Market, LP_object &Problem){
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
	
	// Declare variables for the main loop
	bool break_flag = 0;
	double tol = pow(10., -12.);
	double eps = pow(10., -10.);	
	//double mu = 1. - eps;
	double mu = 0.;
	double dV = .001;
	double deficit_penalty = 0.;
	double error;
	double obj = 0.;
	Eigen::Vector2d obj_temp;
	Eigen::VectorXi price_ID_temp = price_ID;
	Eigen::VectorXd voltage = Eigen::VectorXd::Zero(Market.num_zone);
	Eigen::VectorXd voltage_temp = voltage;
	Eigen::VectorXd quan = Eigen::VectorXd::Zero(Market.num_zone);	
	Eigen::VectorXd quan_temp = quan;	
	Eigen::VectorXd flow = Eigen::VectorXd::Zero(Market.network.num_edges);
	Eigen::VectorXd flow_temp = flow;
	Eigen::VectorXd utility = Eigen::VectorXd::Zero(Market.num_zone);
	Eigen::VectorXd utility_temp = utility;
	Eigen::VectorXd grad(Market.num_zone);
	Eigen::MatrixXd Y_n = (Problem.Constraint.eq_orig_matrix).block(Market.network.num_edges, Market.network.num_edges + Market.network.num_vertice, Market.network.num_vertice, Market.network.num_vertice);
	
	for(int zone_iter = 500; zone_iter < 501; ++ zone_iter){
		for(int dir_iter = 0; dir_iter < 2; ++ dir_iter){
			// Change of voltage, power flow, power source / sink
			voltage_temp(zone_iter) += (double(dir_iter) - .5) * dV;
			quan_temp = Y_n * voltage_temp;
			flow_temp = (Problem.Constraint.eq_orig_matrix).topRightCorner(Market.network.num_edges, Market.network.num_vertice) * voltage_temp;
			std::cout << quan_temp.transpose() << "\n\n";
			
			for(int zone_iter_2 = 0; zone_iter_2 < Market.num_zone; ++ zone_iter_2){
				// Determine the price of each node
				if(quan_temp(zone_iter_2) < bidded_total_aggregated(price_ID_temp(zone_iter_2), zone_iter_2)){
					while(quan_temp(zone_iter_2) < bidded_total_aggregated(price_ID_temp(zone_iter_2), zone_iter_2)){
						if(price_ID_temp(zone_iter_2) == 0){
							break;
						}
						price_ID_temp(zone_iter_2) -= 1;
					}
				}
				else if(quan_temp(zone_iter_2) > bidded_total_aggregated(price_ID_temp(zone_iter_2) + 1, zone_iter_2)){
					while(quan_temp(zone_iter_2) > bidded_total_aggregated(price_ID_temp(zone_iter_2) + 1, zone_iter_2)){
						if(price_ID_temp(zone_iter_2) == Market.price_intervals + 1){
							break;
						}
						price_ID_temp(zone_iter_2) += 1;
					}
				}
				
				// Determine the utility of each node
				if(quan_temp(zone_iter_2) < bidded_total_aggregated(0, zone_iter_2)){
					utility_temp(zone_iter_2) = utility_aggregated(0, zone_iter_2);
					utility_temp(zone_iter_2) -= (quan_temp(zone_iter_2) - bidded_total_aggregated(0, zone_iter_2)) * Market.price_range_inflex(0);
					//std::cout << "Case 1 " << utility(zone_iter_2) << "\n";
				}
				else if(quan_temp(zone_iter_2) > bidded_total_aggregated(Market.price_intervals + 2, zone_iter_2))	{
					utility_temp(zone_iter_2) = utility_aggregated(Market.price_intervals + 2, zone_iter_2);
					utility_temp(zone_iter_2) -= (quan_temp(zone_iter_2) - bidded_total_aggregated(Market.price_intervals + 2, zone_iter_2)) * Market.price_range_inflex(1);
					//std::cout << "Case 2 " << utility(zone_iter_2) << "\n";
				}
				else{
					if(bidded_total_aggregated(price_ID_temp(zone_iter_2) + 1, zone_iter_2) - bidded_total_aggregated(price_ID_temp(zone_iter_2), zone_iter_2) != 0){
						utility_temp(zone_iter_2) = (bidded_total_aggregated(price_ID_temp(zone_iter_2) + 1, zone_iter_2) - quan_temp(zone_iter_2)) * utility_aggregated(price_ID_temp(zone_iter_2), zone_iter_2)
							+ (quan_temp(zone_iter_2) - bidded_total_aggregated(price_ID_temp(zone_iter_2), zone_iter_2)) * utility_aggregated(price_ID_temp(zone_iter_2) + 1, zone_iter_2);
						utility_temp(zone_iter_2) /= bidded_total_aggregated(price_ID_temp(zone_iter_2) + 1, zone_iter_2) - bidded_total_aggregated(price_ID_temp(zone_iter_2), zone_iter_2);	
					}
					else{
						utility_temp(zone_iter_2) = utility_aggregated(price_ID_temp(zone_iter_2), zone_iter_2);
					}
					//std::cout << "Case 3 " << utility(zone_iter_2) << "\n";
				}			
			}
//			std::cout << price_ID_temp.transpose() << "\n\n";
			std::cout << utility_temp.transpose() << "\n\n";
			
			// Determine the error
			error = 0.;
			// Power flow errors
			for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
				error += pow(std::min(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 0), 0.) + std::max(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 1), 0.), 2.);
			}
			// Power errors
			for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
				error += pow(std::min(quan_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + node_iter, 0), 0.) + std::max(quan_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + node_iter, 1), 0.), 2.);
			}
			// Voltage errors
			for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
				error += pow(std::min(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 0), 0.) + std::max(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 1), 0.), 2.);
			}
			
			// Calculate objective and update solution when improved
			obj_temp(dir_iter) = (1. - mu) * utility_temp.sum() - mu * error;
			
			// Return to origin
			voltage_temp = voltage;
			price_ID_temp = price_ID;
			utility_temp = utility;
//			if(obj_temp(dir_iter) >= obj){
//				//std::cout << "Objective updated" << "\n\n";
//				voltage = voltage_temp;
//				flow = flow_temp;
//				quan = quan_temp;
//				price_ID = price_ID_temp;
//				continue;
//			}
//			else{
//				//std::cout << "Objective not updated" << "\n\n";
//				voltage_temp = voltage;
//				price_ID_temp = price_ID;
//			}
		}
		//std::cout << obj_temp.transpose() << "\n\n";
		
		
		//grad(zone_iter) = obj_temp(1) - obj_temp(0);
	}
	
	//grad = 
//	std::cout << voltage.transpose() << "\n\n";
}

void Flow_Based_Market_Optimization_Test_3(int tick, market_inform &Market, LP_object &Problem){
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
	
	// Declare variables for the main loop
	bool divide_flag;
	int node_ref_ID;
	int flow_dir;
	double tol = pow(10., -12.);
	double eps = pow(10., -10.);
	double mu = 1. - eps;
	//double mu = 0.;
	double dS = 100.;
	double error;
	double obj = 0.;
	double obj_temp;
	Eigen::Vector2i obj_temp_vec;
	Eigen::VectorXi price_ID_temp = price_ID;
	Eigen::VectorXd voltage = Eigen::VectorXd::Zero(Market.num_zone);
	Eigen::VectorXd voltage_temp = voltage;
	Eigen::VectorXd quan = Eigen::VectorXd::Zero(Market.num_zone);	
	Eigen::VectorXd quan_temp = quan;	
	Eigen::VectorXd flow = Eigen::VectorXd::Zero(Market.network.num_edges);
	Eigen::VectorXd flow_temp = flow;
	Eigen::VectorXd utility = Eigen::VectorXd::Zero(Market.num_zone);
	Eigen::VectorXd utility_temp = utility;
	Eigen::VectorXd grad = Eigen::VectorXd::Zero(Market.num_zone);
	
	int loop_count = 0;
	while(loop_count < 100){
		loop_count += 1;
		
		divide_flag = 1;
		node_ref_ID = std::rand() % Market.num_zone;
		for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
			if(zone_iter == node_ref_ID){
				continue;
			}			
			
			for(int dir_iter = 0; dir_iter < 2; ++ dir_iter){
				// Update quantity after small increase / decrease
				quan_temp(zone_iter) += 2 * (dir_iter - .5) * dS;
				quan_temp(node_ref_ID) -= 2 * (dir_iter - .5) * dS;
				
				// Update price after small increase / decrease
				if(quan_temp(zone_iter) < bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)){
					while(quan_temp(zone_iter) < bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)){
						if(price_ID_temp(zone_iter) == 0){
							break;
						}
						price_ID_temp(zone_iter) -= 1;
					}
				}
				else if(quan_temp(zone_iter) > bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter)){
					while(quan_temp(zone_iter) > bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter)){
						if(price_ID_temp(zone_iter) == Market.price_intervals + 1){
							break;
						}
						price_ID_temp(zone_iter) += 1;
					}
				}
				if(quan_temp(node_ref_ID) < bidded_total_aggregated(price_ID_temp(node_ref_ID), node_ref_ID)){
					while(quan_temp(node_ref_ID) < bidded_total_aggregated(price_ID_temp(node_ref_ID), node_ref_ID)){
						if(price_ID_temp(node_ref_ID) == 0){
							break;
						}
						price_ID_temp(node_ref_ID) -= 1;
					}
				}
				else if(quan_temp(node_ref_ID) > bidded_total_aggregated(price_ID_temp(node_ref_ID) + 1, node_ref_ID)){
					while(quan_temp(node_ref_ID) > bidded_total_aggregated(price_ID_temp(node_ref_ID) + 1, node_ref_ID)){
						if(price_ID_temp(node_ref_ID) == Market.price_intervals + 1){
							break;
						}
						price_ID_temp(node_ref_ID) += 1;
					}
				}				
				
				// Update utility after small increase / decrease
				if(quan_temp(zone_iter) < bidded_total_aggregated(0, zone_iter)){
					utility_temp(zone_iter) = utility_aggregated(0, zone_iter);
					utility_temp(zone_iter) -= (quan_temp(zone_iter) - bidded_total_aggregated(0, zone_iter)) * Market.price_range_inflex(0);
				}
				else if(quan_temp(zone_iter) > bidded_total_aggregated(Market.price_intervals + 2, zone_iter))	{
					utility_temp(zone_iter) = utility_aggregated(Market.price_intervals + 2, zone_iter);
					utility_temp(zone_iter) -= (quan_temp(zone_iter) - bidded_total_aggregated(Market.price_intervals + 2, zone_iter)) * Market.price_range_inflex(1);
				}
				else{
					if(bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter) != 0){
						utility_temp(zone_iter) = (bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - quan_temp(zone_iter)) * utility_aggregated(price_ID_temp(zone_iter), zone_iter)
							+ (quan_temp(zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)) * utility_aggregated(price_ID_temp(zone_iter) + 1, zone_iter);
						utility_temp(zone_iter) /= bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter);	
					}
					else{
						utility_temp(zone_iter) = utility_aggregated(price_ID_temp(zone_iter), zone_iter);
					}
				}
				if(quan_temp(node_ref_ID) < bidded_total_aggregated(0, node_ref_ID)){
					utility_temp(node_ref_ID) = utility_aggregated(0, node_ref_ID);
					utility_temp(node_ref_ID) -= (quan_temp(node_ref_ID) - bidded_total_aggregated(0, node_ref_ID)) * Market.price_range_inflex(0);
				}
				else if(quan_temp(node_ref_ID) > bidded_total_aggregated(Market.price_intervals + 2, node_ref_ID))	{
					utility_temp(node_ref_ID) = utility_aggregated(Market.price_intervals + 2, node_ref_ID);
					utility_temp(node_ref_ID) -= (quan_temp(node_ref_ID) - bidded_total_aggregated(Market.price_intervals + 2, node_ref_ID)) * Market.price_range_inflex(1);
				}
				else{
					if(bidded_total_aggregated(price_ID_temp(node_ref_ID) + 1, node_ref_ID) - bidded_total_aggregated(price_ID_temp(node_ref_ID), node_ref_ID) != 0){
						utility_temp(node_ref_ID) = (bidded_total_aggregated(price_ID_temp(node_ref_ID) + 1, node_ref_ID) - quan_temp(node_ref_ID)) * utility_aggregated(price_ID_temp(node_ref_ID), node_ref_ID)
							+ (quan_temp(node_ref_ID) - bidded_total_aggregated(price_ID_temp(node_ref_ID), node_ref_ID)) * utility_aggregated(price_ID_temp(node_ref_ID) + 1, node_ref_ID);
						utility_temp(node_ref_ID) /= bidded_total_aggregated(price_ID_temp(node_ref_ID) + 1, node_ref_ID) - bidded_total_aggregated(price_ID_temp(node_ref_ID), node_ref_ID);	
					}
					else{
						utility_temp(node_ref_ID) = utility_aggregated(price_ID_temp(node_ref_ID), node_ref_ID);
					}
				}
				
				// Update voltage and power flow after small increase / decrease
				voltage_temp.tail(Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(quan_temp.tail(Market.network.num_vertice - 1));
				flow_temp = (Problem.Constraint.eq_orig_matrix).topRightCorner(Market.network.num_edges, Market.network.num_vertice) * voltage_temp;		
		
				// Update errors
				error = 0.;
				// Power flow errors
				for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
					error += pow(std::min(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 0), 0.) + std::max(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 1), 0.), 2.);
				}
				// Power errors
				for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
					error += pow(std::min(quan_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + node_iter, 0), 0.) + std::max(quan_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + node_iter, 1), 0.), 2.);
				}
				// Voltage errors
				for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
					error += pow(std::min(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 0), 0.) + std::max(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 1), 0.), 2.);
				}
				
				// Calculate objective and update solution when improved
				obj_temp_vec(dir_iter) = (1. - mu) * utility_temp.sum() - mu * error;
				
				// Return to origin
				quan_temp(zone_iter) = quan(zone_iter);
				quan_temp(node_ref_ID) = quan(node_ref_ID);
				price_ID_temp(zone_iter) = price_ID(zone_iter);
				price_ID_temp(node_ref_ID) = price_ID(node_ref_ID);
				utility_temp(zone_iter) = utility(zone_iter);
				utility_temp(node_ref_ID) = utility(node_ref_ID);																					
			}
			
			// Update gradient at this direction
			grad(zone_iter) = obj_temp_vec(1) - obj_temp_vec(0);
			//grad(zone_iter) = obj_temp_vec.maxCoeff();
			//grad(zone_iter) = std::max(grad(zone_iter), obj) - obj;	
		}
		
		// Normalized gradient
		grad(node_ref_ID) = -grad.sum();
		if(grad.norm() != 0.){
			grad /= grad.norm();
		}
		
		// Update solution and objective
		quan_temp += grad * dS;
		for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
			// Update price after small increase / decrease
			if(quan_temp(zone_iter) < bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)){
				while(quan_temp(zone_iter) < bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)){
					if(price_ID_temp(zone_iter) == 0){
						break;
					}
					price_ID_temp(zone_iter) -= 1;
				}
			}
			else if(quan_temp(zone_iter) > bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter)){
				while(quan_temp(zone_iter) > bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter)){
					if(price_ID_temp(zone_iter) == Market.price_intervals + 1){
						break;
					}
					price_ID_temp(zone_iter) += 1;
				}
			}
			
			// Update utility after small increase / decrease
			if(quan_temp(zone_iter) < bidded_total_aggregated(0, zone_iter)){
				utility_temp(zone_iter) = utility_aggregated(0, zone_iter);
				utility_temp(zone_iter) -= (quan_temp(zone_iter) - bidded_total_aggregated(0, zone_iter)) * Market.price_range_inflex(0);
			}
			else if(quan_temp(zone_iter) > bidded_total_aggregated(Market.price_intervals + 2, zone_iter))	{
				utility_temp(zone_iter) = utility_aggregated(Market.price_intervals + 2, zone_iter);
				utility_temp(zone_iter) -= (quan_temp(zone_iter) - bidded_total_aggregated(Market.price_intervals + 2, zone_iter)) * Market.price_range_inflex(1);
			}
			else{
				if(bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter) != 0){
					utility_temp(zone_iter) = (bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - quan_temp(zone_iter)) * utility_aggregated(price_ID_temp(zone_iter), zone_iter)
						+ (quan_temp(zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)) * utility_aggregated(price_ID_temp(zone_iter) + 1, zone_iter);
					utility_temp(zone_iter) /= bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter);	
				}
				else{
					utility_temp(zone_iter) = utility_aggregated(price_ID_temp(zone_iter), zone_iter);
				}
			}							
		}
		
		// Update voltage and power flow after small increase / decrease
		voltage_temp.tail(Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(quan_temp.tail(Market.network.num_vertice - 1));
		flow_temp = (Problem.Constraint.eq_orig_matrix).topRightCorner(Market.network.num_edges, Market.network.num_vertice) * voltage_temp;		

		// Update errors
		error = 0.;
		// Power flow errors
		for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
			error += pow(std::min(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 0), 0.) + std::max(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 1), 0.), 2.);
		}
		// Power errors
		for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
			error += pow(std::min(quan_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + node_iter, 0), 0.) + std::max(quan_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + node_iter, 1), 0.), 2.);
		}
		// Voltage errors
		for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
			error += pow(std::min(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 0), 0.) + std::max(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 1), 0.), 2.);
		}
		
		// Calculate objective and update solution when improved
		obj_temp = (1. - mu) * utility_temp.sum() - mu * error;
			
		// Check if objective should be updated
		if(obj_temp > obj){
			// Update if objective has been improved
			quan = quan_temp;
			voltage = voltage_temp;
			flow = flow_temp;
			price_ID = price_ID_temp;
			utility = utility_temp;
			obj = obj_temp;	
			divide_flag = 0;
			//break;
		}
		else{
			//std::cout << (1. - mu) * utility_temp.sum() << " " << mu * error << " " << "\n\n";
			// Return to origin
			dS /= 2.;
			quan_temp = quan;
			price_ID_temp = price_ID;
			utility_temp = utility;			
		}
	}
	
	std::cout << "\n" << obj << "\n\n";
	//std::cout << "\n" << grad.transpose() << "\n\n";
	std::cout << "\n" << quan.transpose() << "\n\n";
	//std::cout << "\n" << voltage.transpose() << "\n\n";	
	//std::cout << "\n" << flow.transpose() << "\n\n";
}

void Flow_Based_Market_Optimization_Test_4(int tick, market_inform &Market, LP_object &Problem){
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
	
	// Declare variables for the main loop
	bool divide_flag;
	int node_ref_ID;
	int zone_iter;
	double tol = pow(10., -12.);
	double eps = pow(10., -10.);
	double mu = 1 - eps;
	double dS = bidded_total_aggregated.lpNorm <Eigen::Infinity> ();
	double error;
	double obj = 0;
	double obj_temp;
	Eigen::VectorXi price_ID_temp = price_ID;
	Eigen::VectorXd voltage = Eigen::VectorXd::Zero(Market.num_zone);
	Eigen::VectorXd voltage_temp = voltage;
	Eigen::VectorXd quan = Eigen::VectorXd::Zero(Market.num_zone);
//	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
//		quan(zone_iter) = 100. * (zone_iter < Market.num_zone / 2) - 100. * (zone_iter >= Market.num_zone / 2);
//	}	
	Eigen::VectorXd quan_temp = quan;	
	Eigen::VectorXd flow = Eigen::VectorXd::Zero(Market.network.num_edges);
	Eigen::VectorXd flow_temp = flow;
	Eigen::VectorXd utility = Eigen::VectorXd::Zero(Market.num_zone);
	Eigen::VectorXd utility_temp = utility;
	Eigen::VectorXd grad = Eigen::VectorXd::Zero(Market.num_zone);
	
	int loop_count = 0;
	//while(loop_count < 1){
	while(dS > .01 * bidded_total_aggregated.lpNorm <Eigen::Infinity> ()){
		loop_count += 1;
		divide_flag = 1;
		//for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
		for(int rand_draw = 0; rand_draw < 5000; ++ rand_draw){
			zone_iter = std::rand() % Market.num_zone;
			node_ref_ID = std::rand() % Market.num_zone;
			if(zone_iter == node_ref_ID){
				continue;
			}			
			
			for(int dir_iter = 0; dir_iter < 2; ++ dir_iter){
				//Check if direction is plausible
//				if(price_ID_temp(zone_iter) < price_ID_temp(zone_iter) && dir_iter == 0){
//					continue;
//				}
//				else if(price_ID_temp(zone_iter) > price_ID_temp(zone_iter) && dir_iter == 1){
//					continue;
//				}
				
				// Update quantity after small increase / decrease
				quan_temp(zone_iter) += 2 * (dir_iter - .5) * dS;
				quan_temp(node_ref_ID) -= 2 * (dir_iter - .5) * dS;
				
				// Update price after small increase / decrease
				if(quan_temp(zone_iter) < bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)){
					while(quan_temp(zone_iter) < bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)){
						if(price_ID_temp(zone_iter) == 0){
							break;
						}
						price_ID_temp(zone_iter) -= 1;
					}
				}
				else if(quan_temp(zone_iter) > bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter)){
					while(quan_temp(zone_iter) > bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter)){
						if(price_ID_temp(zone_iter) == Market.price_intervals + 1){
							break;
						}
						price_ID_temp(zone_iter) += 1;
					}
				}
				if(quan_temp(node_ref_ID) < bidded_total_aggregated(price_ID_temp(node_ref_ID), node_ref_ID)){
					while(quan_temp(node_ref_ID) < bidded_total_aggregated(price_ID_temp(node_ref_ID), node_ref_ID)){
						if(price_ID_temp(node_ref_ID) == 0){
							break;
						}
						price_ID_temp(node_ref_ID) -= 1;
					}
				}
				else if(quan_temp(node_ref_ID) > bidded_total_aggregated(price_ID_temp(node_ref_ID) + 1, node_ref_ID)){
					while(quan_temp(node_ref_ID) > bidded_total_aggregated(price_ID_temp(node_ref_ID) + 1, node_ref_ID)){
						if(price_ID_temp(node_ref_ID) == Market.price_intervals + 1){
							break;
						}
						price_ID_temp(node_ref_ID) += 1;
					}
				}				
				
				// Update utility after small increase / decrease
				if(quan_temp(zone_iter) < bidded_total_aggregated(0, zone_iter)){
					utility_temp(zone_iter) = utility_aggregated(0, zone_iter);
					utility_temp(zone_iter) -= (quan_temp(zone_iter) - bidded_total_aggregated(0, zone_iter)) * Market.price_range_inflex(0);
				}
				else if(quan_temp(zone_iter) > bidded_total_aggregated(Market.price_intervals + 2, zone_iter))	{
					utility_temp(zone_iter) = utility_aggregated(Market.price_intervals + 2, zone_iter);
					utility_temp(zone_iter) -= (quan_temp(zone_iter) - bidded_total_aggregated(Market.price_intervals + 2, zone_iter)) * Market.price_range_inflex(1);
				}
				else{
					if(bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter) != 0){
						utility_temp(zone_iter) = (bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - quan_temp(zone_iter)) * utility_aggregated(price_ID_temp(zone_iter), zone_iter)
							+ (quan_temp(zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)) * utility_aggregated(price_ID_temp(zone_iter) + 1, zone_iter);
						utility_temp(zone_iter) /= bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter);	
					}
					else{
						utility_temp(zone_iter) = utility_aggregated(price_ID_temp(zone_iter), zone_iter);
					}
				}
				if(quan_temp(node_ref_ID) < bidded_total_aggregated(0, node_ref_ID)){
					utility_temp(node_ref_ID) = utility_aggregated(0, node_ref_ID);
					utility_temp(node_ref_ID) -= (quan_temp(node_ref_ID) - bidded_total_aggregated(0, node_ref_ID)) * Market.price_range_inflex(0);
				}
				else if(quan_temp(node_ref_ID) > bidded_total_aggregated(Market.price_intervals + 2, node_ref_ID))	{
					utility_temp(node_ref_ID) = utility_aggregated(Market.price_intervals + 2, node_ref_ID);
					utility_temp(node_ref_ID) -= (quan_temp(node_ref_ID) - bidded_total_aggregated(Market.price_intervals + 2, node_ref_ID)) * Market.price_range_inflex(1);
				}
				else{
					if(bidded_total_aggregated(price_ID_temp(node_ref_ID) + 1, node_ref_ID) - bidded_total_aggregated(price_ID_temp(node_ref_ID), node_ref_ID) != 0){
						utility_temp(node_ref_ID) = (bidded_total_aggregated(price_ID_temp(node_ref_ID) + 1, node_ref_ID) - quan_temp(node_ref_ID)) * utility_aggregated(price_ID_temp(node_ref_ID), node_ref_ID)
							+ (quan_temp(node_ref_ID) - bidded_total_aggregated(price_ID_temp(node_ref_ID), node_ref_ID)) * utility_aggregated(price_ID_temp(node_ref_ID) + 1, node_ref_ID);
						utility_temp(node_ref_ID) /= bidded_total_aggregated(price_ID_temp(node_ref_ID) + 1, node_ref_ID) - bidded_total_aggregated(price_ID_temp(node_ref_ID), node_ref_ID);	
					}
					else{
						utility_temp(node_ref_ID) = utility_aggregated(price_ID_temp(node_ref_ID), node_ref_ID);
					}
				}
				
				// Update voltage and power flow after small increase / decrease
				voltage_temp.tail(Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(quan_temp.tail(Market.network.num_vertice - 1));
				flow_temp = (Problem.Constraint.eq_orig_matrix).topRightCorner(Market.network.num_edges, Market.network.num_vertice) * voltage_temp;		
		
				// Update errors
				error = 0.;
				// Power flow errors
				for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
					error += pow(std::min(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 0), 0.) + std::max(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 1), 0.), 2.);
				}
				// Power errors
				for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
					error += pow(std::min(quan_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + node_iter, 0), 0.) + std::max(quan_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + node_iter, 1), 0.), 2.);
				}
				// Voltage errors
				for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
					error += pow(std::min(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 0), 0.) + std::max(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 1), 0.), 2.);
				}
				
				// Calculate objective 
				obj_temp = (1. - mu) * utility_temp.sum() - mu * error;
				//std::cout << zone_iter << " " << node_ref_ID << " " << obj_temp << "\n";
				
				// Check if objective should be updated
				if(obj_temp > obj){
					// Update if objective has been improved
					quan(zone_iter) = quan_temp(zone_iter);
					quan(node_ref_ID) = quan_temp(node_ref_ID);
					voltage = voltage_temp;
					flow = flow_temp;
					price_ID(zone_iter) = price_ID_temp(zone_iter);
					price_ID(node_ref_ID) = price_ID_temp(node_ref_ID);
					utility(zone_iter) = utility_temp(zone_iter);
					utility(node_ref_ID) = utility_temp(node_ref_ID);
					obj = obj_temp;	
					divide_flag = 0;
					break;
				}
				else{
					// Return to origin
					quan_temp(zone_iter) = quan(zone_iter);
					quan_temp(node_ref_ID) = quan(node_ref_ID);
					price_ID_temp(zone_iter) = price_ID(zone_iter);
					price_ID_temp(node_ref_ID) = price_ID(node_ref_ID);
					utility_temp(zone_iter) = utility(zone_iter);
					utility_temp(node_ref_ID) = utility(node_ref_ID);					
				}																					
			}	
		}		
		
		// Check if increment should be smaller
		if(divide_flag){
			dS /= 2.;
			//std::cout << dS << "\n";
		}
	}

	//std::cout << "\n" << grad.transpose() << "\n\n";
	//std::cout << "\n" << obj << "\n\n";
	std::cout << "\n" << quan.transpose() << "\n\n";
	//std::cout << "\n" << voltage.transpose() << "\n\n";
	//std::cout << "\n" << flow.transpose() << "\n\n";			
}

int main(){
	network_inform Power_network_inform;
	power_network_input_process(Power_network_inform, "../../power_network/");	
	
	market_inform TSO_Market;
	TSO_Market_Set_Test_2(TSO_Market, 1);
	LP_object TSO_Problem;
//	Flow_Based_Market_LP_Set(TSO_Market, TSO_Problem);
//	Flow_Based_Market_Optimization_Test_3(0, TSO_Market, TSO_Problem);

	Flow_Based_Market_LP_Set(TSO_Market, TSO_Problem);
	Flow_Based_Market_Optimization_Test_4(0, TSO_Market, TSO_Problem);
}