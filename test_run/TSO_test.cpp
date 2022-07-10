// Source file for re-dispatch and tertiary control reserve market clearing of TSO in Norway
#include <iostream>
#include <iomanip>
#include <chrono>
#include "../basic/rw_csv.cpp"
#include "../basic/alglib/optimization.h"
#include "../power_network/power_network.h"
#include "../power_market/power_market.h"

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
	TSO_Market.num_zone = 20;
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
	TSO_Market.submitted_supply.leftCols(TSO_Market.num_zone / 2 - 1) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 1.);
	TSO_Market.submitted_demand.rightCols(TSO_Market.num_zone / 2 - 1) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 1.);
}

void Flow_Based_Market_LP_Set_Test(market_inform &Market, alglib::minlpstate &Problem){
	// -------------------------------------------------------------------------------
	// LP Solver initialization for flow-based market optimization
	// Warm-up once and reuse for the rest of time slices
	// Variables are sorted as {V}, {S}, {I}
	// -------------------------------------------------------------------------------

	// -------------------------------------------------------------------------------
	// Set matrix for general constraints
	// -------------------------------------------------------------------------------	
	// Construct node admittance matrix
	std::vector <Trip> Y_n_trip;
	Y_n_trip.reserve(2 * Market.network.num_edges + Market.network.num_vertice);
	Eigen::SparseMatrix <double, Eigen::RowMajor> Y_n(Market.network.num_vertice, Market.network.num_vertice);
	Eigen::VectorXd Y_n_diag = Eigen::VectorXd::Zero(Market.network.num_vertice);
	Eigen::VectorXpd Connection_num = Eigen::VectorXpd::Ones(Market.network.num_vertice);
	for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
		// Equality constraints of voltage - source / sink at the nodes, off-diagonal terms
		Y_n_trip.push_back(Trip(Market.network.incidence_matrix(edge_iter, 0), Market.network.incidence_matrix(edge_iter, 1), -Market.network.admittance_vector(edge_iter)));
		Y_n_trip.push_back(Trip(Market.network.incidence_matrix(edge_iter, 1), Market.network.incidence_matrix(edge_iter, 0), -Market.network.admittance_vector(edge_iter)));
		Connection_num(Market.network.incidence_matrix(edge_iter, 0)) += 1;
		Connection_num(Market.network.incidence_matrix(edge_iter, 1)) += 1;
		
		// Equality constraints of voltage - source / sink at the nodes, diagonal terms
		Y_n_diag(Market.network.incidence_matrix(edge_iter, 0)) += Market.network.admittance_vector(edge_iter);
		Y_n_diag(Market.network.incidence_matrix(edge_iter, 1)) += Market.network.admittance_vector(edge_iter);
	}
	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		// Equality constraints of voltage - source / sink at the nodes, summed diagonal terms
		Y_n_trip.push_back(Trip(node_iter, node_iter, Y_n_diag(node_iter)));
	}
	Y_n.setFromTriplets(Y_n_trip.begin(), Y_n_trip.end());

	// Generate sparse matrix for general (equality) constraints of voltage, power flow, and source / sink summation
	int constrant_num = Market.network.num_vertice + Market.network.num_edges + 1;
	int variable_num = 2 * Market.network.num_vertice + Market.network.num_edges;	
	Eigen::VectorXpd non_zero_num(constrant_num);
	non_zero_num << Connection_num + Eigen::VectorXpd::Ones(Market.network.num_vertice), Eigen::VectorXpd::Constant(Market.network.num_edges, 3), 1;
	alglib::integer_1d_array row_sizes_general;
	row_sizes_general.setcontent(non_zero_num.size(), non_zero_num.data());
	alglib::sparsematrix constraint_general;
	alglib::sparsecreatecrs(constrant_num, variable_num, row_sizes_general, constraint_general);
	// Rows for voltage - source / sink equalities
	for(int row_iter = 0; row_iter < Y_n.outerSize(); ++ row_iter){
		for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator inner_iter(Y_n, row_iter); inner_iter; ++ inner_iter){
			alglib::sparseset(constraint_general, inner_iter.row(), inner_iter.col(), inner_iter.value());
		}
		// Update the columns right to the node admittance block
		alglib::sparseset(constraint_general, row_iter, Market.network.num_vertice + row_iter, -1.);
	}
	// Rows for voltage - power flow equalities
	for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
		if(Market.network.incidence_matrix(edge_iter, 0) < Market.network.incidence_matrix(edge_iter, 1)){
			alglib::sparseset(constraint_general, Y_n.rows() + edge_iter, Market.network.incidence_matrix(edge_iter, 0), 1.);
			alglib::sparseset(constraint_general, Y_n.rows() + edge_iter, Market.network.incidence_matrix(edge_iter, 1), -1.);
		}
		else{
			alglib::sparseset(constraint_general, Y_n.rows() + edge_iter, Market.network.incidence_matrix(edge_iter, 1), -1.);
			alglib::sparseset(constraint_general, Y_n.rows() + edge_iter, Market.network.incidence_matrix(edge_iter, 0), 1.);						
		}
		alglib::sparseset(constraint_general, Y_n.rows() + edge_iter, 2 * Market.network.num_vertice + edge_iter, -1.);
	}
	
	// Voltage at reference node
	alglib::sparseset(constraint_general, non_zero_num.size() - 1, 0, 1.);
	
//	// Check if the sparse matrix is correct
//	std::cout << Y_n << "\n\n";
//	double value;
//	for(int row_iter = 0; row_iter < constrant_num; ++ row_iter){
//		for(int col_iter = 0; col_iter < variable_num; ++ col_iter){
//			value = sparseget(constraint_general, row_iter, col_iter);
//			std::cout << value << "\t";
//		}
//		std::cout << "\n";
//	}

	// -------------------------------------------------------------------------------
	// Set bounds for general and box constraints
	// -------------------------------------------------------------------------------
	Eigen::MatrixXd bound_general = Eigen::MatrixXd::Zero(constrant_num, 2);
	
	// Bounds of general constraints
	alglib::real_1d_array lb_general;
	alglib::real_1d_array ub_general;
	lb_general.setcontent(bound_general.rows(), bound_general.col(0).data());
	ub_general.setcontent(bound_general.rows(), bound_general.col(1).data());		
	
	// -------------------------------------------------------------------------------
	// Set the LP problem object
	// -------------------------------------------------------------------------------
	alglib::minlpcreate(variable_num, Problem);
	alglib::minlpsetlc2(Problem, constraint_general, lb_general, ub_general, constrant_num);
	//alglib::minlpsetalgoipm(Problem);
	alglib::minlpsetalgodss(Problem, 0);
}

void Flow_Based_Market_Optimization_Test(int tick, market_inform &Market, alglib::minlpstate &Problem){
	// -------------------------------------------------------------------------------
	// Initial optimization of each node as isolated trading zones
	// -------------------------------------------------------------------------------	
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

	// -------------------------------------------------------------------------------
	// Set scale of variables
	// -------------------------------------------------------------------------------
	int variable_num = 2 * Market.network.num_vertice + Market.network.num_edges;
	Eigen::VectorXd scale_vec(variable_num);
	scale_vec.head(Market.network.num_vertice) = Market.network.voltage_constraint.col(1) - Market.network.voltage_constraint.col(0);
	scale_vec.segment(Market.network.num_vertice, Market.network.num_vertice) = (bidded_total_aggregated.bottomRows(1) - bidded_total_aggregated.bottomRows(0)).transpose();
	scale_vec.tail(Market.network.num_edges) = Market.network.power_constraint.col(1) - Market.network.power_constraint.col(0);
	alglib::real_1d_array scale;
	scale.setcontent(scale_vec.size(), scale_vec.data());
	alglib::minlpsetscale(Problem, scale);

	// -------------------------------------------------------------------------------
	// Initilize source / sink bounds with isolated market clearing prices
	// -------------------------------------------------------------------------------
	int row_ID;
	Eigen::MatrixXd bound_box(variable_num, 2);
	bound_box.topRows(Market.network.num_vertice) = Market.network.voltage_constraint;
	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		row_ID = Market.network.num_vertice + node_iter;
		bound_box(row_ID, 0) = bidded_total_aggregated(price_ID(node_iter), node_iter);
		bound_box(row_ID, 1) = bidded_total_aggregated(price_ID(node_iter) + 1, node_iter);
	}
	bound_box.bottomRows(Market.network.num_edges) = Market.network.power_constraint;
	
	// Bounds of box constraints
	alglib::real_1d_array lb_box;
	alglib::real_1d_array ub_box;
	lb_box.setcontent(bound_box.rows(), bound_box.col(0).data());
	ub_box.setcontent(bound_box.rows(), bound_box.col(1).data());
	alglib::minlpsetbc(Problem, lb_box, ub_box);

	// -------------------------------------------------------------------------------
	// Initialize objective coefficients of variables
	// -------------------------------------------------------------------------------
	Eigen::VectorXd obj_vec = Eigen::VectorXd::Zero(variable_num);
	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		row_ID = Market.network.num_vertice + node_iter;
		obj_vec(row_ID) = Market.bidded_price(price_ID(node_iter));
	}
	alglib::real_1d_array obj_coeff;
	obj_coeff.setcontent(obj_vec.size(), obj_vec.data());
	alglib::minlpsetcost(Problem, obj_coeff);
	
	// -------------------------------------------------------------------------------
	// Main loop for iterative solver
	// -------------------------------------------------------------------------------	
	// Declare variables for the loop
	bool stop_flag = 0;
	double error;
	double tol = 1E-12;
	alglib::real_1d_array sol;
	alglib::real_1d_array sol_prev;
	alglib::minlpreport rep;
	
	// Solve the problem for the 1st time
	alglib::minlpoptimize(Problem);
	alglib::minlpresults(Problem, sol, rep);

	while(!stop_flag){
		stop_flag = 1;
		sol_prev = sol;
		
		// Update price at each node
		for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
			row_ID = Market.network.num_vertice + node_iter;
			if(rep.stats[row_ID] > 0){
				// Increase the price at the node if possible
				price_ID(node_iter) += (price_ID(node_iter) < Market.price_intervals + 1);
				while(bidded_total_aggregated(price_ID(node_iter), node_iter) == bidded_total_aggregated(price_ID(node_iter) + 1, node_iter)){
					if(price_ID(node_iter) == Market.price_intervals + 1){
						break;
					}
					price_ID(node_iter) += 1;
				}
				
				// Update box boundary for source / sink at the node
				minlpsetbci(Problem, row_ID, bidded_total_aggregated(price_ID(node_iter), node_iter), bidded_total_aggregated(price_ID(node_iter) + 1, node_iter));
				
				// Update objective coefficient at the node
				obj_coeff[row_ID] = Market.bidded_price(price_ID(node_iter));
				
				stop_flag = 0;
			}
			else if(rep.stats[row_ID] < 0){
				// Decrease the price at the node if possible
				price_ID(node_iter) -= (price_ID(node_iter) > 0);
				while(bidded_total_aggregated(price_ID(node_iter), node_iter) == bidded_total_aggregated(price_ID(node_iter) + 1, node_iter)){
					if(price_ID(node_iter) == 0){
						break;
					}
					price_ID(node_iter) -= 1;
				}
				
				// Update box boundary for source / sink at the node			
				minlpsetbci(Problem, row_ID, bidded_total_aggregated(price_ID(node_iter), node_iter), bidded_total_aggregated(price_ID(node_iter) + 1, node_iter));
				
				// Update objective coefficient at the node
				obj_coeff[row_ID] = Market.bidded_price(price_ID(node_iter));
				
				stop_flag = 0;								
			}
		}
		
		// Solve the problem
		alglib::minlpsetcost(Problem, obj_coeff);
		alglib::minlpoptimize(Problem);
		alglib::minlpresults(Problem, sol, rep);
		//std::cout << sol[3] << "\t" << sol[4] << "\t" << sol[5] << "\n\n";
		
		// Check if solution changes
		error = 0.;
		for(int var_iter; var_iter < variable_num; ++ var_iter){
			error += pow(sol[var_iter] - sol_prev[var_iter], 2.);
		}
		if(!stop_flag){
			stop_flag = (error < tol);
		}		
	}
	
	// Update confirmed prices and bids
	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		row_ID = Market.network.num_vertice + node_iter;
		std::cout << sol[row_ID] << "\t";
	}
	std::cout << "\n";
}

int main(){
	network_inform Power_network_inform;
	power_network_input_process(Power_network_inform, "../power_network/");	
	
	market_inform TSO_Market;
	TSO_Market_Set_Test_1(TSO_Market, 1);
	alglib::minlpstate TSO_Problem;

	auto start = std::chrono::high_resolution_clock::now();
	Flow_Based_Market_LP_Set_Test(TSO_Market, TSO_Problem);
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast <std::chrono::microseconds> (stop - start);
	std::cout << "Set time: " << duration.count() << " microseconds" << "\n\n";
	start = std::chrono::high_resolution_clock::now();
	Flow_Based_Market_Optimization_Test(0, TSO_Market, TSO_Problem);
	stop = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast <std::chrono::microseconds> (stop - start);
	std::cout << "Optimization time: " << duration.count() << " microseconds" << "\n\n";
}