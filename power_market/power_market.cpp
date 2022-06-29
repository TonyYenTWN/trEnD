// Source file for the complete procedure of the power market
#include "power_market.h"

// ------------------------------------------------------------------------------------------------
// Generic functions for all market types
// ------------------------------------------------------------------------------------------------
void Market_Initialization(market_inform &Market){
	// Initialization of process variables
	// Should re-initialize for every time slice
	Market.submitted_supply = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Market.submitted_demand = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
}

void Market_clearing_nodal(int tick, market_inform &Market, Eigen::VectorXi &default_price_ID, Eigen::MatrixXd &bidded_supply, Eigen::MatrixXd &bidded_demand){
	Eigen::VectorXi price_demand_ID = (Market.price_intervals + 1) * Eigen::VectorXi::Ones(Market.num_zone);
	Eigen::VectorXi price_supply_ID = Eigen::VectorXi::Zero(Market.num_zone);
	
	double trade_quantity;
	//#pragma omp parallel for
	for(int zone_ID = 0; zone_ID < Market.num_zone; ++ zone_ID){
		while(price_demand_ID(zone_ID) > price_supply_ID(zone_ID)){
			// Check if there are demand bids at current price interval
			while(bidded_demand(price_demand_ID(zone_ID), zone_ID) == 0){
				if(price_demand_ID(zone_ID) > 0){
					price_demand_ID(zone_ID) -= 1;
				}
				else{
					// No available buyer left to buy electricity
					default_price_ID(zone_ID) = price_supply_ID(zone_ID);
					break;
				}
			}
			
			// Check if there are supply bids at current price interval
			while(bidded_supply(price_supply_ID(zone_ID), zone_ID) == 0){
				if(price_supply_ID(zone_ID) < Market.bidded_price.size() - 1){
					price_supply_ID(zone_ID) += 1;
				}
				else{
					// No available seller left to sell electricity
					default_price_ID(zone_ID) = price_demand_ID(zone_ID);
					break;				
				}			
			}
			
			if(price_demand_ID(zone_ID) > price_supply_ID(zone_ID)){
				trade_quantity = std::min(bidded_supply(price_supply_ID(zone_ID), zone_ID), bidded_demand(price_demand_ID(zone_ID), zone_ID));
				Market.confirmed_supply(tick, zone_ID) += trade_quantity;
				Market.confirmed_demand(tick, zone_ID) += trade_quantity;
				bidded_supply(price_supply_ID(zone_ID), zone_ID) -= trade_quantity;
				bidded_demand(price_demand_ID(zone_ID), zone_ID) -= trade_quantity;
			}
			else{
				if(bidded_supply(price_supply_ID(zone_ID), zone_ID) > 0){
					default_price_ID(zone_ID) = price_supply_ID(zone_ID);
				}
				else{
					default_price_ID(zone_ID) = price_demand_ID(zone_ID);
				}
			}
		}
		
		// Record default market clearing prices
		Market.confirmed_price(tick, zone_ID) = Market.bidded_price(default_price_ID(zone_ID));		
	}	
}

// ------------------------------------------------------------------------------------------------
// Functions involving all markets
// ------------------------------------------------------------------------------------------------
//void Submitted_bid_calculation(market_inform &IMO_Market)

// ------------------------------------------------------------------------------------------------
// Specific functions for for flow-based markets
// ------------------------------------------------------------------------------------------------
void Flow_Based_Market_LP_Set(market_inform &Market, LP_object &Problem){
	// LP Solver initialization for flow-based market optimization
	// Warm-up once and reuse for the rest of time slices
	
	// -----------------------------------------------------------------------------------
	// Set Flow-based market object
	// -----------------------------------------------------------------------------------
	// Initialize metric tensor solver for moving around the degree of freedoms
	std::vector <Trip> Metric_trip;
	Metric_trip.reserve(3 * Market.num_zone - 5);
	for(int zone_iter = 0; zone_iter < Market.num_zone - 1; ++ zone_iter){
		Metric_trip.push_back(Trip(zone_iter, zone_iter, 2.));
		if(zone_iter != 0){
			Metric_trip.push_back(Trip(zone_iter, zone_iter - 1, -1.));
		}
		if(zone_iter != Market.num_zone - 2){
			Metric_trip.push_back(Trip(zone_iter, zone_iter + 1, -1.));
		}
	}
	Eigen::SparseMatrix <double> Metric_matrix (Market.num_zone - 1, Market.num_zone - 1);
	Metric_matrix.setFromTriplets(Metric_trip.begin(), Metric_trip.end());
	Market.dof_metric.compute(Metric_matrix);

	// -----------------------------------------------------------------------------------
	// Set LP problem object
	// -----------------------------------------------------------------------------------
	// Set dimension of the problem
	Problem.Constraints_eq_num = Market.network.num_edges + Market.network.num_vertice + 1;
	Problem.Constraints_ie_num = 0;
	Problem.Variables_num = Market.network.num_edges + 2 * Market.network.num_vertice;

	// Set objective vector
	// Since submitted bids not yet updated, will set to 0
	// The variables are ordered as {{I}, {S}, {V}}
	Problem.Objective.orig_vector = Eigen::VectorXd::Zero(Problem.Variables_num);
	Problem.Objective.varying_vector = Eigen::VectorXd::Zero(Problem.Variables_num);
	Problem.Objective.varying_vector.segment(Market.network.num_edges, Market.network.num_vertice) = Eigen::VectorXd::Ones(Market.network.num_vertice);

	// Set boudary values for equality and inequality constraints 
	// Since submitted bids not yet updated, inequality constraints for {S} will be initialized as 0
	Problem.Boundary.eq_vector = Eigen::VectorXd::Zero(Problem.Constraints_eq_num);
	Problem.Boundary.ie_orig_matrix = Eigen::MatrixXd::Zero(Problem.Variables_num + Problem.Constraints_ie_num, 2);
	Problem.Boundary.ie_orig_matrix.topRows(Market.network.num_edges) = Market.network.power_constraint;
	Problem.Boundary.ie_orig_matrix.bottomRows(Market.network.num_vertice) = Market.network.voltage_constraint;
	
	// Set sparse matrix for equality constraints
	Eigen::VectorXd Y_n_diag = Eigen::VectorXd::Zero(Market.network.num_vertice);
	std::vector<Trip> Constraint_eq_trip;
	Constraint_eq_trip.reserve(pow(Market.network.num_edges, 2) + Market.network.num_edges + 2 * Market.network.num_vertice);
	for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
		// Equality constraints of voltage at the nodes, off-diagonal terms
		Constraint_eq_trip.push_back(Trip(Market.network.num_edges + Market.network.incidence_matrix(edge_iter, 0), Market.network.num_edges + Market.network.num_vertice + Market.network.incidence_matrix(edge_iter, 1), -Market.network.admittance_vector(edge_iter)));
		Constraint_eq_trip.push_back(Trip(Market.network.num_edges + Market.network.incidence_matrix(edge_iter, 1), Market.network.num_edges + Market.network.num_vertice + Market.network.incidence_matrix(edge_iter, 0), -Market.network.admittance_vector(edge_iter)));
		Y_n_diag(Market.network.incidence_matrix(edge_iter, 0)) += Market.network.admittance_vector(edge_iter);
		Y_n_diag(Market.network.incidence_matrix(edge_iter, 1)) += Market.network.admittance_vector(edge_iter);
		
		// Equality constraints of power flows at the edges
		Constraint_eq_trip.push_back(Trip(edge_iter, Market.network.num_edges + Market.network.num_vertice + Market.network.incidence_matrix(edge_iter, 0), Market.network.admittance_vector(edge_iter)));
		Constraint_eq_trip.push_back(Trip(edge_iter, Market.network.num_edges + Market.network.num_vertice + Market.network.incidence_matrix(edge_iter, 1), -Market.network.admittance_vector(edge_iter)));
		Constraint_eq_trip.push_back(Trip(edge_iter, edge_iter, -1));	
	}
	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		// Equality constraints of voltage at the nodes, diagonal terms
		Constraint_eq_trip.push_back(Trip(Market.network.num_edges + node_iter, Market.network.num_edges + Market.network.num_vertice + node_iter, Y_n_diag(node_iter)));
		Constraint_eq_trip.push_back(Trip(Market.network.num_edges + node_iter, Market.network.num_edges + node_iter, -1));
	}
	// Equality constraint for the reference bus
	Constraint_eq_trip.push_back(Trip(Market.network.num_vertice + Market.network.num_edges, Market.network.num_vertice + Market.network.num_edges, 1));
	Problem.Constraint.eq_orig_matrix = Eigen::SparseMatrix <double> (Problem.Constraints_eq_num, Problem.Variables_num);
	Problem.Constraint.eq_orig_matrix.setFromTriplets(Constraint_eq_trip.begin(), Constraint_eq_trip.end());
	
	// Set initial feasible solution
	Problem.Solution.orig_vector = Eigen::VectorXd::Zero(Problem.Variables_num);
	
	// Initialize ldlt solver from source / sink to voltage
	Problem.Solver.ldlt.compute((Problem.Constraint.eq_orig_matrix).block(Market.network.num_edges + 1, Market.network.num_edges + Market.network.num_vertice + 1, Market.network.num_vertice - 1, Market.network.num_vertice - 1));		
}

void Flow_Based_Market_Optimization(int tick, market_inform &Market, LP_object &Problem){
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
	bool break_flag = 0;
	double tol = pow(10., -12.);
	double eps = pow(10., -10.);
	double mu = 1. - eps;
	double dS = 10.;
	double deficit_penalty = 100.;
	double error;
	double obj = 0.;
	double obj_temp;
	Eigen::VectorXi flow_dir(Market.num_zone - 1);
	Eigen::VectorXi price_ID_temp = price_ID;
	// Initial quantity should be randomize to avoid degeneracy
	Eigen::VectorXd quan = Eigen::VectorXd::Zero(Market.num_zone);
	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
		quan(zone_iter) = (std::rand() % 100 - 50) * eps;
	}
	quan -= Eigen::VectorXd::Constant(Market.num_zone, quan.sum() / Market.num_zone);
	Eigen::VectorXd quan_temp = quan;
	Eigen::VectorXd utility = Eigen::VectorXd::Zero(Market.num_zone);
	Eigen::VectorXd utility_temp = Eigen::VectorXd::Zero(Market.num_zone);
	Eigen::VectorXd increment_dot(Market.num_zone - 1);
	Eigen::VectorXd increment_coeff(Market.num_zone - 1);
	Eigen::VectorXd grad(Market.num_zone);
	Eigen::VectorXd voltage_temp = Eigen::VectorXd::Zero(Market.network.num_vertice);
	Eigen::VectorXd flow_temp(Market.network.num_edges);
	Eigen::MatrixXd Boundary_gap(Problem.Variables_num + Problem.Constraints_ie_num, 2);

	int loop_count = 0;
	//while(loop_count < 10000){
	while(!break_flag){
		//std::cout << "---------------------------------------------------------------------------\n";
		//std::cout << loop_count << "\n";
		//std::cout << "---------------------------------------------------------------------------\n";
		loop_count += 1;
		break_flag = 1;
		//std::cout << "Search all directions\n";
		for(int zone_iter = 0; zone_iter < Market.num_zone - 1; ++ zone_iter){
			// Determine flow direction
			//std::cout << "Determine flow direction\n";
			if(price_ID(zone_iter) > price_ID(zone_iter + 1)){
				flow_dir(zone_iter) = -1;
				//std::cout << "Negative flow\n";
				//std::cout << price_ID(zone_iter) << " " << price_ID(zone_iter + 1) << "\n";
			}
			else if(price_ID(zone_iter) < price_ID(zone_iter + 1)){
				flow_dir(zone_iter) = 1;
				//std::cout << "Positive flow\n";
				//std::cout << price_ID(zone_iter) << " " << price_ID(zone_iter + 1) << "\n";
			}
			else{
				increment_dot(zone_iter) = 0.;
				//std::cout << "No flow\n\n";
				//flow_dir(zone_iter) = 0;
				continue;
			}		
					
			// Update quantity after small increase / decrease
			//std::cout << "Update quantity after small increase / decrease\n";
			quan_temp(zone_iter) = quan(zone_iter) + flow_dir(zone_iter) * dS;
			quan_temp(zone_iter + 1) = quan(zone_iter + 1) - flow_dir(zone_iter) * dS;
			//std::cout << quan_temp.transpose() << "\n";
			
			// Update price after small increase / decrease
			//std::cout << "Update price after small increase / decrease\n";
			if(price_ID_temp(zone_iter) > -1 && price_ID_temp(zone_iter) < Market.price_intervals + 2){
				if(quan_temp(zone_iter) < bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)){				
					while(quan_temp(zone_iter) < bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)){
						price_ID_temp(zone_iter) -= 1;
						if(price_ID_temp(zone_iter) == -1){
							break;
						}					
					}
				}
				else if(quan_temp(zone_iter) > bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter)){
					while(quan_temp(zone_iter) > bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter)){
						if(price_ID_temp(zone_iter) == Market.price_intervals + 2){
							break;
						}
						price_ID_temp(zone_iter) += 1;					
					}
				}				
			}
			else if(quan_temp(zone_iter) < bidded_total_aggregated(0, zone_iter)){
				price_ID_temp(zone_iter) = -1;
			}
			else if(quan_temp(zone_iter) > bidded_total_aggregated(Market.price_intervals + 2, zone_iter)){
				price_ID_temp(zone_iter) = Market.price_intervals + 2;
			}
			if(price_ID_temp(zone_iter + 1) > -1 && price_ID_temp(zone_iter + 1) < Market.price_intervals + 2){
				if(quan_temp(zone_iter + 1) < bidded_total_aggregated(price_ID_temp(zone_iter + 1), zone_iter + 1)){
					while(quan_temp(zone_iter + 1) < bidded_total_aggregated(price_ID_temp(zone_iter + 1), zone_iter + 1)){
						if(price_ID_temp(zone_iter + 1) == -1){
							break;
						}
						price_ID_temp(zone_iter + 1) -= 1;					
					}
				}
				else if(quan_temp(zone_iter + 1) > bidded_total_aggregated(price_ID_temp(zone_iter + 1) + 1, zone_iter + 1)){
					while(quan_temp(zone_iter + 1) > bidded_total_aggregated(price_ID_temp(zone_iter + 1) + 1, zone_iter + 1)){
						if(price_ID_temp(zone_iter + 1) == Market.price_intervals + 2){
							break;
						}
						price_ID_temp(zone_iter + 1) += 1;					
					}
				}		
			}
			else if(quan_temp(zone_iter + 1) < bidded_total_aggregated(0, zone_iter + 1)){
				price_ID_temp(zone_iter + 1) = -1;
			}
			else if(quan_temp(zone_iter + 1) > bidded_total_aggregated(Market.price_intervals + 2, zone_iter + 1)){
				price_ID_temp(zone_iter + 1) = Market.price_intervals + 2;
			}			
			
			// Update utility function after small increase / decrease
			//std::cout << "Update utility function after small increase / decrease\n";
			if(price_ID_temp(zone_iter) == -1){
				utility_temp(zone_iter) = utility_aggregated(0, zone_iter);
				utility_temp(zone_iter) -= (quan_temp(zone_iter) - bidded_total_aggregated(0, zone_iter)) * ((deficit_penalty + 1) * Market.price_range_inflex(0) - deficit_penalty * Market.price_range_inflex(1));
			}
			else if(price_ID_temp(zone_iter) == Market.price_intervals + 2){
				utility_temp(zone_iter) = utility_aggregated(Market.price_intervals + 2, zone_iter);
				utility_temp(zone_iter) -= (quan_temp(zone_iter) - bidded_total_aggregated(Market.price_intervals + 2, zone_iter)) * ((deficit_penalty + 1) * Market.price_range_inflex(1) - deficit_penalty * Market.price_range_inflex(0));
			}
			else{
				utility_temp(zone_iter) = (bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - quan_temp(zone_iter)) * utility_aggregated(price_ID_temp(zone_iter), zone_iter)
					+ (quan_temp(zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)) * utility_aggregated(price_ID_temp(zone_iter) + 1, zone_iter);
				utility_temp(zone_iter) /= bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter);				
			}
			if(price_ID_temp(zone_iter + 1) == -1){
				utility_temp(zone_iter + 1) = utility_aggregated(0, zone_iter + 1);
				utility_temp(zone_iter + 1) -= (quan_temp(zone_iter + 1) - bidded_total_aggregated(0, zone_iter + 1)) * ((deficit_penalty + 1) * Market.price_range_inflex(0) - deficit_penalty * Market.price_range_inflex(1));				
			}
			else if(price_ID_temp(zone_iter + 1) == Market.price_intervals + 2){
				utility_temp(zone_iter + 1) = utility_aggregated(Market.price_intervals + 2, zone_iter + 1);
				utility_temp(zone_iter + 1) -= (quan_temp(zone_iter + 1) - bidded_total_aggregated(Market.price_intervals + 2, zone_iter + 1)) * ((deficit_penalty + 1) * Market.price_range_inflex(1) - deficit_penalty * Market.price_range_inflex(0));				
			}
			else{
				utility_temp(zone_iter + 1) = (bidded_total_aggregated(price_ID_temp(zone_iter + 1) + 1, zone_iter + 1) - quan_temp(zone_iter + 1)) * utility_aggregated(price_ID_temp(zone_iter + 1), zone_iter + 1)
					+ (quan_temp(zone_iter + 1) - bidded_total_aggregated(price_ID_temp(zone_iter + 1), zone_iter + 1)) * utility_aggregated(price_ID_temp(zone_iter + 1) + 1, zone_iter + 1);
				utility_temp(zone_iter + 1) /= bidded_total_aggregated(price_ID_temp(zone_iter + 1) + 1, zone_iter + 1) - bidded_total_aggregated(price_ID_temp(zone_iter + 1), zone_iter + 1);
			}
			
			// Update source / sink, voltage, and power flow
			//std::cout << "Update source / sink, voltage, and power flow\n";
			voltage_temp.tail(Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(quan_temp.tail(Market.network.num_vertice - 1));
			flow_temp = (Problem.Constraint.eq_orig_matrix).topRightCorner(Market.network.num_edges, Market.network.num_vertice) * voltage_temp;		

			// Update errors
			//std::cout << "Update errors\n";
			error = 0.;
			// Power flow errors
			for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
				error += pow(std::min(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 0), 0.) + std::max(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 1), 0.), 2.);
			}
			// Voltage errors
			for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
				error += pow(std::min(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 0), 0.) + std::max(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 1), 0.), 2.);
			}
			
			// Update objective increment for this degree of freedom
			//std::cout << "Update objective increment for this degree of freedom\n\n";
			obj_temp = (1. - mu) * utility_temp.sum() - mu * error;
			increment_dot(zone_iter) = flow_dir(zone_iter) * (obj_temp - obj);
			//std::cout << (1. - mu) * utility_temp.sum() << "\n";
			//std::cout << mu * error << "\n\n";
			
			// Return to original values
			quan_temp(zone_iter) = quan(zone_iter);
			quan_temp(zone_iter + 1) = quan(zone_iter + 1);			
			price_ID_temp(zone_iter) = price_ID(zone_iter);
			price_ID_temp(zone_iter + 1) = price_ID(zone_iter + 1);
			utility_temp(zone_iter)	= utility(zone_iter);
			utility_temp(zone_iter + 1)	= utility(zone_iter + 1);																
		}
		
		// Update gradient direction
		//std::cout << "Update gradient direction\n";
		increment_coeff = Market.dof_metric.solve(increment_dot);
		grad = Eigen::VectorXd::Zero(Market.num_zone);
		for(int zone_iter = 0; zone_iter < Market.num_zone - 1; ++ zone_iter){
			grad(zone_iter) += increment_coeff(zone_iter);
			grad(zone_iter + 1) -= increment_coeff(zone_iter);
		}
		if(grad.norm() > tol){
			grad /= grad.norm();
		}
		//std::cout << increment_coeff.transpose() << "\n";
		//std::cout << grad.transpose() << "\n\n";
		
		// Update quantity on gradient direction
		//std::cout << "Update along gradient direction\n";
		quan_temp = quan + grad * dS;
		
		// Update price after small increase / decrease
		for(int zone_iter = 0; zone_iter < Market.num_zone - 1; ++ zone_iter){
			//std::cout << "Update price after small increase / decrease\n";
			if(price_ID_temp(zone_iter) > -1 && price_ID_temp(zone_iter) < Market.price_intervals + 2){
				if(quan_temp(zone_iter) < bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)){
					while(quan_temp(zone_iter) < bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)){
						price_ID_temp(zone_iter) -= 1;
						if(price_ID_temp(zone_iter) == -1){
							break;
						}					
					}
				}
				else if(quan_temp(zone_iter) > bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter)){
					while(quan_temp(zone_iter) > bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter)){
						price_ID_temp(zone_iter) += 1;
						if(price_ID_temp(zone_iter) == Market.price_intervals + 2){
							break;
						}					
					}
				}				
			}
			else if(quan_temp(zone_iter) < bidded_total_aggregated(0, zone_iter)){
				price_ID_temp(zone_iter) = -1;
			}
			else if(quan_temp(zone_iter) > bidded_total_aggregated(Market.price_intervals + 2, zone_iter)){
				price_ID_temp(zone_iter) = Market.price_intervals + 2;
			}				
			if(price_ID_temp(zone_iter + 1) > -1 && price_ID_temp(zone_iter + 1) < Market.price_intervals + 2){
				if(quan_temp(zone_iter + 1) < bidded_total_aggregated(price_ID_temp(zone_iter + 1), zone_iter + 1)){
					while(quan_temp(zone_iter + 1) < bidded_total_aggregated(price_ID_temp(zone_iter + 1), zone_iter + 1)){
						price_ID_temp(zone_iter + 1) -= 1;
						if(price_ID_temp(zone_iter + 1) == -1){
							break;
						}					
					}
				}
				else if(quan_temp(zone_iter + 1) > bidded_total_aggregated(price_ID_temp(zone_iter + 1) + 1, zone_iter + 1)){
					while(quan_temp(zone_iter + 1) > bidded_total_aggregated(price_ID_temp(zone_iter + 1) + 1, zone_iter + 1)){
						price_ID_temp(zone_iter + 1) += 1;
						if(price_ID_temp(zone_iter + 1) == Market.price_intervals + 2){
							break;
						}					
					}
				}				
			}
			else if(quan_temp(zone_iter + 1) < bidded_total_aggregated(0, zone_iter + 1)){
				price_ID_temp(zone_iter + 1) = -1;
			}
			else if(quan_temp(zone_iter + 1) > bidded_total_aggregated(Market.price_intervals + 2, zone_iter + 1)){
				price_ID_temp(zone_iter + 1) = Market.price_intervals + 2;
			}
			
			// Update utility function after small increase / decrease
			//std::cout << "Update utility function after small increase / decrease\n";
			if(price_ID_temp(zone_iter) == -1){
				utility_temp(zone_iter) = utility_aggregated(0, zone_iter);
				utility_temp(zone_iter) -= (quan_temp(zone_iter) - bidded_total_aggregated(0, zone_iter)) * ((deficit_penalty + 1) * Market.price_range_inflex(0) - deficit_penalty * Market.price_range_inflex(1));
			}
			else if(price_ID_temp(zone_iter) == Market.price_intervals + 2){
				utility_temp(zone_iter) = utility_aggregated(Market.price_intervals + 2, zone_iter);
				utility_temp(zone_iter) -= (quan_temp(zone_iter) - bidded_total_aggregated(Market.price_intervals + 2, zone_iter)) * ((deficit_penalty + 1) * Market.price_range_inflex(1) - deficit_penalty * Market.price_range_inflex(0));
			}
			else{
				utility_temp(zone_iter) = (bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - quan_temp(zone_iter)) * utility_aggregated(price_ID_temp(zone_iter), zone_iter)
					+ (quan_temp(zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)) * utility_aggregated(price_ID_temp(zone_iter) + 1, zone_iter);
				utility_temp(zone_iter) /= bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter);				
			}
			if(price_ID_temp(zone_iter + 1) == -1){
				utility_temp(zone_iter + 1) = utility_aggregated(0, zone_iter + 1);
				utility_temp(zone_iter + 1) -= (quan_temp(zone_iter + 1) - bidded_total_aggregated(0, zone_iter + 1)) * ((deficit_penalty + 1) * Market.price_range_inflex(0) - deficit_penalty * Market.price_range_inflex(1));				
			}
			else if(price_ID_temp(zone_iter + 1) == Market.price_intervals + 2){
				utility_temp(zone_iter + 1) = utility_aggregated(Market.price_intervals + 2, zone_iter + 1);
				utility_temp(zone_iter + 1) -= (quan_temp(zone_iter + 1) - bidded_total_aggregated(Market.price_intervals + 2, zone_iter + 1)) * ((deficit_penalty + 1) * Market.price_range_inflex(1) - deficit_penalty * Market.price_range_inflex(0));				
			}
			else{
				utility_temp(zone_iter + 1) = (bidded_total_aggregated(price_ID_temp(zone_iter + 1) + 1, zone_iter + 1) - quan_temp(zone_iter + 1)) * utility_aggregated(price_ID_temp(zone_iter + 1), zone_iter + 1)
					+ (quan_temp(zone_iter + 1) - bidded_total_aggregated(price_ID_temp(zone_iter + 1), zone_iter + 1)) * utility_aggregated(price_ID_temp(zone_iter + 1) + 1, zone_iter + 1);
				utility_temp(zone_iter + 1) /= bidded_total_aggregated(price_ID_temp(zone_iter + 1) + 1, zone_iter + 1) - bidded_total_aggregated(price_ID_temp(zone_iter + 1), zone_iter + 1);
			}		
		}
		
		// Update source / sink, voltage, and power flow
		//std::cout << "Update source / sink, voltage, and power flow\n";
		voltage_temp.tail(Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(quan_temp.tail(Market.network.num_vertice - 1));
		flow_temp = (Problem.Constraint.eq_orig_matrix).topRightCorner(Market.network.num_edges, Market.network.num_vertice) * voltage_temp;		

		// Update errors
		//std::cout << "Update errors\n";
		error = 0.;
		// Power flow errors
		for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
			error += pow(std::min(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 0), 0.) + std::max(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 1), 0.), 2.);
		}
		// Voltage errors
		for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
			error += pow(std::min(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 0), 0.) + std::max(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 1), 0.), 2.);
		}
		
		// Update objective
		//std::cout << "Update objective\n";
		obj_temp = (1. - mu) * utility_temp.sum() - mu * error;
		
		// Check if objective function is actually improved
		//std::cout << "Check if objective function is actually improved\n";
		if(obj_temp > obj){			
			quan = quan_temp;
			price_ID = price_ID_temp;
			utility = utility_temp;
			obj = obj_temp;
			//std::cout << "Objective function updated\n";
			//std::cout << price_ID.transpose() << "\n";
			//std::cout << quan.transpose() << "\n";
			//std::cout << obj_temp - obj << "\n\n";			
		}
		else{
			dS /= 2.;
			//std::cout << "Objective function not updated\n";
			//std::cout << quan_temp.transpose() << "\n";
			//std::cout << quan.transpose() << "\n";
			//std::cout << obj_temp << "\n\n";
			//std::cout << dS << "\n\n";
			
			// Return to original values
			quan_temp = quan;		
			price_ID_temp = price_ID;
			utility_temp = utility;
		}
		
		if(dS > .0001){
			break_flag = 0;
		}
	}
	
	// Update final source / sink, voltage, and power flow
	//std::cout << "Update final source / sink, voltage, and power flow\n";
	Problem.Solution.orig_vector.segment(Market.network.num_edges, Market.network.num_vertice) = quan;
	Problem.Solution.orig_vector.tail(Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(Problem.Solution.orig_vector.segment(Market.network.num_edges + 1, Market.network.num_vertice - 1));
	Problem.Solution.orig_vector.head(Market.network.num_edges) = (Problem.Constraint.eq_orig_matrix).topRightCorner(Market.network.num_edges, Market.network.num_vertice) * Problem.Solution.orig_vector.tail(Market.network.num_vertice);		
	//std::cout << Problem.Solution.orig_vector.segment(Market.network.num_edges, Market.network.num_vertice).transpose() << "\n\n";
	//std::cout << Problem.Solution.orig_vector.transpose() << "\n\n";
	std::cout << quan.minCoeff() << " " << quan.maxCoeff() << "\n";
	std::cout << (Problem.Solution.orig_vector.tail(Market.network.num_vertice - 1)).minCoeff() << " " << (Problem.Solution.orig_vector.tail(Market.network.num_vertice - 1)).maxCoeff() << "\n";
	std::cout << (Problem.Solution.orig_vector.head(Market.network.num_edges)).minCoeff() << " " << (Problem.Solution.orig_vector.head(Market.network.num_edges)).maxCoeff() << "\n\n";		
}