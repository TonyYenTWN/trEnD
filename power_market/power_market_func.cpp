// Source file for the functions of the power market
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
			while(bidded_demand(price_demand_ID(zone_ID), zone_ID) == 0.){
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
			while(bidded_supply(price_supply_ID(zone_ID), zone_ID) == 0.){
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
void Submitted_bid_calculation(int tick, DSO_Markets &DSO_Markets, market_inform &TSO_Market, market_inform &International_Market, network_inform &Power_network_inform, std::string fin_point_demand){
	// Calculation of submit bids at the beginning of each time slice
	auto fin_point_demand_dim = get_file_dim(fin_point_demand);
	auto point_demand_inform = read_file(fin_point_demand_dim[0], fin_point_demand_dim[1], fin_point_demand);
	
	// Initialize submit bids of markets
	Market_Initialization(TSO_Market);
	Market_Initialization(International_Market);
	for(int DSO_iter = 0; DSO_iter < DSO_Markets.markets.size(); ++ DSO_iter){
		Market_Initialization(DSO_Markets.markets[DSO_iter]);
	}

	// Declare variables for the loops
	int bz_ID;
	int DSO_ID;	
	int node_ID;
	int x_ID;
	int y_ID;	
	int point_ID;
	double bid_quan;
	Eigen::VectorXd bid_vec;
	
	// Trivial case: demand at each point are 100% inflexible
	for(int point_iter = 0; point_iter < Power_network_inform.points.bidding_zone.size(); ++ point_iter){
		node_ID = Power_network_inform.points.node(point_iter);
		DSO_ID = Power_network_inform.nodes.cluster(node_ID);
		bz_ID = Power_network_inform.points.bidding_zone(point_iter);
		bid_quan = point_demand_inform(point_iter, 0) * Power_network_inform.points.population_density(point_iter); //* Power_network_inform.points.point_area;
		// nominal demand currently wrong in processed files, should change them and then multiply area of a point later
		
		DSO_Markets.markets[DSO_ID].submitted_demand(DSO_Markets.markets[DSO_ID].price_intervals + 1, Power_network_inform.points.in_cluster_ID(point_iter)) = bid_quan;
		TSO_Market.submitted_demand(DSO_Markets.markets[DSO_ID].price_intervals + 1, node_ID) += bid_quan;
		International_Market.submitted_demand(DSO_Markets.markets[DSO_ID].price_intervals + 1, bz_ID) += bid_quan;
	}
	
	// Supply at each point (LV power plants) / node (HV power plants)
	for(int hydro_iter = 0; hydro_iter < Power_network_inform.plants.hydro.node.size(); ++ hydro_iter){
		if(Power_network_inform.plants.hydro.cap(hydro_iter) >= 20.){
			node_ID = Power_network_inform.plants.hydro.node(hydro_iter);
			DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);
			bid_vec = International_Market.merit_order_curve.col(Power_network_inform.nodes.bidding_zone(node_ID)) 
				* Power_network_inform.plants.hydro.cap(hydro_iter) / (International_Market.merit_order_curve.col(Power_network_inform.nodes.bidding_zone(node_ID)).sum());	
			
			DSO_Markets.markets[DSO_ID].submitted_supply.col(Power_network_inform.DSO_cluster[DSO_ID].points_ID.size() + Power_network_inform.nodes.in_cluster_ID(node_ID)) += bid_vec;
		}
		else{
			x_ID = int((Power_network_inform.plants.hydro.x(hydro_iter) - Power_network_inform.points.x.minCoeff()) / Power_network_inform.points.grid_length + .5);
			y_ID = int((Power_network_inform.plants.hydro.y(hydro_iter) - Power_network_inform.points.y.minCoeff()) / Power_network_inform.points.grid_length + .5);
			point_ID = Power_network_inform.points.coordinate_grid(x_ID, y_ID);
			if(point_ID == -1){
				continue;
			}
			node_ID = Power_network_inform.points.node(point_ID);
			DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);
			bid_vec = International_Market.merit_order_curve.col(Power_network_inform.points.bidding_zone(point_ID))
				* Power_network_inform.plants.hydro.cap(hydro_iter) / (International_Market.merit_order_curve.col(Power_network_inform.points.bidding_zone(point_ID)).sum());
			
			DSO_Markets.markets[DSO_ID].submitted_supply.col(Power_network_inform.points.in_cluster_ID(point_ID)) += bid_vec;
		}
		TSO_Market.submitted_supply.col(node_ID) += bid_vec;
		International_Market.submitted_supply.col(bz_ID) += bid_vec;		
	}
}

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
	//Problem.Objective.varying_vector = Eigen::VectorXd::Zero(Problem.Variables_num);
	//Problem.Objective.varying_vector.segment(Market.network.num_edges, Market.network.num_vertice) = Eigen::VectorXd::Ones(Market.network.num_vertice);

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
	//Problem.Solver.qr.compute(Problem.Constraint.eq_orig_matrix * Problem.Constraint.eq_orig_matrix.transpose());
	
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
	Eigen::VectorXi Zero_node = Eigen::VectorXi::Zero(Market.num_zone);
	std::vector <int> Non_zero;
	Non_zero.reserve(Market.num_zone);
	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
		if(abs(bidded_total_aggregated(0, zone_iter) - bidded_total_aggregated(Market.price_intervals + 2, zone_iter)) <= 0.){
			Zero_node(zone_iter) = 1;
		}
		else{
			Non_zero.push_back(zone_iter);
		}
	}
	
	// Declare variables for the main loop
	bool divide_flag;
	int node_ref_ID;
	int zone_iter;
	double tol = pow(10., -12.);
	double eps = pow(10., -10.);
	double mu = 1 - eps;
	double dS_max = .2 * (bidded_total_aggregated.bottomRows(1)).lpNorm <1> () / bidded_total_aggregated.cols();
	//double dS_max = 10.;
	double dS = dS_max;
	double error;
	double obj = 0.;
	double obj_temp;
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
	//while(loop_count < 30){
	while(dS > .01 * dS_max){
		loop_count += 1;
		divide_flag = 1;
		//for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
		for(int rand_draw = 0; rand_draw < 5000; ++ rand_draw){
			zone_iter = std::rand() % Non_zero.size();
			node_ref_ID = std::rand() % Non_zero.size();
			if(zone_iter == node_ref_ID){
				continue;
			}			
			
			for(int dir_iter = 0; dir_iter < 2; ++ dir_iter){
				//Check if direction is plausible
				if(price_ID_temp(Non_zero[zone_iter]) < price_ID_temp(Non_zero[node_ref_ID]) && dir_iter == 0){
					continue;
				}
				else if(price_ID_temp(Non_zero[zone_iter]) > price_ID_temp(Non_zero[node_ref_ID]) && dir_iter == 1){
					continue;
				}
				
				// Update quantity after small increase / decrease
				quan_temp(Non_zero[zone_iter]) += 2 * (dir_iter - .5) * dS;
				quan_temp(Non_zero[node_ref_ID]) -= 2 * (dir_iter - .5) * dS;
				
				// Update price after small increase / decrease
				if(quan_temp(Non_zero[zone_iter]) < bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])){
					while(quan_temp(Non_zero[zone_iter]) < bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])){
						if(price_ID_temp(Non_zero[zone_iter]) == 0){
							break;
						}
						price_ID_temp(Non_zero[zone_iter]) -= 1;
					}
				}
				else if(quan_temp(Non_zero[zone_iter]) > bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter])){
					while(quan_temp(Non_zero[zone_iter]) > bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter])){
						if(price_ID_temp(Non_zero[zone_iter]) == Market.price_intervals + 1){
							break;
						}
						price_ID_temp(Non_zero[zone_iter]) += 1;
					}
				}
				if(quan_temp(Non_zero[node_ref_ID]) < bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID])){
					while(quan_temp(Non_zero[node_ref_ID]) < bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID])){
						if(price_ID_temp(Non_zero[node_ref_ID]) == 0){
							break;
						}
						price_ID_temp(Non_zero[node_ref_ID]) -= 1;
					}
				}
				else if(quan_temp(Non_zero[node_ref_ID]) > bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID])){
					while(quan_temp(Non_zero[node_ref_ID]) > bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID])){
						if(price_ID_temp(Non_zero[node_ref_ID]) == Market.price_intervals + 1){
							break;
						}
						price_ID_temp(Non_zero[node_ref_ID]) += 1;
					}
				}				
				
				// Update utility after small increase / decrease
				if(quan_temp(Non_zero[zone_iter]) < bidded_total_aggregated(0, Non_zero[zone_iter])){
					utility_temp(Non_zero[zone_iter]) = utility_aggregated(0, Non_zero[zone_iter]);
					utility_temp(Non_zero[zone_iter]) -= (quan_temp(Non_zero[zone_iter]) - bidded_total_aggregated(0, Non_zero[zone_iter])) * Market.price_range_inflex(0);
				}
				else if(quan_temp(Non_zero[zone_iter]) > bidded_total_aggregated(Market.price_intervals + 2, Non_zero[zone_iter]))	{
					utility_temp(Non_zero[zone_iter]) = utility_aggregated(Market.price_intervals + 2, Non_zero[zone_iter]);
					utility_temp(Non_zero[zone_iter]) -= (quan_temp(Non_zero[zone_iter]) - bidded_total_aggregated(Market.price_intervals + 2, Non_zero[zone_iter])) * Market.price_range_inflex(1);
				}
				else{
					if(bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]) != 0){
						utility_temp(Non_zero[zone_iter]) = (bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]) - quan_temp(Non_zero[zone_iter])) * utility_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])
							+ (quan_temp(Non_zero[zone_iter]) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])) * utility_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]);
						utility_temp(Non_zero[zone_iter]) /= bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]);	
					}
					else{
						utility_temp(Non_zero[zone_iter]) = utility_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]);
					}
				}
				if(quan_temp(Non_zero[node_ref_ID]) < bidded_total_aggregated(0, Non_zero[node_ref_ID])){
					utility_temp(Non_zero[node_ref_ID]) = utility_aggregated(0, Non_zero[node_ref_ID]);
					utility_temp(Non_zero[node_ref_ID]) -= (quan_temp(Non_zero[node_ref_ID]) - bidded_total_aggregated(0, Non_zero[node_ref_ID])) * Market.price_range_inflex(0);
				}
				else if(quan_temp(Non_zero[node_ref_ID]) > bidded_total_aggregated(Market.price_intervals + 2, Non_zero[node_ref_ID]))	{
					utility_temp(Non_zero[node_ref_ID]) = utility_aggregated(Market.price_intervals + 2, Non_zero[node_ref_ID]);
					utility_temp(Non_zero[node_ref_ID]) -= (quan_temp(Non_zero[node_ref_ID]) - bidded_total_aggregated(Market.price_intervals + 2, Non_zero[node_ref_ID])) * Market.price_range_inflex(1);
				}
				else{
					if(bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID]) - bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID]) != 0){
						utility_temp(Non_zero[node_ref_ID]) = (bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID]) - quan_temp(Non_zero[node_ref_ID])) * utility_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID])
							+ (quan_temp(Non_zero[node_ref_ID]) - bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID])) * utility_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID]);
						utility_temp(Non_zero[node_ref_ID]) /= bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID]) - bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID]);	
					}
					else{
						utility_temp(Non_zero[node_ref_ID]) = utility_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID]);
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
				//std::cout << Non_zero[zone_iter] << " " << Non_zero[node_ref_ID] << " " << obj_temp << "\n";
				
				// Check if objective should be updated
				if(obj_temp > obj){
					// Update if objective has been improved
					quan(Non_zero[zone_iter]) = quan_temp(Non_zero[zone_iter]);
					quan(Non_zero[node_ref_ID]) = quan_temp(Non_zero[node_ref_ID]);
					voltage = voltage_temp;
					flow = flow_temp;
					price_ID(Non_zero[zone_iter]) = price_ID_temp(Non_zero[zone_iter]);
					price_ID(Non_zero[node_ref_ID]) = price_ID_temp(Non_zero[node_ref_ID]);
					utility(Non_zero[zone_iter]) = utility_temp(Non_zero[zone_iter]);
					utility(Non_zero[node_ref_ID]) = utility_temp(Non_zero[node_ref_ID]);
					obj = obj_temp;	
					divide_flag = 0;
					break;
				}
				else{
					// Return to origin
					quan_temp(Non_zero[zone_iter]) = quan(Non_zero[zone_iter]);
					quan_temp(Non_zero[node_ref_ID]) = quan(Non_zero[node_ref_ID]);
					price_ID_temp(Non_zero[zone_iter]) = price_ID(Non_zero[zone_iter]);
					price_ID_temp(Non_zero[node_ref_ID]) = price_ID(Non_zero[node_ref_ID]);
					utility_temp(Non_zero[zone_iter]) = utility(Non_zero[zone_iter]);
					utility_temp(Non_zero[node_ref_ID]) = utility(Non_zero[node_ref_ID]);					
				}																					
			}	
		}
		
		// Check if increment should be smaller
		if(divide_flag){
			dS /= 2.;
			//std::cout << dS << "\n";
		}
	}

	Problem.Solution.orig_vector.segment(Market.network.num_edges, Market.network.num_vertice) = quan;
	Problem.Solution.orig_vector.tail(Market.network.num_vertice) = voltage;
	Problem.Solution.orig_vector.head(Market.network.num_edges) = flow;	
	std::cout << quan.minCoeff() << " " << quan.maxCoeff() << " " << .5 * quan.array().abs().sum() << "\n";
	std::cout << voltage.minCoeff() << " " << voltage.maxCoeff() << "\n";
	std::cout << flow.minCoeff() << " " << flow.maxCoeff() << "\n";	
	//std::cout << "\n" << grad.transpose() << "\n\n";
	//std::cout << "\n" << obj << "\n\n";
	//std::cout << "\n" << quan.transpose() << "\n\n";
	//std::cout << "\n" << voltage.transpose() << "\n\n";
	//std::cout << "\n" << flow.transpose() << "\n\n";
}