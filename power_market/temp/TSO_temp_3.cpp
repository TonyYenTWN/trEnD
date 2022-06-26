// Source File for re-dispatch and tertiary control reserve market clearing of TSO in Norway
#include <iostream>
//#include <chrono>
//#include "../../basic/LP_gpa.cpp"
#include "../../basic/LP_gpa_fast.cpp"
#include "../power_market.cpp"

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

	// Initialize metric tensor solver for degree of freedoms
	Market_Solver_Set(TSO_Market);
	
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
//	TSO_Market.submitted_supply(0, 0) = 15.;
//	TSO_Market.submitted_supply(0, 1) = 35.;
//	TSO_Market.submitted_supply(0, 2) = 45.;
//	TSO_Market.submitted_demand(TSO_Market.price_intervals + 1, 0) = 35.;
//	TSO_Market.submitted_demand(TSO_Market.price_intervals + 1, 1) = 22.5;
//	TSO_Market.submitted_demand(TSO_Market.price_intervals + 1, 2) = 10.;
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
	TSO_Market.network.voltage_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, -10.);
	TSO_Market.network.voltage_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, 10.);
	TSO_Market.network.power_constraint = Eigen::MatrixXd(TSO_Market.network.num_edges, 2);
	TSO_Market.network.power_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, -150.);
	TSO_Market.network.power_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, 150.);

	// Initialization of output variables
	TSO_Market.confirmed_supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_price = Eigen::MatrixXd(Time, TSO_Market.num_zone);
	TSO_Market.network.confirmed_power = Eigen::MatrixXd(Time, TSO_Market.network.num_edges);
	
	// Initialize metric tensor solver for degree of freedoms
	Market_Solver_Set(TSO_Market);
	
	// For the trivial case only: initialize submitted supply and demand bids at each node
	Market_Initialization(TSO_Market);
	TSO_Market.submitted_supply.leftCols(TSO_Market.num_zone / 2) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 1.);
	TSO_Market.submitted_demand.rightCols(TSO_Market.num_zone / 2) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 1.);
}

//double impedence_conversion(Eigen::MatrixXd pu_dc_inform, double voltage){
//	int v_iter = 0;
//	while(pu_dc_inform(v_iter, 0) != voltage){
//		v_iter += 1;
//	}
//	
//	return(pu_dc_inform(v_iter, 1));
//}

void TSO_Market_Set(market_inform &TSO_Market, int Time, std::string fin_node, std::string fin_edge, std::string fin_pu_dc){
	// Read power network data
	auto fin_node_dim = get_file_dim(fin_node);
	auto fin_edge_dim = get_file_dim(fin_edge);
	auto fin_pu_dc_dim = get_file_dim(fin_pu_dc);
	auto node_inform = read_file(fin_node_dim[0], fin_node_dim[1], fin_node);
	auto edge_inform = read_file(fin_edge_dim[0], fin_edge_dim[1], fin_edge);
	auto pu_dc_inform = read_file(fin_pu_dc_dim[0], fin_pu_dc_dim[1], fin_pu_dc);
	//double x_L = 5. * pow(10., -4.);		// Inductance per meter of transmission line
	
	// Input parameters of TSO market
	TSO_Market.num_zone = fin_node_dim[0];
	TSO_Market.time_intervals = Time;
	TSO_Market.price_intervals = 600;
	TSO_Market.price_range_inflex << -500., 3000.;
	TSO_Market.price_range_flex << -100., 500.;
	TSO_Market.bidded_price = Eigen::VectorXd(TSO_Market.price_intervals + 2);
	TSO_Market.bidded_price(0) = TSO_Market.price_range_inflex(0);
	TSO_Market.bidded_price.array().tail(1) = TSO_Market.price_range_inflex(1);
	TSO_Market.bidded_price.array().segment(1, TSO_Market.price_intervals) = Eigen::VectorXd::LinSpaced(TSO_Market.price_intervals, TSO_Market.price_range_flex(0) + .5 * (TSO_Market.price_range_flex(1) - TSO_Market.price_range_flex(0)) / TSO_Market.price_intervals, TSO_Market.price_range_flex(1) - .5 * (TSO_Market.price_range_flex(1) - TSO_Market.price_range_flex(0)) / TSO_Market.price_intervals);
	
	// Set compact incidence matrix and edge admittance matrix
	TSO_Market.network.num_vertice = TSO_Market.num_zone;
	TSO_Market.network.num_edges = fin_edge_dim[0];
	TSO_Market.network.incidence_matrix = Eigen::MatrixXi(TSO_Market.network.num_edges, 2);
	TSO_Market.network.admittance_vector = Eigen::VectorXd(TSO_Market.network.num_edges);
	for(int edge_iter = 0; edge_iter < TSO_Market.network.num_edges; ++ edge_iter){
		TSO_Market.network.incidence_matrix(edge_iter, 0) = edge_inform(edge_iter, 0) - 1;
		TSO_Market.network.incidence_matrix(edge_iter, 1) = edge_inform(edge_iter, 1) - 1;
		//TSO_Market.network.admittance_vector(edge_iter) = impedence_conversion(pu_dc_inform, edge_inform(edge_iter, 4)) / x_L / edge_inform(edge_iter, 5);
		TSO_Market.network.admittance_vector(edge_iter) = edge_inform(edge_iter, 2);
	}

	// Set voltage and power constraints at each edge
	TSO_Market.network.voltage_constraint = Eigen::MatrixXd(TSO_Market.network.num_vertice, 2);
	TSO_Market.network.voltage_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, -pi / 18);
	TSO_Market.network.voltage_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, pi / 18);
	TSO_Market.network.power_constraint = Eigen::MatrixXd(TSO_Market.network.num_edges, 2);
	TSO_Market.network.power_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, -5.);
	TSO_Market.network.power_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, 5.);

	// Initialization of output variables
	TSO_Market.confirmed_supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_price = Eigen::MatrixXd(Time, TSO_Market.num_zone);
	TSO_Market.network.confirmed_power = Eigen::MatrixXd(Time, TSO_Market.network.num_edges);

	// Initialize metric tensor solver for degree of freedoms
	Market_Solver_Set(TSO_Market);

	// Trivial initialization at the nodes
	Market_Initialization(TSO_Market);
	TSO_Market.submitted_supply.leftCols(TSO_Market.num_zone / 2) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 1.);
	TSO_Market.submitted_demand.leftCols(TSO_Market.num_zone / 2) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 2.);
	TSO_Market.submitted_demand.rightCols(TSO_Market.num_zone / 2) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 1.);
	TSO_Market.submitted_supply.rightCols(TSO_Market.num_zone / 2) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 2.);
}

void TSO_LP_Set(market_inform &TSO_Market, LP_object &Problem){
	// Set dimension of the problem
	Problem.Constraints_eq_num = TSO_Market.network.num_edges + TSO_Market.network.num_vertice + 1;
	Problem.Constraints_ie_num = 0;
	Problem.Variables_num = TSO_Market.network.num_edges + 2 * TSO_Market.network.num_vertice;

	// Set objective vector
	// Since submitted bids not yet updated, will set to 0
	// The variables are ordered as {{V}, {S}, {I}} -> changed
	// The variables are ordered as {{I}, {S}, {V}}
	Problem.Objective.orig_vector = Eigen::VectorXd::Zero(Problem.Variables_num);
	Problem.Objective.varying_vector = Eigen::VectorXd::Zero(Problem.Variables_num);
	Problem.Objective.varying_vector.segment(TSO_Market.network.num_edges, TSO_Market.network.num_vertice) = Eigen::VectorXd::Ones(TSO_Market.network.num_vertice);

	// Set boudary values for equality and inequality constraints 
	// Since submitted bids not yet updated, inequality constraints for {S} will set to 0
	Problem.Boundary.eq_vector = Eigen::VectorXd::Zero(Problem.Constraints_eq_num);
	Problem.Boundary.ie_orig_matrix = Eigen::MatrixXd::Zero(Problem.Variables_num + Problem.Constraints_ie_num, 2);
	Problem.Boundary.ie_orig_matrix.topRows(TSO_Market.network.num_edges) = TSO_Market.network.power_constraint;
	Problem.Boundary.ie_orig_matrix.bottomRows(TSO_Market.network.num_vertice) = TSO_Market.network.voltage_constraint;
	
	// Set sparse matrix for equality constraints
	Eigen::VectorXd Y_n_diag = Eigen::VectorXd::Zero(TSO_Market.network.num_vertice);
	std::vector<Trip> Constraint_eq_trip;
	Constraint_eq_trip.reserve(pow(TSO_Market.network.num_edges, 2) + TSO_Market.network.num_edges + 2 * TSO_Market.network.num_vertice);
	for(int edge_iter = 0; edge_iter < TSO_Market.network.num_edges; ++ edge_iter){
		// Equality constraints of voltage at the nodes, off-diagonal terms
		Constraint_eq_trip.push_back(Trip(TSO_Market.network.num_edges + TSO_Market.network.incidence_matrix(edge_iter, 0), TSO_Market.network.num_edges + TSO_Market.network.num_vertice + TSO_Market.network.incidence_matrix(edge_iter, 1), -TSO_Market.network.admittance_vector(edge_iter)));
		Constraint_eq_trip.push_back(Trip(TSO_Market.network.num_edges + TSO_Market.network.incidence_matrix(edge_iter, 1), TSO_Market.network.num_edges + TSO_Market.network.num_vertice + TSO_Market.network.incidence_matrix(edge_iter, 0), -TSO_Market.network.admittance_vector(edge_iter)));
		Y_n_diag(TSO_Market.network.incidence_matrix(edge_iter, 0)) += TSO_Market.network.admittance_vector(edge_iter);
		Y_n_diag(TSO_Market.network.incidence_matrix(edge_iter, 1)) += TSO_Market.network.admittance_vector(edge_iter);
		
		// Equality constraints of power flows at the edges
		Constraint_eq_trip.push_back(Trip(edge_iter, TSO_Market.network.num_edges + TSO_Market.network.num_vertice + TSO_Market.network.incidence_matrix(edge_iter, 0), TSO_Market.network.admittance_vector(edge_iter)));
		Constraint_eq_trip.push_back(Trip(edge_iter, TSO_Market.network.num_edges + TSO_Market.network.num_vertice + TSO_Market.network.incidence_matrix(edge_iter, 1), -TSO_Market.network.admittance_vector(edge_iter)));
		Constraint_eq_trip.push_back(Trip(edge_iter, edge_iter, -1));	
	}
	for(int node_iter = 0; node_iter < TSO_Market.network.num_vertice; ++ node_iter){
		// Equality constraints of voltage at the nodes, diagonal terms
		Constraint_eq_trip.push_back(Trip(TSO_Market.network.num_edges + node_iter, TSO_Market.network.num_edges + TSO_Market.network.num_vertice + node_iter, Y_n_diag(node_iter)));
		Constraint_eq_trip.push_back(Trip(TSO_Market.network.num_edges + node_iter, TSO_Market.network.num_edges + node_iter, -1));
	}
	// Equality constraint for the reference bus
	Constraint_eq_trip.push_back(Trip(TSO_Market.network.num_vertice + TSO_Market.network.num_edges, TSO_Market.network.num_vertice + TSO_Market.network.num_edges, 1));
	Problem.Constraint.eq_orig_matrix = Eigen::SparseMatrix <double> (Problem.Constraints_eq_num, Problem.Variables_num);
	Problem.Constraint.eq_orig_matrix.setFromTriplets(Constraint_eq_trip.begin(), Constraint_eq_trip.end());
	
	// Set initial feasible solution
	Problem.Solution.orig_vector = Eigen::VectorXd::Zero(Problem.Variables_num);
	
	// Initialize ldlt solver from source / sink to voltage
	Problem.Solver.ldlt.compute((Problem.Constraint.eq_orig_matrix).block(TSO_Market.network.num_edges + 1, TSO_Market.network.num_edges + TSO_Market.network.num_vertice + 1, TSO_Market.network.num_vertice - 1, TSO_Market.network.num_vertice - 1));
}

void TSO_Market_Optimization(int tick, market_inform &TSO_Market, LP_object &Problem){
	Eigen::MatrixXd bidded_supply = TSO_Market.submitted_supply;
	Eigen::MatrixXd bidded_demand = TSO_Market.submitted_demand;
	
	// Initial market clearing within each nodes
	Eigen::VectorXi price_ID(TSO_Market.num_zone);
	Market_clearing_nodal(tick, TSO_Market, price_ID, bidded_supply, bidded_demand);
	
	// Initialization of process variables for the main optimization loop
	Eigen::MatrixXd bidded_supply_export = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	Eigen::MatrixXd bidded_demand_export = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	Eigen::MatrixXd bidded_supply_import = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	Eigen::MatrixXd bidded_demand_import = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	Eigen::MatrixXd bidded_supply_aggregated = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 3, TSO_Market.num_zone);
	Eigen::MatrixXd bidded_demand_aggregated = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 3, TSO_Market.num_zone);
	Eigen::MatrixXd bidded_total_aggregated = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 3, TSO_Market.num_zone);
	Eigen::MatrixXd utility_aggregated = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 3, TSO_Market.num_zone);
	#pragma omp parallel
	{
		#pragma omp for
		for(int zone_iter = 0; zone_iter < TSO_Market.num_zone; ++ zone_iter){
			// Supply curves
			bidded_supply_export.col(zone_iter).array().tail(TSO_Market.price_intervals + 1 - price_ID(zone_iter)) = TSO_Market.submitted_supply.col(zone_iter).array().tail(TSO_Market.price_intervals + 1 - price_ID(zone_iter));
			bidded_supply_export(price_ID(zone_iter), zone_iter) = bidded_supply(price_ID(zone_iter), zone_iter);
			bidded_supply_import.col(zone_iter).array().head(price_ID(zone_iter)) = TSO_Market.submitted_supply.col(zone_iter).array().head(price_ID(zone_iter));
			bidded_supply_import(price_ID(zone_iter), zone_iter) = TSO_Market.submitted_supply(price_ID(zone_iter), zone_iter) - bidded_supply(price_ID(zone_iter), zone_iter);
			bidded_supply_aggregated(0, zone_iter) = -bidded_supply_import.col(zone_iter).sum();
			
			// Demand curves
			bidded_demand_export.col(zone_iter).array().tail(TSO_Market.price_intervals + 1 - price_ID(zone_iter)) = TSO_Market.submitted_demand.col(zone_iter).array().tail(TSO_Market.price_intervals + 1 - price_ID(zone_iter));
			bidded_demand_export(price_ID(zone_iter), zone_iter) = TSO_Market.submitted_demand(price_ID(zone_iter), zone_iter) - bidded_demand(price_ID(zone_iter), zone_iter);
			bidded_demand_import.col(zone_iter).array().head(price_ID(zone_iter)) = TSO_Market.submitted_demand.col(zone_iter).array().head(price_ID(zone_iter));
			bidded_demand_import(price_ID(zone_iter), zone_iter) = bidded_demand(price_ID(zone_iter), zone_iter);
			bidded_demand_aggregated(0, zone_iter) = -bidded_demand_import.col(zone_iter).sum();		
		}
	}
	// Aggregated import-export curves for supply and demand
	for(int price_iter = 1; price_iter < TSO_Market.price_intervals + 3; ++ price_iter){
		bidded_supply_aggregated.row(price_iter) = bidded_supply_aggregated.row(price_iter - 1) + bidded_supply_import.row(price_iter - 1) + bidded_supply_export.row(price_iter - 1);
		bidded_demand_aggregated.row(price_iter) = bidded_demand_aggregated.row(price_iter - 1) + bidded_demand_import.row(price_iter - 1) + bidded_demand_export.row(price_iter - 1);
	}
	bidded_total_aggregated = bidded_supply_aggregated + bidded_demand_aggregated;
	// Aggregated utility function for each bidding zone
	for(int zone_iter = 0; zone_iter < TSO_Market.num_zone; ++ zone_iter){	
		// Demand
		utility_aggregated(price_ID(zone_iter), zone_iter) = (bidded_demand_import(price_ID(zone_iter), zone_iter) + bidded_supply_import(price_ID(zone_iter), zone_iter)) * TSO_Market.bidded_price(price_ID(zone_iter));
		if(price_ID(zone_iter) > 0){
			for(int price_iter = price_ID(zone_iter) - 1; price_iter >= 0; -- price_iter){
				utility_aggregated(price_iter, zone_iter) = utility_aggregated(price_iter + 1, zone_iter) + (bidded_demand_import(price_iter, zone_iter) + bidded_supply_import(price_iter, zone_iter)) * TSO_Market.bidded_price(price_iter);
			}
		}
		
		// Supply
		utility_aggregated(price_ID(zone_iter) + 1, zone_iter) = -(bidded_demand_export(price_ID(zone_iter), zone_iter) + bidded_supply_export(price_ID(zone_iter), zone_iter)) * TSO_Market.bidded_price(price_ID(zone_iter));
		if(price_ID(zone_iter) < TSO_Market.price_intervals){
			for(int price_iter = price_ID(zone_iter) + 2; price_iter < TSO_Market.price_intervals + 3; ++ price_iter){
				utility_aggregated(price_iter, zone_iter) = utility_aggregated(price_iter - 1, zone_iter) - (bidded_demand_export(price_iter - 1, zone_iter) + bidded_supply_export(price_iter - 1, zone_iter)) * TSO_Market.bidded_price(price_iter - 1);
			}			
		}
	}

	// Declare variables for the main loop
	bool break_flag = 0;
	double tol = pow(10., -12.);
	double eps = pow(10., -10.);
	double mu = 1. - eps;
	double dS = 10.;
	double error;
	double obj = 0.;
	double obj_temp;
	Eigen::VectorXi flow_dir(TSO_Market.num_zone - 1);
	Eigen::VectorXi price_ID_temp = price_ID;
	// Initial quantity should be randomize to avoid degeneracy
	Eigen::VectorXd quan = Eigen::VectorXd::Zero(TSO_Market.num_zone);
	for(int zone_iter = 0; zone_iter < TSO_Market.num_zone; ++ zone_iter){
		quan(zone_iter) = (std::rand() % 100 - 50) * eps;
	}
	quan -= Eigen::VectorXd::Constant(TSO_Market.num_zone, quan.sum() / TSO_Market.num_zone);
	Eigen::VectorXd quan_temp = quan;
	Eigen::VectorXd utility = Eigen::VectorXd::Zero(TSO_Market.num_zone);
	Eigen::VectorXd utility_temp = Eigen::VectorXd::Zero(TSO_Market.num_zone);
	Eigen::VectorXd increment_dot(TSO_Market.num_zone - 1);
	Eigen::VectorXd increment_coeff(TSO_Market.num_zone - 1);
	Eigen::VectorXd grad(TSO_Market.num_zone);
	Eigen::VectorXd voltage_temp = Eigen::VectorXd::Zero(TSO_Market.network.num_vertice);
	Eigen::VectorXd flow_temp(TSO_Market.network.num_edges);
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
		for(int zone_iter = 0; zone_iter < TSO_Market.num_zone - 1; ++ zone_iter){
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
				increment_dot(zone_iter) = 0;
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
			if(price_ID_temp(zone_iter) > -1 && price_ID_temp(zone_iter) < TSO_Market.price_intervals + 2){
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
						if(price_ID_temp(zone_iter) == TSO_Market.price_intervals + 2){
							break;
						}
						price_ID_temp(zone_iter) += 1;					
					}
				}				
			}
			else if(quan_temp(zone_iter) < bidded_total_aggregated(0, zone_iter)){
				price_ID_temp(zone_iter) = -1;
			}
			else if(quan_temp(zone_iter) > bidded_total_aggregated(TSO_Market.price_intervals + 2, zone_iter)){
				price_ID_temp(zone_iter) = TSO_Market.price_intervals + 2;
			}
			if(price_ID_temp(zone_iter + 1) > -1 && price_ID_temp(zone_iter + 1) < TSO_Market.price_intervals + 2){
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
						if(price_ID_temp(zone_iter + 1) == TSO_Market.price_intervals + 2){
							break;
						}
						price_ID_temp(zone_iter + 1) += 1;					
					}
				}		
			}
			else if(quan_temp(zone_iter + 1) < bidded_total_aggregated(0, zone_iter + 1)){
				price_ID_temp(zone_iter + 1) = -1;
			}
			else if(quan_temp(zone_iter + 1) > bidded_total_aggregated(TSO_Market.price_intervals + 2, zone_iter + 1)){
				price_ID_temp(zone_iter + 1) = TSO_Market.price_intervals + 2;
			}			
			
			// Update utility function after small increase / decrease
			//std::cout << "Update utility function after small increase / decrease\n";
			if(price_ID_temp(zone_iter) == -1){
				utility_temp(zone_iter) = utility_aggregated(0, zone_iter);
				utility_temp(zone_iter) -= (quan_temp(zone_iter) - bidded_total_aggregated(0, zone_iter)) * (2 * TSO_Market.price_range_inflex(0) - TSO_Market.price_range_inflex(1));
			}
			else if(price_ID_temp(zone_iter) == TSO_Market.price_intervals + 2){
				utility_temp(zone_iter) = utility_aggregated(TSO_Market.price_intervals + 2, zone_iter);
				utility_temp(zone_iter) -= (quan_temp(zone_iter) - bidded_total_aggregated(TSO_Market.price_intervals + 2, zone_iter)) * (2 * TSO_Market.price_range_inflex(1) - TSO_Market.price_range_inflex(0));
			}
			else{
				utility_temp(zone_iter) = (bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - quan_temp(zone_iter)) * utility_aggregated(price_ID_temp(zone_iter), zone_iter)
					+ (quan_temp(zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)) * utility_aggregated(price_ID_temp(zone_iter) + 1, zone_iter);
				utility_temp(zone_iter) /= bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter);				
			}
			if(price_ID_temp(zone_iter + 1) == -1){
				utility_temp(zone_iter + 1) = utility_aggregated(0, zone_iter + 1);
				utility_temp(zone_iter + 1) -= (quan_temp(zone_iter + 1) - bidded_total_aggregated(0, zone_iter + 1)) * (2 * TSO_Market.price_range_inflex(0) - TSO_Market.price_range_inflex(1));				
			}
			else if(price_ID_temp(zone_iter + 1) == TSO_Market.price_intervals + 2){
				utility_temp(zone_iter + 1) = utility_aggregated(TSO_Market.price_intervals + 2, zone_iter + 1);
				utility_temp(zone_iter + 1) -= (quan_temp(zone_iter + 1) - bidded_total_aggregated(TSO_Market.price_intervals + 2, zone_iter + 1)) * (2 * TSO_Market.price_range_inflex(1) - TSO_Market.price_range_inflex(0));				
			}
			else{
				utility_temp(zone_iter + 1) = (bidded_total_aggregated(price_ID_temp(zone_iter + 1) + 1, zone_iter + 1) - quan_temp(zone_iter + 1)) * utility_aggregated(price_ID_temp(zone_iter + 1), zone_iter + 1)
					+ (quan_temp(zone_iter + 1) - bidded_total_aggregated(price_ID_temp(zone_iter + 1), zone_iter + 1)) * utility_aggregated(price_ID_temp(zone_iter + 1) + 1, zone_iter + 1);
				utility_temp(zone_iter + 1) /= bidded_total_aggregated(price_ID_temp(zone_iter + 1) + 1, zone_iter + 1) - bidded_total_aggregated(price_ID_temp(zone_iter + 1), zone_iter + 1);
			}
			
			// Update source / sink, voltage, and power flow
			//std::cout << "Update source / sink, voltage, and power flow\n";
			voltage_temp.tail(TSO_Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(quan_temp.tail(TSO_Market.network.num_vertice - 1));
			flow_temp = (Problem.Constraint.eq_orig_matrix).topRightCorner(TSO_Market.network.num_edges, TSO_Market.network.num_vertice) * voltage_temp;		

			// Update errors
			//std::cout << "Update errors\n";
			error = 0.;
			// Power flow errors
			for(int edge_iter = 0; edge_iter < TSO_Market.network.num_edges; ++ edge_iter){
				error += pow(std::min(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 0), 0.) + std::max(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 1), 0.), 2.);
			}
			// Voltage errors
			for(int node_iter = 0; node_iter < TSO_Market.network.num_vertice; ++ node_iter){
				error += pow(std::min(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(TSO_Market.network.num_edges + TSO_Market.network.num_vertice + node_iter, 0), 0.) + std::max(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(TSO_Market.network.num_edges + TSO_Market.network.num_vertice + node_iter, 1), 0.), 2.);
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
		increment_coeff = TSO_Market.dof_metric.solve(increment_dot);
		grad = Eigen::VectorXd::Zero(TSO_Market.num_zone);
		for(int zone_iter = 0; zone_iter < TSO_Market.num_zone - 1; ++ zone_iter){
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
		for(int zone_iter = 0; zone_iter < TSO_Market.num_zone - 1; ++ zone_iter){
			//std::cout << "Update price after small increase / decrease\n";
			if(price_ID_temp(zone_iter) > -1 && price_ID_temp(zone_iter) < TSO_Market.price_intervals + 2){
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
						if(price_ID_temp(zone_iter) == TSO_Market.price_intervals + 2){
							break;
						}					
					}
				}				
			}
			else if(quan_temp(zone_iter) < bidded_total_aggregated(0, zone_iter)){
				price_ID_temp(zone_iter) = -1;
			}
			else if(quan_temp(zone_iter) > bidded_total_aggregated(TSO_Market.price_intervals + 2, zone_iter)){
				price_ID_temp(zone_iter) = TSO_Market.price_intervals + 2;
			}				
			if(price_ID_temp(zone_iter + 1) > -1 && price_ID_temp(zone_iter + 1) < TSO_Market.price_intervals + 2){
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
						if(price_ID_temp(zone_iter + 1) == TSO_Market.price_intervals + 2){
							break;
						}					
					}
				}				
			}
			else if(quan_temp(zone_iter + 1) < bidded_total_aggregated(0, zone_iter + 1)){
				price_ID_temp(zone_iter + 1) = -1;
			}
			else if(quan_temp(zone_iter + 1) > bidded_total_aggregated(TSO_Market.price_intervals + 2, zone_iter + 1)){
				price_ID_temp(zone_iter + 1) = TSO_Market.price_intervals + 2;
			}
			
			// Update utility function after small increase / decrease
			//std::cout << "Update utility function after small increase / decrease\n";
			if(price_ID_temp(zone_iter) == -1){
				utility_temp(zone_iter) = utility_aggregated(0, zone_iter);
				utility_temp(zone_iter) -= (quan_temp(zone_iter) - bidded_total_aggregated(0, zone_iter)) * (2 * TSO_Market.price_range_inflex(0) - TSO_Market.price_range_inflex(1));
			}
			else if(price_ID_temp(zone_iter) == TSO_Market.price_intervals + 2){
				utility_temp(zone_iter) = utility_aggregated(TSO_Market.price_intervals + 2, zone_iter);
				utility_temp(zone_iter) -= (quan_temp(zone_iter) - bidded_total_aggregated(TSO_Market.price_intervals + 2, zone_iter)) * (2 * TSO_Market.price_range_inflex(1) - TSO_Market.price_range_inflex(0));
			}
			else{
				utility_temp(zone_iter) = (bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - quan_temp(zone_iter)) * utility_aggregated(price_ID_temp(zone_iter), zone_iter)
					+ (quan_temp(zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)) * utility_aggregated(price_ID_temp(zone_iter) + 1, zone_iter);
				utility_temp(zone_iter) /= bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter);				
			}
			if(price_ID_temp(zone_iter + 1) == -1){
				utility_temp(zone_iter + 1) = utility_aggregated(0, zone_iter + 1);
				utility_temp(zone_iter + 1) -= (quan_temp(zone_iter + 1) - bidded_total_aggregated(0, zone_iter + 1)) * (2 * TSO_Market.price_range_inflex(0) - TSO_Market.price_range_inflex(1));				
			}
			else if(price_ID_temp(zone_iter + 1) == TSO_Market.price_intervals + 2){
				utility_temp(zone_iter + 1) = utility_aggregated(TSO_Market.price_intervals + 2, zone_iter + 1);
				utility_temp(zone_iter + 1) -= (quan_temp(zone_iter + 1) - bidded_total_aggregated(TSO_Market.price_intervals + 2, zone_iter + 1)) * (2 * TSO_Market.price_range_inflex(1) - TSO_Market.price_range_inflex(0));				
			}
			else{
				utility_temp(zone_iter + 1) = (bidded_total_aggregated(price_ID_temp(zone_iter + 1) + 1, zone_iter + 1) - quan_temp(zone_iter + 1)) * utility_aggregated(price_ID_temp(zone_iter + 1), zone_iter + 1)
					+ (quan_temp(zone_iter + 1) - bidded_total_aggregated(price_ID_temp(zone_iter + 1), zone_iter + 1)) * utility_aggregated(price_ID_temp(zone_iter + 1) + 1, zone_iter + 1);
				utility_temp(zone_iter + 1) /= bidded_total_aggregated(price_ID_temp(zone_iter + 1) + 1, zone_iter + 1) - bidded_total_aggregated(price_ID_temp(zone_iter + 1), zone_iter + 1);
			}		
		}
		
		// Update source / sink, voltage, and power flow
		//std::cout << "Update source / sink, voltage, and power flow\n";
		voltage_temp.tail(TSO_Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(quan_temp.tail(TSO_Market.network.num_vertice - 1));
		flow_temp = (Problem.Constraint.eq_orig_matrix).topRightCorner(TSO_Market.network.num_edges, TSO_Market.network.num_vertice) * voltage_temp;		

		// Update errors
		//std::cout << "Update errors\n";
		error = 0.;
		// Power flow errors
		for(int edge_iter = 0; edge_iter < TSO_Market.network.num_edges; ++ edge_iter){
			error += pow(std::min(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 0), 0.) + std::max(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 1), 0.), 2.);
		}
		// Voltage errors
		for(int node_iter = 0; node_iter < TSO_Market.network.num_vertice; ++ node_iter){
			error += pow(std::min(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(TSO_Market.network.num_edges + TSO_Market.network.num_vertice + node_iter, 0), 0.) + std::max(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(TSO_Market.network.num_edges + TSO_Market.network.num_vertice + node_iter, 1), 0.), 2.);
		}
		
		// Update objective
		//std::cout << "Update objective\n";
		obj_temp = (1. - mu) * utility_temp.sum() - mu * error;
		
		// Check if objective function is actually improved
		//std::cout << "Check if objective function is actually improved\n";
		if(obj_temp >= obj){
			quan = quan_temp;
			price_ID = price_ID_temp;
			utility = utility_temp;
			obj = obj_temp;
			//std::cout << "Objective function updated\n";
			//std::cout << price_ID.transpose() << "\n";
			//std::cout << quan.transpose() << "\n";
			//std::cout << obj << "\n\n";
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
		
		if(dS > .001){
			break_flag = 0;
		}
	}

	// Update final source / sink, voltage, and power flow
	//std::cout << "Update final source / sink, voltage, and power flow\n";
	Problem.Solution.orig_vector.segment(TSO_Market.network.num_edges, TSO_Market.network.num_vertice) = quan;
	Problem.Solution.orig_vector.tail(TSO_Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(Problem.Solution.orig_vector.segment(TSO_Market.network.num_edges + 1, TSO_Market.network.num_vertice - 1));
	Problem.Solution.orig_vector.head(TSO_Market.network.num_edges) = (Problem.Constraint.eq_orig_matrix).topRightCorner(TSO_Market.network.num_edges, TSO_Market.network.num_vertice) * Problem.Solution.orig_vector.tail(TSO_Market.network.num_vertice);		
	std::cout << Problem.Solution.orig_vector.segment(TSO_Market.network.num_edges, TSO_Market.network.num_vertice).transpose() << "\n\n";
	std::cout << Problem.Solution.orig_vector.transpose() << "\n\n";
	std::cout << quan.minCoeff() << " " << quan.maxCoeff() << "\n";
	std::cout << (Problem.Solution.orig_vector.tail(TSO_Market.network.num_vertice - 1)).minCoeff() << " " << (Problem.Solution.orig_vector.tail(TSO_Market.network.num_vertice - 1)).maxCoeff() << "\n";
	std::cout << (Problem.Solution.orig_vector.head(TSO_Market.network.num_edges)).minCoeff() << " " << (Problem.Solution.orig_vector.head(TSO_Market.network.num_edges)).maxCoeff() << "\n\n";
}

int main(){
	market_inform TSO_Market;
	TSO_Market_Set_Test_1(TSO_Market, 1);
	LP_object TSO_Problem;
	TSO_LP_Set(TSO_Market, TSO_Problem);

//	auto fin_node = "../../power_network/input/transmission_nodes.csv";
//	auto fin_edge = "../../power_network/input/transmission_edges_pu_simp.csv";
//	auto fin_pu_dc = "../../power_network/input/transmission_pu_dc.csv";
//	market_inform TSO_Market; 
//	TSO_Market_Set(TSO_Market, 1, fin_node, fin_edge, fin_pu_dc);
//	LP_object TSO_Problem;
//	TSO_LP_Set(TSO_Market, TSO_Problem);
	
	TSO_Market_Optimization(0, TSO_Market, TSO_Problem);
}