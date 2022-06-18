// Source File for re-dispatch and tertiary control reserve market clearing of TSO in Norway
#include <iostream>
//#include <chrono>
#include "../basic/LP_gpa.cpp"
#include "power_market.cpp"

market_inform TSO_Market_Set_Test_1(int Time){
	market_inform TSO_Market;
	
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
	
	// Set node admittance matrix Y_n
	TSO_Market.network.num_vertice = TSO_Market.num_zone;
	TSO_Market.network.num_edges = 3;
	TSO_Market.network.incidence_matrix = Eigen::MatrixXi(TSO_Market.network.num_edges, 2);
	TSO_Market.network.incidence_matrix.row(0) << 0, 1;
	TSO_Market.network.incidence_matrix.row(1) << 0, 2;
	TSO_Market.network.incidence_matrix.row(2) << 1, 2;
	TSO_Market.network.admittance_vector = Eigen::VectorXd(3);
	TSO_Market.network.admittance_vector << 100., 100., 50.;
	
	// Set voltage and power constraints at each edge
	TSO_Market.network.voltage_constraint = Eigen::MatrixXd(TSO_Market.network.num_vertice, 2);
	TSO_Market.network.voltage_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, -.1);
	TSO_Market.network.voltage_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, .1);
	TSO_Market.network.power_constraint = Eigen::MatrixXd(TSO_Market.network.num_edges, 2);
	TSO_Market.network.power_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, -50.);
	TSO_Market.network.power_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, 50.);

	// Initialization of output variables
	TSO_Market.confirmed_supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_price = Eigen::MatrixXd(Time, TSO_Market.num_zone);
	TSO_Market.network.confirmed_power = Eigen::MatrixXd(Time, TSO_Market.network.num_edges);
	
	// For the trivial case only: initialize submitted supply and demand bids at each node
	TSO_Market.submitted_supply = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	TSO_Market.submitted_demand = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
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
	
	return(TSO_Market);
}

market_inform TSO_Market_Set_Test_2(int Time){
	market_inform TSO_Market;
	
	// Input parameters of TSO market
	// A trivial test case with 20 nodes connected as a radial line
	TSO_Market.num_zone = 6;
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
	TSO_Market.network.admittance_vector = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, 2000.);
	
	// Set voltage and power constraints at each edge
	TSO_Market.network.voltage_constraint = Eigen::MatrixXd(TSO_Market.network.num_vertice, 2);
	TSO_Market.network.voltage_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, -1.);
	TSO_Market.network.voltage_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_vertice, 1.);
	TSO_Market.network.power_constraint = Eigen::MatrixXd(TSO_Market.network.num_edges, 2);
	TSO_Market.network.power_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, -50.);
	TSO_Market.network.power_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, 50.);

	// Initialization of output variables
	TSO_Market.confirmed_supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_price = Eigen::MatrixXd(Time, TSO_Market.num_zone);
	TSO_Market.network.confirmed_power = Eigen::MatrixXd(Time, TSO_Market.network.num_edges);
	
	// For the trivial case only: initialize submitted supply and demand bids at each node
	TSO_Market.submitted_supply = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	TSO_Market.submitted_demand = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	TSO_Market.submitted_supply.row(0).head(TSO_Market.num_zone / 2) = Eigen::VectorXd::Constant(TSO_Market.num_zone / 2, 10.);
	TSO_Market.submitted_demand.row(TSO_Market.price_intervals + 1).tail(TSO_Market.num_zone / 2) = Eigen::VectorXd::Constant(TSO_Market.num_zone / 2, 10.);
	
	return(TSO_Market);
}

double impedence_conversion(Eigen::MatrixXd pu_dc_inform, double voltage){
	int v_iter = 0;
	while(pu_dc_inform(v_iter, 0) != voltage){
		v_iter += 1;
	}
	
	return(pu_dc_inform(v_iter, 1));
}

market_inform TSO_Market_Set(int Time, std::string fin_node, std::string fin_edge, std::string fin_pu_dc){
	market_inform TSO_Market;
	
	// Read power network data
	auto fin_node_dim = get_file_dim(fin_node);
	auto fin_edge_dim = get_file_dim(fin_edge);
	auto fin_pu_dc_dim = get_file_dim(fin_pu_dc);
	auto node_inform = read_file(fin_node_dim[0], fin_node_dim[1], fin_node);
	auto edge_inform = read_file(fin_edge_dim[0], fin_edge_dim[1], fin_edge);
	auto pu_dc_inform = read_file(fin_pu_dc_dim[0], fin_pu_dc_dim[1], fin_pu_dc);
	double x_L = 5. * pow(10., -4.);		// Inductance per meter of transmission line
	
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
		TSO_Market.network.admittance_vector(edge_iter) = impedence_conversion(pu_dc_inform, edge_inform(edge_iter, 4)) / x_L / edge_inform(edge_iter, 5);
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

	// Trivial initialization at the nodes
	TSO_Market.submitted_supply = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	TSO_Market.submitted_demand = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	TSO_Market.submitted_supply.row(0).head(10) = Eigen::VectorXd::Constant(10, 10.);
	TSO_Market.submitted_demand.row(TSO_Market.price_intervals + 1).tail(10) = Eigen::VectorXd::Constant(10, 10.);
	
	return(TSO_Market);
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
	
	// Set sparse matrix for original inequality constraints
	Problem.Constraint.ie_orig_matrix = Eigen::SparseMatrix <double> (Problem.Variables_num + Problem.Constraints_ie_num, Problem.Variables_num);
	std::vector<Trip> Constraint_ie_trip;
	Constraint_ie_trip.reserve((Problem.Constraints_ie_num + 1) * Problem.Variables_num);
	for(int var_iter = 0; var_iter < Problem.Variables_num; ++ var_iter){
		Constraint_ie_trip.push_back(Trip(var_iter, var_iter, 1));
	}
	Problem.Constraint.ie_orig_matrix.setFromTriplets(Constraint_ie_trip.begin(), Constraint_ie_trip.end());
	
	// Set initial feasible solution
	Problem.Solution.orig_vector = Eigen::VectorXd::Zero(Problem.Variables_num);
	
	// Initialize the LP problem solver
	LP_process(Problem, "Linear Problem", 0, 0, 1);
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
	//#pragma omp parallel
	{
		//#pragma omp for
		for(int zone_ID = 0; zone_ID < TSO_Market.num_zone; ++ zone_ID){
			// Supply curves
			bidded_supply_export.col(zone_ID).array().tail(TSO_Market.price_intervals + 1 - price_ID(zone_ID)) = TSO_Market.submitted_supply.col(zone_ID).array().tail(TSO_Market.price_intervals + 1 - price_ID(zone_ID));
			bidded_supply_export(price_ID(zone_ID), zone_ID) = bidded_supply(price_ID(zone_ID), zone_ID);
			bidded_supply_import.col(zone_ID).array().head(price_ID(zone_ID)) = TSO_Market.submitted_supply.col(zone_ID).array().head(price_ID(zone_ID));
			bidded_supply_import(price_ID(zone_ID), zone_ID) = TSO_Market.submitted_supply(price_ID(zone_ID), zone_ID) - bidded_supply(price_ID(zone_ID), zone_ID);
			
			// Demand curves
			bidded_demand_export.col(zone_ID).array().tail(TSO_Market.price_intervals + 1 - price_ID(zone_ID)) = TSO_Market.submitted_demand.col(zone_ID).array().tail(TSO_Market.price_intervals + 1 - price_ID(zone_ID));
			bidded_demand_export(price_ID(zone_ID), zone_ID) = TSO_Market.submitted_demand(price_ID(zone_ID), zone_ID) - bidded_demand(price_ID(zone_ID), zone_ID);
			bidded_demand_import.col(zone_ID).array().head(price_ID(zone_ID)) = TSO_Market.submitted_demand.col(zone_ID).array().head(price_ID(zone_ID));
			bidded_demand_import(price_ID(zone_ID), zone_ID) = bidded_demand(price_ID(zone_ID), zone_ID);
			
			//std::cout << price_ID(zone_ID) << " " << bidded_supply_import(price_ID(zone_ID), zone_ID) << " " << bidded_demand_import(price_ID(zone_ID), zone_ID) << "\n";
			//std::cout << price_ID(zone_ID) << " " << bidded_supply_export(price_ID(zone_ID), zone_ID) << " " << bidded_demand_export(price_ID(zone_ID), zone_ID) << "\n\n";		
		}
	}
	
	// Declare variables for the main loop
	double eps = pow(10, -8);
	double tol = pow(10, -12);
	double ratio; 
	double Improvement_Obj;
	Eigen::MatrixXd Updated_Boundary = Problem.Constraint.permutation_matrix * Problem.Boundary.ie_orig_matrix;
	Eigen::VectorXd Updated_Objective = Eigen::VectorXd::Zero(Problem.Variables_num);
	Eigen::VectorXd Previous_Sol = Problem.Solution.orig_vector;
	Eigen::VectorXd trade_quantity;
	
	// Initialize objective values and coefficients and inequality boundary values for the LP problem
	#pragma omp parallel
	{
		#pragma omp for
		for(int node_iter = 0; node_iter < TSO_Market.network.num_vertice; ++ node_iter){
			Updated_Objective(TSO_Market.network.num_edges + node_iter) = -TSO_Market.bidded_price(price_ID(node_iter));
			Updated_Boundary(TSO_Market.network.num_edges + node_iter, 1) = bidded_supply_export(price_ID(node_iter), node_iter) + bidded_demand_export(price_ID(node_iter), node_iter);
			Updated_Boundary(TSO_Market.network.num_edges + node_iter, 0) = -(bidded_supply_import(price_ID(node_iter), node_iter) + bidded_demand_import(price_ID(node_iter), node_iter));
		}
	}
	Problem.Objective.orig_vector = Problem.Constraint.permutation_matrix.transpose() * Updated_Objective;
	Problem.Boundary.ie_orig_matrix = Problem.Constraint.permutation_matrix.transpose() * Updated_Boundary;
	Problem.Objective.orig_value_sum = Problem.Solution.orig_vector.dot(Problem.Objective.orig_vector);
	//std::cout << Problem.Objective.orig_vector.transpose() << "\n\n";
	
	// Main loop
	int loop_count = 0;
	//while(loop_count < 5){
	while(1){
		loop_count += 1;
		
		// Test optimization
		//LP_process(Problem, "TSO Problem", 1, 1, 0, 0, 1, 1);
		//break;
		
		// One iteration of optimization
		LP_process(Problem, "TSO Problem", 1, 1, 1, 0, 1, 1);
		//std::cout << Problem.Objective.update_coeff.transpose() << "\n\n";
		
		// Update equality boundary values for the LP problem and confirmed prices and bids for the market
		trade_quantity = Problem.Solution.orig_vector - Previous_Sol;
		Improvement_Obj = trade_quantity.dot(Problem.Constraint.permutation_matrix * Problem.Objective.orig_vector);
		trade_quantity = trade_quantity.array().abs();
		for(int node_iter = 0; node_iter < TSO_Market.network.num_vertice; ++ node_iter){
			// Update confirmed supply and demand bids at each node
			// Check whether the node is importing or exporting electricity
			if(Problem.Solution.orig_vector(TSO_Market.network.num_edges + node_iter) <= Previous_Sol(TSO_Market.network.num_edges + node_iter)){
				if(bidded_supply_import(price_ID(node_iter), node_iter) + bidded_demand_import(price_ID(node_iter), node_iter) > tol){
					ratio = bidded_supply_import(price_ID(node_iter), node_iter) / (bidded_supply_import(price_ID(node_iter), node_iter) + bidded_demand_import(price_ID(node_iter), node_iter));
					TSO_Market.confirmed_supply(tick, node_iter) -= ratio * trade_quantity(TSO_Market.network.num_edges + node_iter);
					TSO_Market.confirmed_demand(tick, node_iter) += (1 - ratio) * trade_quantity(TSO_Market.network.num_edges + node_iter);
					bidded_supply_import(price_ID(node_iter), node_iter) -= ratio * trade_quantity(TSO_Market.network.num_edges + node_iter);
					bidded_demand_import(price_ID(node_iter), node_iter) -= (1 - ratio) * trade_quantity(TSO_Market.network.num_edges + node_iter);					
				}
				else{
					bidded_supply_import(price_ID(node_iter), node_iter) = 0.;
					bidded_demand_import(price_ID(node_iter), node_iter) = 0.;
				}
//				std::cout << price_ID(node_iter) << " " << bidded_supply_import(price_ID(node_iter), node_iter) << " " << bidded_demand_import(price_ID(node_iter), node_iter) << "\n\n";	
			}
			else{
				if(bidded_supply_export(price_ID(node_iter), node_iter) + bidded_demand_export(price_ID(node_iter), node_iter) > tol){
					ratio = bidded_supply_export(price_ID(node_iter), node_iter) / (bidded_supply_export(price_ID(node_iter), node_iter) + bidded_demand_export(price_ID(node_iter), node_iter));
					TSO_Market.confirmed_supply(tick, node_iter) += ratio * trade_quantity(TSO_Market.network.num_edges + node_iter);
					TSO_Market.confirmed_demand(tick, node_iter) -= (1 - ratio) * trade_quantity(TSO_Market.network.num_edges + node_iter);
					bidded_supply_export(price_ID(node_iter), node_iter) -= ratio * trade_quantity(TSO_Market.network.num_edges + node_iter);
					bidded_demand_export(price_ID(node_iter), node_iter) -= (1 - ratio) * trade_quantity(TSO_Market.network.num_edges + node_iter);					
				}
				else{
					bidded_supply_export(price_ID(node_iter), node_iter) = 0.;
					bidded_demand_export(price_ID(node_iter), node_iter) = 0.;
				}
//				std::cout << price_ID(node_iter) << " " << bidded_supply_export(price_ID(node_iter), node_iter) << " " << bidded_demand_export(price_ID(node_iter), node_iter) << "\n\n";
			}
			
			// Update price_ID and boundary at each node
			if(Problem.Objective.update_coeff(TSO_Market.network.num_edges + node_iter) == -1.){
				if(price_ID(node_iter) > 0){
					while(bidded_supply_import(price_ID(node_iter), node_iter) + bidded_demand_import(price_ID(node_iter), node_iter) < tol && price_ID(node_iter) > 0){
						price_ID(node_iter) -= 1;
					}
				}
				Updated_Boundary(TSO_Market.network.num_edges + node_iter, 1) = Updated_Boundary(TSO_Market.network.num_edges + node_iter, 0);
				Updated_Boundary(TSO_Market.network.num_edges + node_iter, 0) -= bidded_supply_import(price_ID(node_iter), node_iter) + bidded_demand_import(price_ID(node_iter), node_iter);
			}
			else if(Problem.Objective.update_coeff(TSO_Market.network.num_edges + node_iter) == 1.){
				if(price_ID(node_iter) < TSO_Market.price_intervals + 1){
					while(bidded_supply_export(price_ID(node_iter), node_iter) + bidded_demand_export(price_ID(node_iter), node_iter) < tol && price_ID(node_iter) < TSO_Market.price_intervals + 1){
						price_ID(node_iter) += 1;
					}				
				}
				Updated_Boundary(TSO_Market.network.num_edges + node_iter, 0) = Updated_Boundary(TSO_Market.network.num_edges + node_iter, 1);				
				Updated_Boundary(TSO_Market.network.num_edges + node_iter, 1) += bidded_supply_export(price_ID(node_iter), node_iter) + bidded_demand_export(price_ID(node_iter), node_iter);			
			}
			Updated_Objective(TSO_Market.network.num_edges + node_iter) = -TSO_Market.bidded_price(price_ID(node_iter));						
		}
		
		Problem.Objective.orig_vector = Problem.Constraint.permutation_matrix.transpose() * Updated_Objective;
		Problem.Boundary.ie_orig_matrix = Problem.Constraint.permutation_matrix.transpose() * Updated_Boundary;
		
		// Check optimality
		if(Problem.Objective.orig_value > tol && Improvement_Obj > tol){
			Previous_Sol = Problem.Solution.orig_vector;
			Problem.Objective.orig_value_sum += Improvement_Obj;
		}
		else{
			std::cout << Problem.Solution.orig_vector.transpose() << "\n\n";
			std::cout << Problem.Objective.orig_value_sum << "\n\n";
			break;
		}		
	}
}

int main(){
	market_inform TSO_Market = TSO_Market_Set_Test_2(1);
	LP_object TSO_Problem;
	TSO_LP_Set(TSO_Market, TSO_Problem);

//	auto fin_node = "../power_network/input/transmission_nodes.csv";
//	auto fin_edge = "../power_network/input/transmission_edges.csv";
//	auto fin_pu_dc = "../power_network/input/transmission_pu_dc.csv";
//	auto TSO_Market = TSO_Market_Set(1, fin_node, fin_edge, fin_pu_dc);
//	LP_object TSO_Problem;
//	TSO_LP_Set(TSO_Market, TSO_Problem);
	
	TSO_Market_Optimization(0, TSO_Market, TSO_Problem);
}