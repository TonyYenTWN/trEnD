// Source File for re-dispatch and tertiary control reserve market clearing of TSO in Norway
#include <iostream>
//#include <chrono>
//#include "../../basic/LP_gpa.cpp"
#include "../../basic/LP_gpa_fast.cpp"
#include "../power_market.cpp"

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
	TSO_Market.network.admittance_vector << 100., 100., 100.;
	
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
	TSO_Market.num_zone = 4;
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
	TSO_Market.network.power_constraint.col(0) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, -1000.);
	TSO_Market.network.power_constraint.col(1) = Eigen::VectorXd::Constant(TSO_Market.network.num_edges, 1000.);

	// Initialization of output variables
	TSO_Market.confirmed_supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed_price = Eigen::MatrixXd(Time, TSO_Market.num_zone);
	TSO_Market.network.confirmed_power = Eigen::MatrixXd(Time, TSO_Market.network.num_edges);
	
	// For the trivial case only: initialize submitted supply and demand bids at each node
	TSO_Market.submitted_supply = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	TSO_Market.submitted_demand = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	TSO_Market.submitted_supply.leftCols(TSO_Market.num_zone / 2) = Eigen::VectorXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 1.);
	TSO_Market.submitted_demand.rightCols(TSO_Market.num_zone / 2) = Eigen::VectorXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 1.);
	
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

	// Trivial initialization at the nodes
	TSO_Market.submitted_supply = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	TSO_Market.submitted_demand = Eigen::MatrixXd::Zero(TSO_Market.price_intervals + 2, TSO_Market.num_zone);
	TSO_Market.submitted_supply.leftCols(TSO_Market.num_zone / 2) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 1.);
	TSO_Market.submitted_demand.leftCols(TSO_Market.num_zone / 2) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 2.);
	TSO_Market.submitted_demand.rightCols(TSO_Market.num_zone / 2) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 1.);
	TSO_Market.submitted_supply.rightCols(TSO_Market.num_zone / 2) = Eigen::MatrixXd::Constant(TSO_Market.price_intervals + 2, TSO_Market.num_zone / 2, 2.);
//	TSO_Market.submitted_supply.row(0).head(100) = Eigen::VectorXd::Constant(100, 1.);
//	TSO_Market.submitted_supply.row(1).head(100) = Eigen::VectorXd::Constant(100, 1.);
//	TSO_Market.submitted_supply.row(2).head(100) = Eigen::VectorXd::Constant(100, 1.);
//	TSO_Market.submitted_demand.row(TSO_Market.price_intervals + 1).tail(100) = Eigen::VectorXd::Constant(100, 3.);
	
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
	// LP_process(Problem, "Linear Problem", 0, 0, 1);
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
//	std::cout << price_ID.transpose() << "\n\n";
//	std::cout << utility_aggregated << "\n\n";
//	std::cout << (bidded_supply_aggregated + bidded_demand_aggregated).topRows(10) << "\n\n";
//	std::cout << (bidded_supply_aggregated + bidded_demand_aggregated).bottomRows(10) << "\n\n";
//	std::cout << (bidded_demand_import + bidded_supply_import).bottomRows(10) << "\n\n";

	// Declare variables for the main loop
	double tol = pow(10., -12.);
	double eps = 1.;
	double mu = 1. - eps;
	double obj_max = -std::numeric_limits<double>::infinity();
	double obj_temp_mixed = 0.;
	double obj_temp_orig = 0.;
	Eigen::Vector2i direction_ID;
	Eigen::VectorXi price_ID_temp = price_ID;
	Eigen::VectorXd error(Problem.Constraints_eq_num + Problem.Variables_num);
	Eigen::MatrixXd Boundary_gap(Problem.Variables_num + Problem.Constraints_ie_num, 2);
	
	int loop_count = 0;
	while(eps > .001){
		while(1){
			loop_count += 1;
			//std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
			//std::cout << "Main Loop: " << loop_count << std::endl;
			//std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
					
			for(int zone_iter = 0; zone_iter < TSO_Market.num_zone; ++ zone_iter){
				//std::cout << "---------------------------------------------------------------------------" << std::endl;
				//std::cout << "New loop" << std::endl;
				//std::cout << "---------------------------------------------------------------------------" << std::endl;			
				// Increase price of bidding zone
				if(price_ID_temp(zone_iter) < TSO_Market.price_intervals + 2){
					obj_temp_orig -= utility_aggregated(price_ID_temp(zone_iter), zone_iter);
					price_ID_temp(zone_iter) = price_ID(zone_iter) + 1;
					while(bidded_supply_export(price_ID_temp(zone_iter), zone_iter) + bidded_demand_export(price_ID_temp(zone_iter), zone_iter) + bidded_supply_import(price_ID_temp(zone_iter), zone_iter) + bidded_demand_import(price_ID_temp(zone_iter), zone_iter) < tol && price_ID_temp(zone_iter) < TSO_Market.price_intervals + 2){
						price_ID_temp(zone_iter) += 1;
					}				
					Problem.Solution.orig_vector(TSO_Market.network.num_edges + zone_iter) = bidded_supply_aggregated(price_ID_temp(zone_iter), zone_iter) + bidded_demand_aggregated(price_ID_temp(zone_iter), zone_iter);
					obj_temp_orig += utility_aggregated(price_ID_temp(zone_iter), zone_iter);
					//std::cout << Problem.Solution.orig_vector.transpose() << "\n";
					
					// Update voltage
					Problem.Solver.ldlt.compute((Problem.Constraint.eq_orig_matrix).block(TSO_Market.network.num_edges + 1, TSO_Market.network.num_edges + TSO_Market.network.num_vertice + 1, TSO_Market.network.num_vertice - 1, TSO_Market.network.num_vertice - 1));
					Problem.Solution.orig_vector.tail(TSO_Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(Problem.Solution.orig_vector.segment(TSO_Market.network.num_edges + 1, TSO_Market.network.num_vertice - 1));
					//std::cout << (Problem.Constraint.eq_orig_matrix).block(TSO_Market.network.num_edges + 1, TSO_Market.network.num_edges + TSO_Market.network.num_vertice + 1, TSO_Market.network.num_vertice - 1, TSO_Market.network.num_vertice - 1) << "\n\n";
					//std::cout << Problem.Constraint.eq_orig_matrix << "\n\n";
					
					// Update power flow
					Problem.Solution.orig_vector.head(TSO_Market.network.num_edges) = (Problem.Constraint.eq_orig_matrix).topRightCorner(TSO_Market.network.num_edges, TSO_Market.network.num_vertice) * Problem.Solution.orig_vector.tail(TSO_Market.network.num_vertice);
					
					// Update error
					error = Eigen::VectorXd::Zero(Problem.Constraints_eq_num + Problem.Variables_num);
					error.head(Problem.Constraints_eq_num) = Problem.Constraint.eq_orig_matrix * Problem.Solution.orig_vector - Problem.Boundary.eq_vector;
					Boundary_gap.col(0) = Problem.Constraint.ie_orig_matrix * Problem.Solution.orig_vector - Problem.Boundary.ie_orig_matrix.col(0);
					Boundary_gap.col(1) = Problem.Boundary.ie_orig_matrix.col(1) - Problem.Constraint.ie_orig_matrix * Problem.Solution.orig_vector;
					for(int constraint_iter = 0; constraint_iter < Boundary_gap.rows(); ++ constraint_iter){
						error(Problem.Constraints_eq_num + constraint_iter) = -std::min(0., std::min(Boundary_gap(constraint_iter, 0), Boundary_gap(constraint_iter, 1)));
					}
					//std::cout << error.transpose() << "\n\n";
					
					// Rewrite current optimal direction, if applicable
					obj_temp_mixed = (1 - mu) * obj_temp_orig - mu * error.squaredNorm();
					if(obj_temp_mixed > obj_max){
						obj_max = obj_temp_mixed;
						direction_ID << zone_iter, 1;
					}
					//std::cout << (1 - mu) * utility_aggregated(price_ID_temp(zone_iter), zone_iter) << "\n\n";
					//std::cout << -mu * error.squaredNorm() << "\n\n";
					
					// Return price and objective to origin
					obj_temp_orig -= utility_aggregated(price_ID_temp(zone_iter), zone_iter);
					price_ID_temp(zone_iter) = price_ID(zone_iter);
					obj_temp_orig += utility_aggregated(price_ID_temp(zone_iter), zone_iter);									
				}
				
				// Decrease price of bidding zone
				if(price_ID_temp(zone_iter) > 0){
					obj_temp_orig -= utility_aggregated(price_ID_temp(zone_iter), zone_iter);
					price_ID_temp(zone_iter) = price_ID(zone_iter) - 1;
					while(bidded_supply_export(price_ID_temp(zone_iter), zone_iter) + bidded_demand_export(price_ID_temp(zone_iter), zone_iter) + bidded_supply_import(price_ID_temp(zone_iter), zone_iter) + bidded_demand_import(price_ID_temp(zone_iter), zone_iter) < tol && price_ID_temp(zone_iter) > 0){
						price_ID_temp(zone_iter) -= 1;
	//					std::cout << bidded_supply_export(price_ID_temp(zone_iter), zone_iter) + bidded_demand_export(price_ID_temp(zone_iter), zone_iter) << "\n\n"
	//					std::cout << price_ID_temp(zone_iter) << "\n\n";
	//					std::cout << "loop decrease" << "\n\n";
					}
					Problem.Solution.orig_vector(TSO_Market.network.num_edges + zone_iter) = bidded_supply_aggregated(price_ID_temp(zone_iter), zone_iter) + bidded_demand_aggregated(price_ID_temp(zone_iter), zone_iter);
					obj_temp_orig += utility_aggregated(price_ID_temp(zone_iter), zone_iter);
					//std::cout << Problem.Solution.orig_vector.transpose() << "\n";
					
					// Update voltage
					Problem.Solver.ldlt.compute((Problem.Constraint.eq_orig_matrix).block(TSO_Market.network.num_edges + 1, TSO_Market.network.num_edges + TSO_Market.network.num_vertice + 1, TSO_Market.network.num_vertice - 1, TSO_Market.network.num_vertice - 1));
					Problem.Solution.orig_vector.tail(TSO_Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(Problem.Solution.orig_vector.segment(TSO_Market.network.num_edges + 1, TSO_Market.network.num_vertice - 1));
			
					// Update power flow
					Problem.Solution.orig_vector.head(TSO_Market.network.num_edges) = (Problem.Constraint.eq_orig_matrix).topRightCorner(TSO_Market.network.num_edges, TSO_Market.network.num_vertice) * Problem.Solution.orig_vector.tail(TSO_Market.network.num_vertice);
			
					// Update error
					error = Eigen::VectorXd::Zero(Problem.Constraints_eq_num + Problem.Variables_num);
					error.head(Problem.Constraints_eq_num) = Problem.Constraint.eq_orig_matrix * Problem.Solution.orig_vector - Problem.Boundary.eq_vector;
					Boundary_gap.col(0) = Problem.Constraint.ie_orig_matrix * Problem.Solution.orig_vector - Problem.Boundary.ie_orig_matrix.col(0);
					Boundary_gap.col(1) = Problem.Boundary.ie_orig_matrix.col(1) - Problem.Constraint.ie_orig_matrix * Problem.Solution.orig_vector;
					for(int constraint_iter = 0; constraint_iter < Boundary_gap.rows(); ++ constraint_iter){
						error(Problem.Constraints_eq_num + constraint_iter) = -std::min(0., std::min(Boundary_gap(constraint_iter, 0), Boundary_gap(constraint_iter, 1)));
					}
					//std::cout << error.transpose() << "\n\n";
					
					// Rewrite current optimal direction, if applicable
					obj_temp_mixed = (1. - mu) * obj_temp_orig - mu * error.squaredNorm();
					if(obj_temp_mixed > obj_max){
						obj_max = obj_temp_mixed;
						direction_ID << zone_iter, -1;
					}
					//std::cout << (1 - mu) * utility_aggregated(price_ID_temp(zone_iter), zone_iter) << "\n\n";
					//std::cout << -mu * error.squaredNorm() << "\n\n";
					
					// Return price and objective to origin
					obj_temp_orig -= utility_aggregated(price_ID_temp(zone_iter), zone_iter);
					price_ID_temp(zone_iter) = price_ID(zone_iter);
					obj_temp_orig += utility_aggregated(price_ID_temp(zone_iter), zone_iter);					
				}		
			}
			
			// Update optimal direction and default objective
			if(direction_ID(1) == 0){
				break;
			}
			
			switch(direction_ID(1)){
				case -1:
					obj_temp_orig -= utility_aggregated(price_ID(direction_ID(0)), direction_ID(0));
					price_ID(direction_ID(0)) -= 1;
					obj_temp_orig += utility_aggregated(price_ID(direction_ID(0)), direction_ID(0));
					break;
				case 1:
					obj_temp_orig -= utility_aggregated(price_ID(direction_ID(0)), direction_ID(0));
					price_ID(direction_ID(0)) += 1;
					obj_temp_orig += utility_aggregated(price_ID(direction_ID(0)), direction_ID(0));
					break;
				case 0:
					break;
			}
			price_ID_temp = price_ID;
			direction_ID(1) = 0;		
		}
		
		eps /= 2.;
		mu = 1. - eps;
	}
	
	std::cout << Problem.Solution.orig_vector.transpose() << "\n";	
}

int main(){
	market_inform TSO_Market = TSO_Market_Set_Test_1(1);
	LP_object TSO_Problem;
	TSO_LP_Set(TSO_Market, TSO_Problem);

//	auto fin_node = "../power_network/input/transmission_nodes.csv";
//	auto fin_edge = "../power_network/input/transmission_edges_pu_simp.csv";
//	auto fin_pu_dc = "../power_network/input/transmission_pu_dc.csv";
//	auto TSO_Market = TSO_Market_Set(1, fin_node, fin_edge, fin_pu_dc);
//	LP_object TSO_Problem;
//	TSO_LP_Set(TSO_Market, TSO_Problem);
	
	TSO_Market_Optimization(0, TSO_Market, TSO_Problem);
}