// Source file for the functions of the power market
#include "power_market.h"

// ------------------------------------------------------------------------------------------------
// Generic functions for all market types
// ------------------------------------------------------------------------------------------------
void power_market::parameters::bidded_price(price_ID_bimap &obj){
	int price_intervals = price_interval();
	// Range of lowest and highest possible bidding prices.
	Eigen::Vector2d price_range_inflex = Eigen::Vector2d(-500., 3000.);
	// Range of bidding prices for flexible supply and demand in the model.
	Eigen::Vector2d price_range_flex = Eigen::Vector2d(-100., 500.);

	Eigen::VectorXd bidded_price = Eigen::VectorXd(price_intervals + 2);
	bidded_price(0) = price_range_inflex(0);
	bidded_price.array().tail(1) = price_range_inflex(1);
	bidded_price.segment(1, price_intervals) = Eigen::VectorXd::LinSpaced(price_intervals, price_range_flex(0) + .5 * (price_range_flex(1) - price_range_flex(0)) / price_intervals, price_range_flex(1) - .5 * (price_range_flex(1) - price_range_flex(0)) / price_intervals);

	std::map <double, int> price_ID;
	for(int price_iter = 0; price_iter < bidded_price.size(); ++ price_iter){
		price_ID.insert(std::pair <double, int> (bidded_price(price_iter), price_iter));
	}

	obj.bidded_price = bidded_price;
	obj.price_ID = price_ID;
}

void power_market::Market_Initialization(market_inform &Market){
	// Initialization of process variables
	// Should re-initialize for every time slice
	Market.submitted_supply = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Market.submitted_demand = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Market.reference_price = Eigen::VectorXd::Zero(Market.num_zone);
	Market.confirmed.ratio_supply = Eigen::VectorXd::Zero(Market.num_zone);
	Market.confirmed.ratio_demand = Eigen::VectorXd::Zero(Market.num_zone);
	Market.actual.ratio_supply = Eigen::VectorXd::Zero(Market.num_zone);
	Market.actual.ratio_demand = Eigen::VectorXd::Zero(Market.num_zone);
}

void power_market::Operation_Initialization(market_inform &Market, int Time){
    //int Time = configuration::parameters::Time();

	Market.operation.cross_border.balancing = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.end_user.balancing = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.industrial.balancing = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.hydro.balancing = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.slack.balancing = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.wind.balancing = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.pump_storage.balancing = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.cross_border.EOM = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.end_user.EOM = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.industrial.EOM = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.hydro.EOM = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.slack.EOM = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.wind.EOM = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.pump_storage.EOM = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.cross_border.redispatch = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.end_user.redispatch = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.industrial.redispatch = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.hydro.redispatch = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.slack.redispatch = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.wind.redispatch = Eigen::MatrixXd::Zero(Time, Market.num_zone);
	Market.operation.pump_storage.redispatch = Eigen::MatrixXd::Zero(Time, Market.num_zone);
}

// ------------------------------------------------------------------------------------------------
// Specific functions for for flow-based markets
// ------------------------------------------------------------------------------------------------
void power_market::Flow_Based_Market_LP_Set(market_inform &Market){
	// -------------------------------------------------------------------------------
	// LP Solver initialization for flow-based market optimization
	// Warm-up once and reuse for the rest of time slices
	// Variables are sorted as {V}, {S}, {S}_0, {S}_1, {S}_2, ..., {I}
	// {S}_i is the set of source / sink at node #i at different price levels
	// -------------------------------------------------------------------------------

	// -------------------------------------------------------------------------------
	// Set matrix for general constraints
	// -------------------------------------------------------------------------------
	// Construct node admittance matrix
	std::vector <Eigen::TripletXd> Y_n_trip;
	Y_n_trip.reserve(2 * Market.network.num_edges + Market.network.num_vertice);
	Eigen::SparseMatrix <double, Eigen::RowMajor> Y_n(Market.network.num_vertice, Market.network.num_vertice);
	Eigen::VectorXd Y_n_diag = Eigen::VectorXd::Zero(Market.network.num_vertice);
	Eigen::VectorXpd Connection_num = Eigen::VectorXpd::Ones(Market.network.num_vertice);
	for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
		double y_edge = Market.network.admittance[edge_iter];

		// Equality constraints of voltage - source / sink at the nodes, off-diagonal terms
		Y_n_trip.push_back(Eigen::TripletXd(Market.network.incidence[edge_iter](0), Market.network.incidence[edge_iter](1), -y_edge));
		Y_n_trip.push_back(Eigen::TripletXd(Market.network.incidence[edge_iter](1), Market.network.incidence[edge_iter](0), -y_edge));
		Connection_num(Market.network.incidence[edge_iter](0)) += 1;
		Connection_num(Market.network.incidence[edge_iter](1)) += 1;

		// Equality constraints of voltage - source / sink at the nodes, diagonal terms
		Y_n_diag(Market.network.incidence[edge_iter](0)) += y_edge;
		Y_n_diag(Market.network.incidence[edge_iter](1)) += y_edge;
	}
	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		// Equality constraints of voltage - source / sink at the nodes, summed diagonal terms
		Y_n_trip.push_back(Eigen::TripletXd(node_iter, node_iter, Y_n_diag(node_iter)));
	}
	Y_n.setFromTriplets(Y_n_trip.begin(), Y_n_trip.end());

	// Generate sparse matrix for general (equality) constraints of voltage, power flow, and source / sink summation
	int constrant_num = 2 * Market.network.num_vertice + Market.network.num_edges;
	int variable_num = Market.network.num_vertice * (Market.price_intervals + 4) + Market.network.num_edges;
	Eigen::VectorXpd non_zero_num(constrant_num);
	non_zero_num << Connection_num + Eigen::VectorXpd::Ones(Market.network.num_vertice), Eigen::VectorXpd::Constant(Market.network.num_edges, 3), Eigen::VectorXpd::Constant(Market.network.num_vertice, Market.price_intervals + 3);
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
		double y_edge = Market.network.admittance[edge_iter];

		if(Market.network.incidence[edge_iter](0) < Market.network.incidence[edge_iter](1)){
			alglib::sparseset(constraint_general, Y_n.rows() + edge_iter, Market.network.incidence[edge_iter](0), y_edge);
			alglib::sparseset(constraint_general, Y_n.rows() + edge_iter, Market.network.incidence[edge_iter](1), -y_edge);
		}
		else{
			alglib::sparseset(constraint_general, Y_n.rows() + edge_iter, Market.network.incidence[edge_iter](0), -y_edge);
			alglib::sparseset(constraint_general, Y_n.rows() + edge_iter, Market.network.incidence[edge_iter](1), y_edge);
		}
		alglib::sparseset(constraint_general, Y_n.rows() + edge_iter, Market.network.num_vertice * (Market.price_intervals + 4) + edge_iter, -1.);
	}

	// Rows for source / sink summation at each node
	for(int node_iter = 0; node_iter < Y_n.rows(); ++ node_iter){
		int row_ID = Y_n.rows() + Market.network.num_edges + node_iter;
		alglib::sparseset(constraint_general, row_ID, Y_n.rows() + node_iter, 1.);
		int col_start = 2 * Y_n.rows() + node_iter * (Market.price_intervals + 2);
		int col_end = 2 * Y_n.rows() + (node_iter + 1) * (Market.price_intervals + 2) - 1;
		for(int col_iter = col_start; col_iter <= col_end; ++ col_iter){
			alglib::sparseset(constraint_general, row_ID, col_iter, -1.);
		}
	}

	// Check if the sparse matrix is correct
//	std::cout << Y_n << "\n\n";
//	double value;
//	for(int row_iter = 0; row_iter < non_zero_num.size(); ++ row_iter){
//		for(int col_iter = 0; col_iter < 2 * Y_n.cols() + 3; ++ col_iter){
//			value = sparseget(constraint_general, row_iter, col_iter);
//			std::cout << value << "\t";
//		}
//		std::cout << "\t";
//		for(int col_iter = Market.network.num_vertice * (Market.price_intervals + 4); col_iter < Market.network.num_vertice * (Market.price_intervals + 4) + 1; ++ col_iter){
//			value = sparseget(constraint_general, row_iter, col_iter);
//			std::cout << value << "\t";
//		}
//		std::cout << "\n";
//	}

	// -------------------------------------------------------------------------------
	// Set bounds for general and box constraints
	// -------------------------------------------------------------------------------
	Eigen::MatrixXd bound_general = Eigen::MatrixXd::Zero(constrant_num, 2);
	Eigen::MatrixXd bound_box = Eigen::MatrixXd::Zero(variable_num, 2);
	bound_box.topRows(Market.network.num_vertice) = Market.network.voltage_constraint;
	bound_box.middleRows(Market.network.num_vertice, Market.network.num_vertice).col(0) = Eigen::VectorXd::Constant(Market.network.num_vertice, -std::numeric_limits<double>::infinity());
	bound_box.middleRows(Market.network.num_vertice, Market.network.num_vertice).col(1) = Eigen::VectorXd::Constant(Market.network.num_vertice, std::numeric_limits<double>::infinity());
	bound_box.bottomRows(Market.network.num_edges) = Market.network.power_constraint;

	// Bounds of general constraints
	alglib::real_1d_array lb_general;
	alglib::real_1d_array ub_general;
	lb_general.setcontent(bound_general.rows(), bound_general.col(0).data());
	ub_general.setcontent(bound_general.rows(), bound_general.col(1).data());

	// -------------------------------------------------------------------------------
	// Set scale of variables
	// -------------------------------------------------------------------------------
	Eigen::VectorXd scale_vec = Eigen::VectorXd::Ones(variable_num);
	scale_vec.head(Market.network.num_vertice) = bound_box.col(1).head(Market.network.num_vertice) - bound_box.col(0).head(Market.network.num_vertice);
	scale_vec.tail(Market.network.num_edges) = bound_box.col(1).tail(Market.network.num_edges) - bound_box.col(0).tail(Market.network.num_edges);
	alglib::real_1d_array scale;
	scale.setcontent(scale_vec.size(), scale_vec.data());

	// -------------------------------------------------------------------------------
	// Set objective coefficients of variables
	// -------------------------------------------------------------------------------
	Eigen::VectorXd obj_vec = Eigen::VectorXd::Zero(variable_num);
	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		int row_start = 2 * Market.network.num_vertice + node_iter * (Market.price_intervals + 2);
		obj_vec.segment(row_start, Market.price_intervals + 2) = Market.bidded_price_map.bidded_price;
	}
	alglib::real_1d_array obj_coeff;
	obj_coeff.setcontent(obj_vec.size(), obj_vec.data());

	// -------------------------------------------------------------------------------
	// Set the LP problem object
	// -------------------------------------------------------------------------------
	alglib::minlpcreate(variable_num, Market.Problem);
	alglib::minlpsetcost(Market.Problem, obj_coeff);
	alglib::minlpsetlc2(Market.Problem, constraint_general, lb_general, ub_general, constrant_num);
	alglib::minlpsetscale(Market.Problem, scale);
	alglib::minlpsetalgodss(Market.Problem, 0.);
}

void power_market::Flow_Based_Market_Optimization(market_inform &Market){
	// -------------------------------------------------------------------------------
	// Update bounds for box constraints
	// -------------------------------------------------------------------------------
	int variable_num = Market.network.num_vertice * (Market.price_intervals + 4) + Market.network.num_edges;
	Eigen::MatrixXd bound_box(variable_num, 2);
	bound_box.topRows(Market.network.num_vertice) = Market.network.voltage_constraint;
	bound_box.middleRows(Market.network.num_vertice, Market.network.num_vertice).col(0) = Eigen::VectorXd::Constant(Market.network.num_vertice, -std::numeric_limits<double>::infinity());
	bound_box.middleRows(Market.network.num_vertice, Market.network.num_vertice).col(1) = Eigen::VectorXd::Constant(Market.network.num_vertice, std::numeric_limits<double>::infinity());
	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		int row_start = 2 * Market.network.num_vertice + node_iter * (Market.price_intervals + 2);
		bound_box.middleRows(row_start, Market.price_intervals + 2).col(0) = -Market.submitted_demand.col(node_iter);
		bound_box.middleRows(row_start, Market.price_intervals + 2).col(1) = Market.submitted_supply.col(node_iter);
	}
	bound_box.bottomRows(Market.network.num_edges) = Market.network.power_constraint;

	// Bounds of box constraints
	alglib::real_1d_array lb_box;
	alglib::real_1d_array ub_box;
	lb_box.setcontent(bound_box.rows(), bound_box.col(0).data());
	ub_box.setcontent(bound_box.rows(), bound_box.col(1).data());
	alglib::minlpsetbc(Market.Problem, lb_box, ub_box);

	// -------------------------------------------------------------------------------
	// Set objective coefficients of variables
	// -------------------------------------------------------------------------------
	Eigen::VectorXd obj_vec = Eigen::VectorXd::Zero(variable_num);
	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		int row_start = 2 * Market.network.num_vertice + node_iter * (Market.price_intervals + 2);
		obj_vec.segment(row_start, Market.price_intervals + 2) = Market.bidded_price_map.bidded_price;
		obj_vec.segment(row_start, Market.price_intervals + 2) -= Eigen::VectorXd::Constant(Market.price_intervals + 2, Market.reference_price(node_iter));
	}
	alglib::real_1d_array obj_coeff;
	obj_coeff.setcontent(obj_vec.size(), obj_vec.data());
	alglib::minlpsetcost(Market.Problem, obj_coeff);

	// -------------------------------------------------------------------------------
	// Solve the problem
	// -------------------------------------------------------------------------------
	alglib::minlpoptimize(Market.Problem);
}
