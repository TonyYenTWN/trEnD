// Source file for the functions of the power market
//#include "../basic/alglib/optimization.h"
//#include "../power_network/power_network.h"
#include "power_market.h"

// ------------------------------------------------------------------------------------------------
// Generic functions for all market types
// ------------------------------------------------------------------------------------------------
void power_market::Market_Initialization(market_inform &Market){
	// Initialization of process variables
	// Should re-initialize for every time slice
	Market.submitted_supply = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Market.submitted_demand = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
}

void power_market::Market_clearing_nodal(int tick, market_inform &Market, Eigen::VectorXi &default_price_ID, Eigen::MatrixXd &bidded_supply, Eigen::MatrixXd &bidded_demand){
	Eigen::VectorXi price_demand_ID = (Market.price_intervals + 1) * Eigen::VectorXi::Ones(Market.num_zone);
	Eigen::VectorXi price_supply_ID = Eigen::VectorXi::Zero(Market.num_zone);

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
				double trade_quantity = std::min(bidded_supply(price_supply_ID(zone_ID), zone_ID), bidded_demand(price_demand_ID(zone_ID), zone_ID));
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
void power_market::Submitted_bid_calculation(int tick, markets_inform &DSO_Markets, market_inform &TSO_Market, market_inform &International_Market, power_network::network_inform &Power_network_inform, std::string fin_point_demand){
	// Calculation of submit bids at the beginning of each time slice
	auto fin_point_demand_dim = basic::get_file_dim(fin_point_demand);
	auto point_demand_inform = basic::read_file(fin_point_demand_dim[0], fin_point_demand_dim[1], fin_point_demand);

	// Initialize submit bids of markets
	Market_Initialization(TSO_Market);
	Market_Initialization(International_Market);
	for(int DSO_iter = 0; DSO_iter < DSO_Markets.size(); ++ DSO_iter){
		Market_Initialization(DSO_Markets[DSO_iter]);
	}

	// Declare variables for the loops
	Eigen::VectorXd bid_vec;

	// Trivial case: demand at each point are 100% inflexible
	for(int point_iter = 0; point_iter < Power_network_inform.points.bidding_zone.size(); ++ point_iter){
		int node_ID = Power_network_inform.points.node(point_iter);
		int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
		int bz_ID = Power_network_inform.points.bidding_zone(point_iter);
		double bid_quan = point_demand_inform(point_iter, 0) * Power_network_inform.points.population_density(point_iter); //* Power_network_inform.points.point_area;
		// nominal demand currently wrong in processed files, should change them and then multiply area of a point later

		DSO_Markets[DSO_ID].submitted_demand(DSO_Markets[DSO_ID].price_intervals + 1, Power_network_inform.points.in_cluster_ID(point_iter)) = bid_quan;
		TSO_Market.submitted_demand(DSO_Markets[DSO_ID].price_intervals + 1, node_ID) += bid_quan;
		International_Market.submitted_demand(DSO_Markets[DSO_ID].price_intervals + 1, bz_ID) += bid_quan;
	}

	// Supply at each point (LV power plants) / node (HV power plants)
	for(int hydro_iter = 0; hydro_iter < Power_network_inform.plants.hydro.node.size(); ++ hydro_iter){
		int bz_ID;
		int DSO_ID;
		int node_ID;

		if(Power_network_inform.plants.hydro.cap(hydro_iter) >= 20.){
			node_ID = Power_network_inform.plants.hydro.node(hydro_iter);
			DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);
			bid_vec = International_Market.merit_order_curve.col(Power_network_inform.nodes.bidding_zone(node_ID))
				* Power_network_inform.plants.hydro.cap(hydro_iter) / (International_Market.merit_order_curve.col(Power_network_inform.nodes.bidding_zone(node_ID)).sum());

			DSO_Markets[DSO_ID].submitted_supply.col(Power_network_inform.DSO_cluster[DSO_ID].points_ID.size() + Power_network_inform.nodes.in_cluster_ID(node_ID)) += bid_vec;
		}
		else{
			int x_ID = int((Power_network_inform.plants.hydro.x(hydro_iter) - Power_network_inform.points.x.minCoeff()) / Power_network_inform.points.grid_length + .5);
			int y_ID = int((Power_network_inform.plants.hydro.y(hydro_iter) - Power_network_inform.points.y.minCoeff()) / Power_network_inform.points.grid_length + .5);
			int point_ID = Power_network_inform.points.coordinate_grid(x_ID, y_ID);
			if(point_ID == -1){
				continue;
			}
			node_ID = Power_network_inform.points.node(point_ID);
			DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);
			bid_vec = International_Market.merit_order_curve.col(Power_network_inform.points.bidding_zone(point_ID))
				* Power_network_inform.plants.hydro.cap(hydro_iter) / (International_Market.merit_order_curve.col(Power_network_inform.points.bidding_zone(point_ID)).sum());

			DSO_Markets[DSO_ID].submitted_supply.col(Power_network_inform.points.in_cluster_ID(point_ID)) += bid_vec;
		}
		TSO_Market.submitted_supply.col(node_ID) += bid_vec;
		International_Market.submitted_supply.col(bz_ID) += bid_vec;
	}
}

// ------------------------------------------------------------------------------------------------
// Specific functions for for flow-based markets
// ------------------------------------------------------------------------------------------------
void power_market::Flow_Based_Market_LP_Set(market_inform &Market, alglib::minlpstate &Problem){
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
	double tol = 1E-16;
	//double tol = 0.;
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
	int constrant_num = 2 * Market.network.num_vertice + Market.network.num_edges + 1;
	int variable_num = Market.network.num_vertice * (Market.price_intervals + 4) + Market.network.num_edges;
	Eigen::VectorXpd non_zero_num(constrant_num);
	non_zero_num << Connection_num + Eigen::VectorXpd::Ones(Market.network.num_vertice), Eigen::VectorXpd::Constant(Market.network.num_edges, 3), Eigen::VectorXpd::Constant(Market.network.num_vertice, Market.price_intervals + 3), 1;
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

	// Voltage at reference node
	alglib::sparseset(constraint_general, non_zero_num.size() - 1, 0, 1.);

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

	// Bounds of box constraints
//	alglib::real_1d_array lb_box;
//	alglib::real_1d_array ub_box;
//	lb_box.setcontent(bound_box.rows(), bound_box.col(0).data());
//	ub_box.setcontent(bound_box.rows(), bound_box.col(1).data());

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
		obj_vec.segment(row_start, Market.price_intervals + 2) = Market.bidded_price;
	}
	alglib::real_1d_array obj_coeff;
	obj_coeff.setcontent(obj_vec.size(), obj_vec.data());

	// -------------------------------------------------------------------------------
	// Set the LP problem object
	// -------------------------------------------------------------------------------
	alglib::minlpcreate(variable_num, Problem);
	alglib::minlpsetcost(Problem, obj_coeff);
//	alglib::minlpsetbc(Problem, lb_box, ub_box);
	alglib::minlpsetlc2(Problem, constraint_general, lb_general, ub_general, constrant_num);
	alglib::minlpsetscale(Problem, scale);
	//alglib::minlpsetalgoipm(Problem);
	alglib::minlpsetalgodss(Problem, 0.);
	//Market.Problem = Problem;
}

void power_market::Flow_Based_Market_Optimization(int tick, market_inform &Market, alglib::minlpstate &Problem){
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
	alglib::minlpsetbc(Problem, lb_box, ub_box);

	// -------------------------------------------------------------------------------
	// Solve the problem
	// -------------------------------------------------------------------------------
	alglib::real_1d_array sol;
	alglib::minlpreport rep;
	alglib::minlpoptimize(Problem);
	alglib::minlpresults(Problem, sol, rep);

	// -------------------------------------------------------------------------------
	// Store the solution
	// -------------------------------------------------------------------------------
	Eigen::VectorXd sol_vec = Eigen::Map <Eigen::VectorXd> (sol.getcontent(), bound_box.rows());
	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		// Store power source / sink
		int row_start = 2 * Market.network.num_vertice + node_iter * (Market.price_intervals + 2);
		Market.confirmed_supply(tick, node_iter) = (sol_vec.segment(row_start, Market.price_intervals + 2).array().max(0)).sum();
		Market.confirmed_demand(tick, node_iter) = -(sol_vec.segment(row_start, Market.price_intervals + 2).array().min(0)).sum();
		Market.confirmed_price(tick, node_iter) = Market.bidded_price(0) + rep.lagbc[row_start];
		Market.confirmed_price(tick, node_iter) = std::min(Market.confirmed_price(tick, node_iter), Market.price_range_inflex(1));
		Market.confirmed_price(tick, node_iter) = std::max(Market.confirmed_price(tick, node_iter), Market.price_range_inflex(0));

		// Store voltage and power flow
		Market.network.confirmed_voltage.row(tick) = sol_vec.head(Market.network.num_vertice);
		Market.network.confirmed_power.row(tick) = sol_vec.tail(Market.network.num_edges);
	}
	//std::cout << Market.confirmed_price.row(tick) << "\n\n";

	//std::cout << sol_vec.segment(Market.network.num_vertice, Market.network.num_vertice).transpose() << "\n\n";
	std::cout << "Number of variables: " << variable_num << "\n\n";
	std::cout << sol_vec.segment(Market.network.num_vertice, Market.network.num_vertice).minCoeff() << " " << sol_vec.segment(Market.network.num_vertice, Market.network.num_vertice).maxCoeff() << " " << .5 * sol_vec.segment(Market.network.num_vertice, Market.network.num_vertice).array().abs().sum() << "\n";
	std::cout << sol_vec.head(Market.network.num_vertice).minCoeff() << " " << sol_vec.head(Market.network.num_vertice).maxCoeff()  << "\n";
	std::cout << sol_vec.tail(Market.network.num_edges).minCoeff() << " " << sol_vec.tail(Market.network.num_edges).maxCoeff() << "\n\n";
}
