// Source file for the functions of the power market
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

// ------------------------------------------------------------------------------------------------
// Functions involving multiple markets
// ------------------------------------------------------------------------------------------------
void power_market::Submitted_bid_calculation(int tick, market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, bool DSO_filter_flag){
	int point_num = Power_network_inform.points.bidding_zone.size();
	int sample_num = Power_market_inform.agent_profiles.end_users[0].size();

	// Initialize submit bids of markets
	Market_Initialization(Power_market_inform.International_Market);
	Market_Initialization(Power_market_inform.TSO_Market);
	for(int DSO_iter = 0; DSO_iter < Power_market_inform.DSO_Markets.size(); ++ DSO_iter){
		Market_Initialization(Power_market_inform.DSO_Markets[DSO_iter]);
		Power_market_inform.DSO_Markets[DSO_iter].filtered_supply = Power_market_inform.DSO_Markets[DSO_iter].submitted_supply;
		Power_market_inform.DSO_Markets[DSO_iter].filtered_demand = Power_market_inform.DSO_Markets[DSO_iter].submitted_demand;
	}
	Power_market_inform.industrial_submitted_demand = Eigen::MatrixXd::Zero(Power_market_inform.TSO_Market.price_intervals + 2, point_num);

	// Demand at each point (residential) / node (industrial)
	// Trivial case: residential demand at each point is 100% inflexible; 5% of industrial demand is flexible
	for(int point_iter = 0; point_iter < point_num; ++ point_iter){
		int point_ID = Power_network_inform.points.in_cluster_ID(point_iter);
		int node_ID = Power_network_inform.points.node(point_iter);
		int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
		int bz_ID = Power_network_inform.points.bidding_zone(point_iter);
		//std::cout << bz_ID << "\t" << node_ID << "\t" << Power_network_inform.nodes.bidding_zone(node_ID) << "\n";

		// Residential demand; connect to distribution power network
		for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
//			double bid_inflex_quan = Power_market_inform.end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile(0);
//			bid_inflex_quan *= Power_market_inform.end_user_profiles[point_iter][sample_iter].operation.weight;
//			bid_inflex_quan *= Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
//
//			double bid_flex_quan = Power_market_inform.end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile(0);
//			bid_flex_quan *= Power_market_inform.end_user_profiles[point_iter][sample_iter].operation.weight;
//			bid_flex_quan *= Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
//
//			// Updating inflexible bids; if quantity is positive (negative -> update supply bids)
//			int price_inflex_ID;
//			if(bid_inflex_quan >= 0.){
//				price_inflex_ID = Power_market_inform.end_user_profiles[point_iter][sample_iter].operation.demand_inflex_price_ID;
//				Power_market_inform.DSO_Markets[DSO_ID].submitted_demand(price_inflex_ID, point_ID) += bid_inflex_quan;
//				Power_market_inform.International_Market.submitted_demand(price_inflex_ID, bz_ID) += bid_inflex_quan;
//
//				// If DSOs do not filter local bids, the demand is added directly on the nodes of the TSO
//				if(!DSO_filter_flag){
//					Power_market_inform.TSO_Market.submitted_demand(price_inflex_ID, node_ID) += bid_inflex_quan;
//				}
//			}
//			else{
//				price_inflex_ID = Power_market_inform.end_user_profiles[point_iter][sample_iter].operation.supply_inflex_price_ID;
//				Power_market_inform.DSO_Markets[DSO_ID].submitted_supply(price_inflex_ID, point_ID) -= bid_inflex_quan;
//				Power_market_inform.International_Market.submitted_supply(price_inflex_ID, bz_ID) -= bid_inflex_quan;
//
//				// If DSOs do not filter local bids, the demand is added directly on the nodes of the TSO
//				if(!DSO_filter_flag){
//					Power_market_inform.TSO_Market.submitted_supply(price_inflex_ID, node_ID) -= bid_inflex_quan;
//				}
//			}
//
//			// Updating flexible bids; if quantity is positive (negative -> update supply bids)
//			int price_flex_ID;
//			if(bid_flex_quan >= 0.){
//				price_flex_ID = Power_market_inform.end_user_profiles[point_iter][sample_iter].operation.demand_flex_price_ID;
//				Power_market_inform.DSO_Markets[DSO_ID].submitted_demand(price_flex_ID, point_ID) += bid_flex_quan;
//				Power_market_inform.International_Market.submitted_demand(price_flex_ID, bz_ID) += bid_flex_quan;
//
//				// If DSOs do not filter local bids, the demand is added directly on the nodes of the TSO
//				if(!DSO_filter_flag){
//					Power_market_inform.TSO_Market.submitted_demand(price_flex_ID, node_ID) += bid_flex_quan;
//				}
//			}
//			else{
//				price_flex_ID = Power_market_inform.end_user_profiles[point_iter][sample_iter].operation.supply_flex_price_ID;
//				Power_market_inform.DSO_Markets[DSO_ID].submitted_supply(price_flex_ID, point_ID) -= bid_flex_quan;
//				Power_market_inform.International_Market.submitted_supply(price_flex_ID, bz_ID) -= bid_flex_quan;
//
//				// If DSOs do not filter local bids, the demand is added directly on the nodes of the TSO
//				if(!DSO_filter_flag){
//					Power_market_inform.TSO_Market.submitted_supply(price_flex_ID, node_ID) -= bid_flex_quan;
//				}
//			}
		}

		// Industrial demand; connect to transmission power network
		double bid_inflex_industrial = Power_network_inform.points.nominal_mean_demand_field(point_iter, tick);
		bid_inflex_industrial *= Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
		bid_inflex_industrial *= 1. - agent::parameters::residential_ratio();
		double bid_flex_industrial = bid_inflex_industrial;
		bid_inflex_industrial *= 1. - agent::industrial::flexible_ratio();
		bid_flex_industrial *= agent::industrial::flexible_ratio();
		Power_market_inform.industrial_submitted_demand(Power_market_inform.TSO_Market.price_intervals + 1, point_iter) += bid_inflex_industrial;
		Power_market_inform.industrial_submitted_demand.col(point_iter).segment(1, Power_market_inform.TSO_Market.price_intervals) += Eigen::VectorXd::Constant(Power_market_inform.TSO_Market.price_intervals, (double) bid_flex_industrial / Power_market_inform.TSO_Market.price_intervals);
		Power_market_inform.TSO_Market.submitted_demand(Power_market_inform.TSO_Market.price_intervals + 1, node_ID) += bid_inflex_industrial;
		Power_market_inform.TSO_Market.submitted_demand.col(node_ID).segment(1, Power_market_inform.TSO_Market.price_intervals) += Eigen::VectorXd::Constant(Power_market_inform.TSO_Market.price_intervals, (double) bid_flex_industrial / Power_market_inform.TSO_Market.price_intervals);
		Power_market_inform.International_Market.submitted_demand(Power_market_inform.TSO_Market.price_intervals + 1, bz_ID) += bid_inflex_industrial;
		Power_market_inform.International_Market.submitted_demand.col(bz_ID).segment(1, Power_market_inform.International_Market.price_intervals) += Eigen::VectorXd::Constant(Power_market_inform.International_Market.price_intervals, (double) bid_flex_industrial / Power_market_inform.International_Market.price_intervals);
	}

	// Supply at each point (LV power plants) / node (HV power plants)
	for(int hydro_iter = 0; hydro_iter < Power_network_inform.plants.hydro.node.size(); ++ hydro_iter){
		int bz_ID;
		int DSO_ID;
		int node_ID;
		Eigen::VectorXd bid_vec;

		// High voltage power plants connect directly to transmission network
		if(Power_network_inform.plants.hydro.cap(hydro_iter) >= 20.){
			node_ID = Power_network_inform.plants.hydro.node(hydro_iter);
			DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);
			bid_vec = Power_market_inform.International_Market.merit_order_curve.col(Power_network_inform.nodes.bidding_zone(node_ID));
			bid_vec *= Power_network_inform.plants.hydro.cap(hydro_iter);
			bid_vec /= (Power_market_inform.International_Market.merit_order_curve.col(Power_network_inform.nodes.bidding_zone(node_ID)).sum());

			Power_market_inform.TSO_Market.submitted_supply.col(node_ID) += bid_vec;
		}
		// Low voltage power plants feed into distribution network
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
			bid_vec = Power_market_inform.International_Market.merit_order_curve.col(Power_network_inform.points.bidding_zone(point_ID));
			bid_vec	*= Power_network_inform.plants.hydro.cap(hydro_iter);
			bid_vec	/= (Power_market_inform.International_Market.merit_order_curve.col(Power_network_inform.points.bidding_zone(point_ID)).sum());

			Power_market_inform.DSO_Markets[DSO_ID].submitted_supply.col(Power_network_inform.points.in_cluster_ID(point_ID)) += bid_vec;

			// If DSOs do not filter local bids, supply from LV buses is directly added to the nodes in the TSO
			if(!DSO_filter_flag){
				Power_market_inform.TSO_Market.submitted_supply.col(node_ID) += bid_vec;
			}
		}
		Power_market_inform.International_Market.submitted_supply.col(bz_ID) += bid_vec;
	}

	for(int wind_iter = 0; wind_iter < Power_network_inform.plants.wind.node.size(); ++ wind_iter){
		int bz_ID;
		int DSO_ID;
		int node_ID;
		int x_ID = int((Power_network_inform.plants.wind.x(wind_iter) - Power_network_inform.points.x.minCoeff()) / Power_network_inform.points.grid_length + .5);
		int y_ID = int((Power_network_inform.plants.wind.y(wind_iter) - Power_network_inform.points.y.minCoeff()) / Power_network_inform.points.grid_length + .5);
		int point_ID = Power_network_inform.points.coordinate_grid(x_ID, y_ID);
		if(point_ID == -1){
			continue;
		}
		double bid_quan = Power_network_inform.points.wind_on_cf(point_ID, tick) * Power_network_inform.plants.wind.cap(wind_iter);

		// High voltage power plants connect directly to transmission network
		if(Power_network_inform.plants.wind.cap(wind_iter) >= 20.){
			node_ID = Power_network_inform.plants.wind.node(wind_iter);
			DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);

			Power_market_inform.TSO_Market.submitted_supply(0, node_ID) += bid_quan;
		}
		// Low voltage power plants feed into distribution network
		else{
			node_ID = Power_network_inform.points.node(point_ID);
			DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);

			Power_market_inform.DSO_Markets[DSO_ID].submitted_supply(0, Power_network_inform.points.in_cluster_ID(point_ID)) += bid_quan;

			// If DSOs do not filter local bids, supply from LV buses is directly added to the nodes in the TSO
			if(!DSO_filter_flag){
				Power_market_inform.TSO_Market.submitted_supply(0, node_ID) += bid_quan;
			}
		}
		Power_market_inform.International_Market.submitted_supply(0, bz_ID) += bid_quan;
	}

	if(DSO_filter_flag){
		double demand_sum = 0;
		double supply_sum = 0;
		for(int DSO_iter = 0; DSO_iter < Power_market_inform.DSO_Markets.size(); ++ DSO_iter){
			demand_sum += Power_market_inform.DSO_Markets[DSO_iter].submitted_demand.sum();
			supply_sum +=  Power_market_inform.DSO_Markets[DSO_iter].submitted_supply.sum();
		}
		std::cout << "Total submitted:\t" << demand_sum << "\t" << supply_sum << "\n\n";
	}
}

void power_market::TSO_boundary_update(int tick, market_inform &TSO_Market, market_inform &International_Market, power_network::network_inform &Power_network_inform){
	// Store cross-border transmission flow as inflexible supply / demand at entry nodes of TSO
	for(int edge_iter = 0; edge_iter < International_Market.network.num_edges; ++ edge_iter){
		if(International_Market.network.incidence[edge_iter](1) < International_Market.cross_border_zone_start){
			continue;
		}
		bool sink_flag = International_Market.network.confirmed_power(tick, edge_iter) >= 0.;
		int price_ID = sink_flag * (International_Market.price_intervals + 1);
		double source = -(1 - sink_flag) * International_Market.network.confirmed_power(tick, edge_iter);
		double sink = sink_flag * International_Market.network.confirmed_power(tick, edge_iter);

		int zone_ID = International_Market.network.incidence[edge_iter](1) - International_Market.cross_border_zone_start;
		int node_num = Power_network_inform.cbt.entry_node_num[zone_ID];
		for(int node_iter = 0; node_iter < node_num; ++ node_iter){
			int node_ID = Power_network_inform.cbt.entry_nodes(zone_ID, node_iter);
			TSO_Market.submitted_supply(price_ID, node_ID) += source / node_num;
			TSO_Market.submitted_demand(price_ID, node_ID) += sink / node_num;
		}
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
}

void power_market::Flow_Based_Market_Optimization(market_inform &Market, alglib::minlpstate &Problem){
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
	alglib::minlpoptimize(Problem);

//	std::cout << "function end\n";
}

void power_market::Filtered_bid_calculation(int tick, markets_inform &DSO_Markets, market_inform &TSO_Market,  power_network::network_inform &Power_network_inform, std::vector <alglib::minlpstate> &DSO_Problems){
//	std::cout << TSO_Market.submitted_demand.sum() << "\t";
//	std::cout << TSO_Market.submitted_supply.sum() << "\n";

	// Correct formula
	for(int DSO_iter = 0; DSO_iter < DSO_Markets.size(); ++ DSO_iter){
		//std::cout << DSO_Markets[DSO_iter].submitted_demand.sum() << "\t";
		//std::cout << DSO_Markets[DSO_iter].submitted_supply.sum() << "\n";

		// Find merit order curve for filtered demand
		Eigen::MatrixXd submitted_supply = DSO_Markets[DSO_iter].submitted_supply;
		DSO_Markets[DSO_iter].submitted_supply = Eigen::MatrixXd::Zero(DSO_Markets[DSO_iter].price_intervals + 2, DSO_Markets[DSO_iter].num_zone);
		power_market::Source_Node_Set(DSO_Markets[DSO_iter], Power_network_inform.DSO_cluster[DSO_iter]);
		power_market::Flow_Based_Market_Optimization(DSO_Markets[DSO_iter], DSO_Problems[DSO_iter]);
		power_market::DSO_Market_Results_Get(tick, DSO_Markets[DSO_iter], DSO_Problems[DSO_iter], Power_network_inform.DSO_cluster[DSO_iter], 0);

		// Find merit order curve for filtered supply
		DSO_Markets[DSO_iter].submitted_supply = submitted_supply;
		DSO_Markets[DSO_iter].submitted_demand = Eigen::MatrixXd::Zero(DSO_Markets[DSO_iter].price_intervals + 2, DSO_Markets[DSO_iter].num_zone);
		power_market::Sink_Node_Set(DSO_Markets[DSO_iter], Power_network_inform.DSO_cluster[DSO_iter]);
		power_market::Flow_Based_Market_Optimization(DSO_Markets[DSO_iter], DSO_Problems[DSO_iter]);
		power_market::DSO_Market_Results_Get(tick, DSO_Markets[DSO_iter], DSO_Problems[DSO_iter], Power_network_inform.DSO_cluster[DSO_iter], 1);

		// Store the filtered results in the submitted bids of TSO
		for(int point_iter = 0; point_iter < Power_network_inform.DSO_cluster[DSO_iter].points_ID.size(); ++ point_iter){
			int point_ID = Power_network_inform.DSO_cluster[DSO_iter].points_ID[point_iter];
			int node_ID = Power_network_inform.points.node(point_ID);
			TSO_Market.submitted_supply.col(node_ID) += DSO_Markets[DSO_iter].filtered_supply.col(point_iter);
			TSO_Market.submitted_demand.col(node_ID) += DSO_Markets[DSO_iter].filtered_demand.col(point_iter);
		}
		//std::cout << "function 1 end\n";
	}
//	std::cout << TSO_Market.submitted_demand.sum() << "\t";
//	std::cout << TSO_Market.submitted_supply.sum() << "\n\n";

	double demand_sum = 0;
	double supply_sum = 0;
	for(int DSO_iter = 0; DSO_iter < DSO_Markets.size(); ++ DSO_iter){
		demand_sum += DSO_Markets[DSO_iter].filtered_demand.sum();
		supply_sum +=  DSO_Markets[DSO_iter].filtered_supply.sum();
	}
	std::cout << "Total filtered:\t" << demand_sum << "\t" << supply_sum << "\n\n";

//	# pragma omp parallel
//	{
//		# pragma omp for
//		for(int DSO_iter = 0; DSO_iter < 10; ++ DSO_iter){
//			auto Market = DSO_Markets[0];
//			auto Problem = DSO_Problems[0];
//			power_market::Flow_Based_Market_Optimization(Market, Problem);
//			std::cout << DSO_iter << "\n";
//		}
//	}
}
