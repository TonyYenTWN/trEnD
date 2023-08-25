// Source File for market clearing of the international market operator of the energy-only-market
#include <iostream>
//#include <chrono>
#include "src/basic/rw_csv.h"
#include "power_market.h"

namespace{
	void International_Market_Submitted_bid_calculation(int tick, power_market::market_inform &International_Market, power_network::network_inform &Power_network_inform){
		// Initialize submit bids of markets
		power_market::Market_Initialization(International_Market);

		// Demand at each point are 100% inflexible
		for(int point_iter = 0; point_iter < Power_network_inform.points.bidding_zone.size(); ++ point_iter){
			int node_ID = Power_network_inform.points.node(point_iter);
			//int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			int bz_ID = Power_network_inform.points.bidding_zone(point_iter);
			double bid_quan = Power_network_inform.points.nominal_mean_demand_field(point_iter, tick);
			bid_quan *= Power_network_inform.points.population_density(point_iter);
			bid_quan *= Power_network_inform.points.point_area / 1000.;

			International_Market.submitted_demand(International_Market.price_intervals + 1, bz_ID) += bid_quan;
		}

		// Supply at each power plant
		//double cutoff_power = agent::power_supplier::parameters::cutoff_power();
		for(int hydro_iter = 0; hydro_iter < Power_network_inform.plants.hydro.node.size(); ++ hydro_iter){
			int x_ID = int((Power_network_inform.plants.hydro.x(hydro_iter) - Power_network_inform.points.x.minCoeff()) / Power_network_inform.points.grid_length + .5);
			int y_ID = int((Power_network_inform.plants.hydro.y(hydro_iter) - Power_network_inform.points.y.minCoeff()) / Power_network_inform.points.grid_length + .5);
			int point_ID = Power_network_inform.points.coordinate_grid(x_ID, y_ID);
			if(point_ID == -1){
				continue;
			}
			int node_ID = Power_network_inform.points.node(point_ID);
			//int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);
			Eigen::VectorXd bid_vec = International_Market.merit_order_curve.col(bz_ID);
			bid_vec /= bid_vec.sum();
			bid_vec *= Power_network_inform.plants.hydro.cap(hydro_iter);
			International_Market.submitted_supply.col(bz_ID) += bid_vec;
		}

		for(int wind_iter = 0; wind_iter < Power_network_inform.plants.wind.node.size(); ++ wind_iter){
			int x_ID = int((Power_network_inform.plants.wind.x(wind_iter) - Power_network_inform.points.x.minCoeff()) / Power_network_inform.points.grid_length + .5);
			int y_ID = int((Power_network_inform.plants.wind.y(wind_iter) - Power_network_inform.points.y.minCoeff()) / Power_network_inform.points.grid_length + .5);
			int point_ID = Power_network_inform.points.coordinate_grid(x_ID, y_ID);
			if(point_ID == -1){
				continue;
			}
			int node_ID = Power_network_inform.points.node(point_ID);
			//int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);
			double bid_quan = Power_network_inform.points.wind_on_cf(point_ID, tick) * Power_network_inform.plants.wind.cap(wind_iter);

			International_Market.submitted_supply(0, bz_ID) += bid_quan;
		}
	}
}

void power_market::International_Market_Set(market_inform &International_Market, power_network::network_inform &Power_network_inform, int Time, fin_market fin_market){
	// Input Parameters of international market
	International_Market.num_zone = Power_network_inform.cbt.bz_names.size();
	International_Market.cross_border_zone_start = Power_network_inform.points.bidding_zone.maxCoeff() + 1;
	International_Market.time_intervals = Time;
	International_Market.zone_names = Power_network_inform.cbt.bz_names;
	parameters::bidded_price(International_Market.bidded_price_map);
	International_Market.network.num_vertice = International_Market.num_zone;
	International_Market.network.line_capacity_matrix = Power_network_inform.cbt.flow_constraint;

	// Construct incidence vector matrix
	International_Market.network.incidence.reserve(International_Market.network.num_vertice * International_Market.network.num_vertice);
	for(int row_iter = 0; row_iter < International_Market.network.num_vertice - 1; ++ row_iter){
		for(int col_iter = row_iter + 1; col_iter < International_Market.network.num_vertice; ++ col_iter){
			bool add_flag = 1 - (International_Market.network.line_capacity_matrix(row_iter, col_iter) < 0.) * (International_Market.network.line_capacity_matrix(col_iter, row_iter) < 0.);
			if(add_flag){
				International_Market.network.incidence.push_back(Eigen::Vector2i(row_iter, col_iter));
			}
		}
	}
	International_Market.network.line_capacity_matrix = Eigen::MatrixXd::Ones(International_Market.num_zone, International_Market.num_zone).array() * International_Market.network.line_capacity_matrix.array().max(0);

	// Construct power constraint matrix
	International_Market.network.num_edges = International_Market.network.incidence.size();
	International_Market.network.power_constraint = Eigen::MatrixXd(International_Market.network.num_edges, 2);
	for(int edge_iter = 0; edge_iter < International_Market.network.num_edges; ++ edge_iter){
		int row_ID = International_Market.network.incidence[edge_iter](0);
		int col_ID = International_Market.network.incidence[edge_iter](1);
		International_Market.network.power_constraint(edge_iter, 0) = International_Market.network.line_capacity_matrix(row_ID, col_ID);
		International_Market.network.power_constraint(edge_iter, 1) = International_Market.network.line_capacity_matrix(col_ID, row_ID);
	}

	// Quantity density at each price
	// Read inferred merit order curve data
	auto moc_dim = basic::get_file_dim(fin_market.moc);
	Eigen::MatrixXd merit_order_curve_q = basic::read_file(moc_dim[0], moc_dim[1], fin_market.moc);
	Eigen::MatrixXd diff_merit_order_curve_q = merit_order_curve_q.bottomRows(moc_dim[0] - 1) - merit_order_curve_q.topRows(moc_dim[0] - 1);
	International_Market.merit_order_curve = merit_order_curve_q;
	International_Market.merit_order_curve.bottomRows(moc_dim[0] - 1) = diff_merit_order_curve_q;
	International_Market.merit_order_curve = Eigen::MatrixXd::Ones(moc_dim[0], moc_dim[1]).array() * International_Market.merit_order_curve.array().max(0);

	// Default (residual) demand time series
	// Read default demand data
	auto demand_ts_dim = basic::get_file_dim(fin_market.demand);
	International_Market.demand_default = basic::read_file(demand_ts_dim[0], demand_ts_dim[1], fin_market.demand, 1);

	// Read cbt data
	auto cbt_ts_dim = basic::get_file_dim(fin_market.cbt);
	auto cbt_ts = basic::read_file(cbt_ts_dim[0], cbt_ts_dim[1], fin_market.cbt, 1);

	// Adjusted demand profile after cbt is taken out
	for(int edge_iter = 0; edge_iter < International_Market.network.incidence.size(); ++ edge_iter){
		int start_zone_ID = International_Market.network.incidence[edge_iter](0);
		int end_zone_ID = International_Market.network.incidence[edge_iter](1);
		int start_edge_ID = 2 * edge_iter;
		int end_edge_ID = 2 * edge_iter + 1;
		International_Market.demand_default.col(start_zone_ID) -= cbt_ts.col(start_edge_ID);
		International_Market.demand_default.col(start_zone_ID) += cbt_ts.col(end_edge_ID);
		International_Market.demand_default.col(end_zone_ID) += cbt_ts.col(start_edge_ID);
		International_Market.demand_default.col(end_zone_ID) -= cbt_ts.col(end_edge_ID);
	}

	// Calculate default residual demand profile from subtracting VRE for nations on the boundary
	auto solar_ts_dim = basic::get_file_dim(fin_market.solar);
	auto solar_ts = basic::read_file(solar_ts_dim[0], solar_ts_dim[1], fin_market.solar, 1);
	auto wind_on_ts_dim = basic::get_file_dim(fin_market.wind_on);
	auto wind_on_ts = basic::read_file(wind_on_ts_dim[0], wind_on_ts_dim[1], fin_market.wind_on, 1);
	auto wind_off_ts_dim = basic::get_file_dim(fin_market.wind_off);
	auto wind_off_ts = basic::read_file(wind_off_ts_dim[0], wind_off_ts_dim[1], fin_market.wind_off, 1);
	int boundary_num = International_Market.num_zone - International_Market.cross_border_zone_start;
	International_Market.demand_default.rightCols(boundary_num) -= solar_ts.rightCols(boundary_num) + wind_on_ts.rightCols(boundary_num) + wind_off_ts.rightCols(boundary_num);

	// Initialization of process variables
	power_market::Market_Initialization(International_Market);

	// Initialization of output variables
	International_Market.confirmed.supply = Eigen::MatrixXd::Zero(Time, International_Market.num_zone);
	International_Market.confirmed.demand = Eigen::MatrixXd::Zero(Time, International_Market.num_zone);
	International_Market.confirmed.price = Eigen::MatrixXd(Time, International_Market.num_zone);
	International_Market.confirmed.ratio_supply = Eigen::VectorXd::Zero(International_Market.num_zone);
	International_Market.confirmed.ratio_demand = Eigen::VectorXd::Zero(International_Market.num_zone);
	International_Market.network.confirmed_power = Eigen::MatrixXd::Zero(Time, International_Market.network.num_edges);
	International_Market.redispatch.price_demand = Eigen::MatrixXd::Zero(Time, International_Market.num_zone);
	International_Market.redispatch.price_supply = Eigen::MatrixXd::Zero(Time, International_Market.num_zone);
	Operation_Initialization(International_Market, Time);

	// Update alglib object for optimization
	// -------------------------------------------------------------------------------
	// Set matrix for general constraints
	// -------------------------------------------------------------------------------
	// Variables are sorted as {F}, {S}, {S}_0, {S}_1, {S}_2, ..., {I}
	// {F} are power flows across the edges
	// {S}_i is the set of source / sink at node #i at different price levels
	std::vector <Eigen::TripletXd> Constraint_trip;
	Constraint_trip.reserve(International_Market.network.num_vertice * (International_Market.price_intervals + 4)  + 2 * International_Market.network.num_edges);
	Eigen::VectorXpd Connection_num = Eigen::VectorXpd::Zero(International_Market.network.num_vertice);
	// Constraint for energy conservation
	for(int node_iter = 0; node_iter < International_Market.network.num_vertice; ++ node_iter){
		int col_ID =  International_Market.network.num_edges + node_iter;
		Constraint_trip.push_back(Eigen::TripletXd(node_iter, col_ID, 1.));
	}
	for(int edge_iter = 0; edge_iter < International_Market.network.num_edges; ++ edge_iter){
		int from_ID = International_Market.network.incidence[edge_iter](0);
		int to_ID = International_Market.network.incidence[edge_iter](1);
		Constraint_trip.push_back(Eigen::TripletXd(from_ID, edge_iter, -1.));
		Constraint_trip.push_back(Eigen::TripletXd(to_ID, edge_iter, 1.));
		Connection_num(from_ID) += 1;
		Connection_num(to_ID) += 1;
	}

	// Rows for source / sink summation at each node
	for(int node_iter = 0; node_iter < International_Market.network.num_vertice; ++ node_iter){
		int row_ID =  International_Market.network.num_vertice + node_iter;
		int col_ID =  International_Market.network.num_edges + node_iter;
		Constraint_trip.push_back(Eigen::TripletXd(row_ID, col_ID, 1.));

		int col_start =  International_Market.network.num_edges + International_Market.network.num_vertice + node_iter * (International_Market.price_intervals + 2);
		int col_end = col_start + International_Market.price_intervals + 1;
		for(int col_iter = col_start; col_iter <= col_end; ++ col_iter){
			Constraint_trip.push_back(Eigen::TripletXd(row_ID, col_iter, -1.));
		}
	}

	int constrant_num = 2 * International_Market.network.num_vertice;
	int variable_num =  International_Market.network.num_edges + International_Market.network.num_vertice * (International_Market.price_intervals + 3);
	Eigen::SparseMatrix <double, Eigen::RowMajor> Constraint(constrant_num, variable_num);
	Constraint.setFromTriplets(Constraint_trip.begin(), Constraint_trip.end());

	// Generate sparse matrix for general (equality) constraints nodal energy conservation and source / sink summation
	Eigen::VectorXpd non_zero_num(constrant_num);
	non_zero_num << Connection_num + Eigen::VectorXpd::Ones(International_Market.network.num_vertice), Eigen::VectorXpd::Constant(International_Market.network.num_vertice, International_Market.price_intervals + 3);
	alglib::integer_1d_array row_sizes_general;
	row_sizes_general.setcontent(non_zero_num.size(), non_zero_num.data());
	alglib::sparsematrix constraint_general;
	alglib::sparsecreatecrs(constrant_num, variable_num, row_sizes_general, constraint_general);
	for(int row_iter = 0; row_iter < Constraint.outerSize(); ++ row_iter){
		for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator inner_iter(Constraint, row_iter); inner_iter; ++ inner_iter){
			alglib::sparseset(constraint_general, inner_iter.row(), inner_iter.col(), inner_iter.value());
		}
	}

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
	// Set scale of variables
	// -------------------------------------------------------------------------------
	Eigen::VectorXd scale_vec = Eigen::VectorXd::Ones(variable_num);
	alglib::real_1d_array scale;
	scale.setcontent(scale_vec.size(), scale_vec.data());

	// -------------------------------------------------------------------------------
	// Set objective coefficients of variables
	// -------------------------------------------------------------------------------
	Eigen::VectorXd obj_vec = Eigen::VectorXd::Zero(variable_num);
	for(int node_iter = 0; node_iter < International_Market.network.num_vertice; ++ node_iter){
		int row_start = International_Market.network.num_edges + International_Market.network.num_vertice + node_iter * (International_Market.price_intervals + 2);
		obj_vec.segment(row_start, International_Market.price_intervals + 2) = International_Market.bidded_price_map.bidded_price;
	}
	alglib::real_1d_array obj_coeff;
	obj_coeff.setcontent(obj_vec.size(), obj_vec.data());

	// -------------------------------------------------------------------------------
	// Set the LP problem object
	// -------------------------------------------------------------------------------
	alglib::minlpcreate(variable_num, International_Market.Problem);
	alglib::minlpsetcost(International_Market.Problem, obj_coeff);
	alglib::minlpsetlc2(International_Market.Problem, constraint_general, lb_general, ub_general, constrant_num);
	alglib::minlpsetscale(International_Market.Problem, scale);
	alglib::minlpsetalgodss(International_Market.Problem, 0.);
}

void power_market::Submitted_bid_calculation(market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
	// Initialize submit bids of markets
	Market_Initialization(Power_market_inform.International_Market);

	// Demand at each point (residential) / node (industrial)
	int point_num = Power_network_inform.points.bidding_zone.size();
	int sample_num = Power_market_inform.agent_profiles.end_users[0].size();
	for(int point_iter = 0; point_iter < point_num; ++ point_iter){
		int bz_ID = Power_network_inform.points.bidding_zone(point_iter);

		// Residential demand; connect to distribution power network
		for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
			Power_market_inform.International_Market.submitted_supply.col(bz_ID) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_inflex;
			Power_market_inform.International_Market.submitted_supply.col(bz_ID)  += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_flex;
			Power_market_inform.International_Market.submitted_demand.col(bz_ID)  += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_inflex;
			Power_market_inform.International_Market.submitted_demand.col(bz_ID) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_flex;
		}

		// Industrial demand; connect to transmission power network
		Power_market_inform.International_Market.submitted_demand.col(bz_ID) += Power_market_inform.agent_profiles.industrial.HV[point_iter].bids.submitted_demand_flex;
	}

	// Supply at each point (LV power plants) / node (HV power plants)
	int hydro_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].point_ID;
		int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
		Power_market_inform.International_Market.submitted_supply.col(bz_ID) += Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].bids.submitted_supply_flex;
	}
	int hydro_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].point_ID;
		int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
		Power_market_inform.International_Market.submitted_supply.col(bz_ID) += Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids.submitted_supply_flex;
	}
	int wind_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].point_ID;
		int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
		Power_market_inform.International_Market.submitted_supply.col(bz_ID) += Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids.submitted_supply_flex;
	}
	int wind_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].point_ID;
		int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
		Power_market_inform.International_Market.submitted_supply.col(bz_ID) += Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.submitted_supply_flex;
	}
	int pump_HV_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV.size();
	for(int agent_iter = 0; agent_iter < pump_HV_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].point_ID;
		int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
		Power_market_inform.International_Market.submitted_supply.col(bz_ID) += Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].bids.submitted_supply_flex;
		Power_market_inform.International_Market.submitted_demand.col(bz_ID) += Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].bids.submitted_demand_flex;
	}
	int pump_LV_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
	for(int agent_iter = 0; agent_iter < pump_LV_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;
		int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
		Power_market_inform.International_Market.submitted_supply.col(bz_ID) += Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.submitted_supply_flex;
		Power_market_inform.International_Market.submitted_demand.col(bz_ID) += Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.submitted_demand_flex;
	}
	int slack_LV_num = Power_market_inform.agent_profiles.power_supplier.slack.LV_plant.size();
	for(int agent_iter = 0; agent_iter < slack_LV_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].point_ID;
		int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
		Power_market_inform.International_Market.submitted_supply.col(bz_ID) += Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].bids.submitted_supply_flex;
		Power_market_inform.International_Market.submitted_demand.col(bz_ID) += Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].bids.submitted_demand_flex;
	}
}

void power_market::International_Market_Optimization(int tick, market_inform &Market){
	// Update of submitted supply and demand bids at other nations
	// Bidding zones of other nations: assume inflexible supply or demand
	Market.submitted_supply.rightCols(Market.num_zone - Market.cross_border_zone_start) = Market.merit_order_curve.rightCols(Market.num_zone - Market.cross_border_zone_start);
	for(int zone_ID = Market.cross_border_zone_start; zone_ID < Market.num_zone; ++ zone_ID){
		if(Market.demand_default(tick, zone_ID) >= 0){
			Market.submitted_demand(Market.price_intervals, zone_ID) += Market.demand_default(tick, zone_ID);
		}
		else{
			Market.submitted_supply(0, zone_ID) += -Market.demand_default(tick, zone_ID);
		}
	}
//	std::cout << Market.submitted_supply.col(0).sum() << "\t";
//	std::cout << Market.submitted_supply.col(1).sum() << "\t";
//	std::cout << Market.submitted_supply.col(2).sum() << "\t";
//	std::cout << Market.submitted_supply.col(3).sum() << "\t";
//	std::cout << Market.submitted_supply.col(4).sum() << "\n";
//	std::cout << Market.submitted_demand.col(0).sum() << "\t";
//	std::cout << Market.submitted_demand.col(1).sum() << "\t";
//	std::cout << Market.submitted_demand.col(2).sum() << "\t";
//	std::cout << Market.submitted_demand.col(3).sum() << "\t";
//	std::cout << Market.submitted_demand.col(4).sum() << "\n\n";

	// -------------------------------------------------------------------------------
	// Update bounds for box constraints
	// -------------------------------------------------------------------------------
	int variable_num =  Market.network.num_edges + Market.network.num_vertice * (Market.price_intervals + 3);
	Eigen::MatrixXd bound_box = Eigen::MatrixXd::Zero(variable_num, 2);
	bound_box.col(0).head(Market.network.num_edges) = -Market.network.power_constraint.col(1);
	bound_box.col(1).head(Market.network.num_edges) = Market.network.power_constraint.col(0);
	bound_box.col(0).segment(Market.network.num_edges, Market.network.num_vertice) = Eigen::VectorXd::Constant(Market.network.num_vertice, -std::numeric_limits<double>::infinity());
	bound_box.col(1).segment(Market.network.num_edges, Market.network.num_vertice) = Eigen::VectorXd::Constant(Market.network.num_vertice, std::numeric_limits<double>::infinity());
	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		int row_start = Market.network.num_edges + Market.network.num_vertice + node_iter * (Market.price_intervals + 2);
		bound_box.middleRows(row_start, Market.price_intervals + 2).col(0) = -Market.submitted_demand.col(node_iter);
		bound_box.middleRows(row_start, Market.price_intervals + 2).col(1) = Market.submitted_supply.col(node_iter);
	}

	// Bounds of box constraints
	alglib::real_1d_array lb_box;
	alglib::real_1d_array ub_box;
	lb_box.setcontent(bound_box.rows(), bound_box.col(0).data());
	ub_box.setcontent(bound_box.rows(), bound_box.col(1).data());
	alglib::minlpsetbc(Market.Problem, lb_box, ub_box);

	// -------------------------------------------------------------------------------
	// Solve the problem and get results
	// -------------------------------------------------------------------------------
	alglib::minlpoptimize(Market.Problem);
	alglib::real_1d_array sol;
	alglib::minlpreport rep;
	alglib::minlpresults(Market.Problem, sol, rep);
	Eigen::VectorXd sol_vec = Eigen::Map <Eigen::VectorXd> (sol.getcontent(), sol.length());

	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		// Store power source / sink
		int row_start = Market.network.num_edges + Market.network.num_vertice + node_iter * (Market.price_intervals + 2);
		Market.confirmed.supply(tick, node_iter) = (sol_vec.segment(row_start, Market.price_intervals + 2).array().max(0)).sum();
		Market.confirmed.demand(tick, node_iter) = -(sol_vec.segment(row_start, Market.price_intervals + 2).array().min(0)).sum();

		// Store nodal prices
		Market.confirmed.price(tick, node_iter) = Market.bidded_price_map.bidded_price(0) + rep.lagbc[row_start];
		Market.confirmed.price(tick, node_iter) = int(Market.confirmed.price(tick, node_iter)) + .5;
		if(Market.confirmed.price(tick, node_iter) < Market.bidded_price_map.bidded_price(1)){
			Market.confirmed.price(tick, node_iter) = Market.bidded_price_map.bidded_price(0);
		}
		else if(Market.confirmed.price(tick, node_iter) > Market.bidded_price_map.bidded_price(Market.price_intervals)){
			Market.confirmed.price(tick, node_iter) = Market.bidded_price_map.bidded_price(Market.price_intervals + 1);
		}

		// Store ratio at nodes
		for(int price_iter = 0; price_iter < Market.price_intervals + 2; ++ price_iter){
			if(Market.bidded_price_map.bidded_price(price_iter) >= Market.confirmed.price(tick, node_iter) || price_iter == Market.price_intervals + 1){
				if(sol[row_start + price_iter] >= 0.){
					Market.confirmed.ratio_demand(node_iter) = std::min(Market.submitted_demand(price_iter, node_iter), Market.submitted_supply(price_iter, node_iter) - sol[row_start + price_iter]);
					Market.confirmed.ratio_supply(node_iter) = Market.confirmed.ratio_demand(node_iter) + sol[row_start + price_iter];
					Market.confirmed.ratio_demand(node_iter) /= Market.submitted_demand(price_iter, node_iter) + 1E-12;
					Market.confirmed.ratio_supply(node_iter) /= Market.submitted_supply(price_iter, node_iter) + 1E-12;
				}
				else{
					Market.confirmed.ratio_supply(node_iter) = std::min(Market.submitted_supply(price_iter, node_iter), Market.submitted_demand(price_iter, node_iter) + sol[row_start + price_iter]);
					Market.confirmed.ratio_demand(node_iter) = Market.confirmed.ratio_supply(node_iter) - sol[row_start + price_iter];
					Market.confirmed.ratio_demand(node_iter) /= Market.submitted_demand(price_iter, node_iter) + 1E-12;
					Market.confirmed.ratio_supply(node_iter) /= Market.submitted_supply(price_iter, node_iter) + 1E-12;
				}
				break;
			}
		}
	}
//	std::cout << Market.confirmed.supply.row(tick) << "\n";
//	std::cout << Market.confirmed.demand.row(tick) << "\n\n";

	// Store cross-border transmission flow
	Market.network.confirmed_power.row(tick) = sol_vec.head(Market.network.num_edges);
}

void power_market::International_Market_Price_Estimation(int tick, market_inform &International_Market, power_network::network_inform &Power_network_inform, configuration::process_config &process_par){
	int foresight_time = agent::aggregator::parameters::foresight_time();

	// Initialization of forecast market clearing price
	if(tick == process_par.time_boundary[0]){
		for(int tock = 0; tock < foresight_time; ++ tock){
			International_Market_Submitted_bid_calculation(tick + tock, International_Market, Power_network_inform);
			International_Market_Optimization(tick + tock, International_Market);
		}
	}
	// Find the profile one time step further
	else{
		International_Market_Submitted_bid_calculation(tick + foresight_time - 1, International_Market, Power_network_inform);
		International_Market_Optimization(tick + foresight_time - 1, International_Market);
	}
}

//int main(){
//	// Input variables
//	int Time = 8760;
//	std::string fin_name_moc = "input/merit_order_curve_q_assimilated_2021.csv";
//	std::string fin_name_demand = "input/residual_load_default_forecast_2021.csv";
//	market_inform International_Market;
//	International_Market_Set(International_Market, Time, fin_name_moc, fin_name_demand);
//
//	// Naive market clearing
//	for(int tick = 0; tick < 10; ++ tick){
//		Market_Initialization(International_Market);
//		International_Market_Optimization(tick, International_Market, 1);
//	}
//
//	// Write csv file
//	// International_Market_Output(International_Market);
//}
