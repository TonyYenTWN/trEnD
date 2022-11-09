// Source file for power flow analysis
#include "power_flow_analysis.h"

// HELM method
// P-Q buses
// Set {V}(s) = \sum {a}_n * s^n
// {\hat V}(s) = \sum {b}_n * s^n
// {1. / V}(s) = \sum {c}_n * s^n
// {1. / \hat V}(s) = \sum {d}_n * s^n
// Equations:
// [Y] {V}(s) = s * {Conj(S)} . {1. / \hat V}(s)
// [Conj(Y)] {\hat V}(s) = s * {S} . {1. / V}(s)
//
// P-U buses
// Set {a}, {b}, {c}, {d} the same as P-Q bus
// {Q}(s) = \sum {e}_n * s^n
// Equations:
// [Y] {V}(s) = s * ({P} - j {Q}(s)) . {1. / \hat V}(s)
// [Conj(Y)] {\hat V}(s) = s * ({P} + j {Q}(s)) . {1. / V}(s)
// {V}(s) * {\hat V}(s) = {U_0^2} + s * ({U^2} - {U_0^2})
//
// Conservation law for currents
// {S} . {1. / V}(s)  = 0.
// {Conj(S)} . {1. / \hat V}(s) = 0.

void power_network::HELM_Set(network_inform &Power_network_inform, power_market::market_whole_inform &Power_market_inform){
	int node_num = Power_network_inform.nodes.bidding_zone.size();
	int point_num = Power_network_inform.points.bidding_zone.size();
	int edge_trans_num = Power_network_inform.edges.distance.size();
	int edge_distr_num = 0;
	int DSO_num = Power_network_inform.DSO_cluster.size();
	for(int DSO_iter = 0; DSO_iter < DSO_num; ++ DSO_iter){
		edge_distr_num += Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() * (Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() - 1) / 2;
	}
	double pi = boost::math::constants::pi<double>();

	// -------------------------------------------------------------------------------
	// Set the nodal admittance matrix
	// -------------------------------------------------------------------------------
	std::vector<Eigen::TripletXcd> Y_n_trip;
	Y_n_trip.reserve(node_num + point_num + 2 * (edge_trans_num + edge_distr_num + node_num));
	Eigen::VectorXcd Y_n_Diag = Eigen::VectorXcd::Zero(node_num + point_num);
//	Eigen::MatrixXcs Y_n_dense = Eigen::MatrixXcd::Zero(node_num + point_num, node_num + point_num);

	// Transmission level
	for(int edge_iter = 0; edge_iter < edge_trans_num; ++ edge_iter){
		int from_ID = Power_network_inform.edges.from(edge_iter);
		int to_ID = Power_network_inform.edges.to(edge_iter);
		int voltage = Power_network_inform.edges.voltage_base(edge_iter);

		// Series admittance
		std::complex <double> y_series(1., 0.);
		y_series /= Power_network_inform.edges.distance(edge_iter);
		y_series /= Power_network_inform.tech_parameters.z_trans_series;
		y_series *= Power_network_inform.tech_parameters.impedenace_base_levels[voltage];

		// Shunt admittance
		std::complex <double> y_shunt(1., 0.);
		y_shunt *= Power_network_inform.edges.distance(edge_iter);
		y_shunt *= Power_network_inform.tech_parameters.y_trans_shunt;
		y_shunt *= Power_network_inform.tech_parameters.impedenace_base_levels[voltage];

//		// Update nodal admittance matrix
//		Y_n_dense(from_ID, to_ID) = -y_series;
//		Y_n_dense(to_ID, from_ID) = -y_series;
//		Y_n_dense(from_ID, from_ID) += y_series + .5 * y_shunt;
//		Y_n_dense(to_ID, to_ID) += y_series + .5 * y_shunt;

		// Triplet for series impedence
		Y_n_trip.push_back(Eigen::TripletXcd(from_ID, to_ID, -y_series));
		Y_n_trip.push_back(Eigen::TripletXcd(to_ID, from_ID, -y_series));

		// Update diagonal terms
		Y_n_Diag(from_ID) += y_series + .5 * y_shunt;
		Y_n_Diag(to_ID) += y_series + .5 * y_shunt;
	}

	// Distribution level
	double z_base_low = Power_network_inform.tech_parameters.z_base_distr();
	double z_base_high = Power_network_inform.tech_parameters.z_base_conn();
	for(int DSO_iter = 0; DSO_iter < DSO_num; ++ DSO_iter){
		int DSO_point_num = Power_network_inform.DSO_cluster[DSO_iter].points_ID.size();
		int DSO_node_num = Power_network_inform.DSO_cluster[DSO_iter].nodes_ID.size();

		double partition_func = 0.;
		Eigen::MatrixXd num_line = Eigen::MatrixXd::Ones(DSO_point_num, DSO_point_num);
		Eigen::MatrixXd distance = Eigen::MatrixXd::Zero(DSO_point_num, DSO_point_num);

		// Connection between points
		for(int row_iter = 0; row_iter < DSO_point_num - 1; ++ row_iter){
			for(int col_iter = row_iter + 1; col_iter < DSO_point_num ; ++ col_iter){
				int point_ID_1 = Power_network_inform.DSO_cluster[DSO_iter].points_ID[row_iter];
				int point_ID_2 = Power_network_inform.DSO_cluster[DSO_iter].points_ID[col_iter];

				distance(row_iter, col_iter) = Power_network_inform.points.distance(point_ID_1, point_ID_2);
				num_line(row_iter, col_iter) = 1. / pow(distance(row_iter, col_iter) * 1E-3, 1. + Power_network_inform.tech_parameters.fraction_dim_distr);
				partition_func += num_line(row_iter, col_iter);
			}
		}
		num_line /= partition_func;
		num_line *= Power_network_inform.tech_parameters.line_density_distr * DSO_point_num;

		for(int row_iter = 0; row_iter < DSO_point_num - 1; ++ row_iter){
			for(int col_iter = row_iter + 1; col_iter < DSO_point_num ; ++ col_iter){
				int point_ID_1 = Power_network_inform.DSO_cluster[DSO_iter].points_ID[row_iter];
				int point_ID_2 = Power_network_inform.DSO_cluster[DSO_iter].points_ID[col_iter];

				// Series admittance
				std::complex <double> y_series(1., 0.);
				y_series /= distance(row_iter, col_iter);
				y_series /= Power_network_inform.tech_parameters.z_distr_series;
				y_series *= z_base_low;
				y_series *= num_line(row_iter, col_iter);

				// Shunt admittance
				std::complex <double> y_shunt(1., 0.);
				y_shunt *= distance(row_iter, col_iter);
				y_shunt *= Power_network_inform.tech_parameters.y_distr_shunt;
				y_shunt *= z_base_low;
				y_shunt *= num_line(row_iter, col_iter);

//				// Update nodal admittance matrix
//				Y_n_dense(node_num + point_ID_1, node_num + point_ID_2) = -y_series;
//				Y_n_dense(node_num + point_ID_2, node_num + point_ID_1) = -y_series;
//				Y_n_dense(node_num + point_ID_1, node_num + point_ID_1) += y_series + .5 * y_shunt;
//				Y_n_dense(node_num + point_ID_2, node_num + point_ID_2) += y_series + .5 * y_shunt;

				// Triplet for series impedence
				Y_n_trip.push_back(Eigen::TripletXcd(node_num + point_ID_1, node_num + point_ID_2, -y_series));
				Y_n_trip.push_back(Eigen::TripletXcd(node_num + point_ID_2, node_num + point_ID_1, -y_series));

				// Update diagonal terms
				Y_n_Diag(node_num + point_ID_1) += y_series + .5 * y_shunt;
				Y_n_Diag(node_num + point_ID_2) += y_series + .5 * y_shunt;
			}
		}

		// Connection between nodes and points
		for(int row_iter = 0; row_iter < DSO_node_num; ++ row_iter){
			int node_ID = Power_network_inform.DSO_cluster[DSO_iter].nodes_ID[row_iter];
			int min_point_ID;
			double distance_min = std::numeric_limits<double>::infinity();

			for(int col_iter = 0; col_iter < DSO_point_num ; ++ col_iter){
				int point_ID = Power_network_inform.DSO_cluster[DSO_iter].points_ID[col_iter];
				Eigen::Vector2d point_coor = Eigen::Vector2d(Power_network_inform.points.lon(point_ID), Power_network_inform.points.lat(point_ID));
				Eigen::Vector2d node_coor = Eigen::Vector2d(Power_network_inform.nodes.lon(node_ID), Power_network_inform.nodes.lat(node_ID));
				point_coor *= pi / 180.;
				node_coor *= pi / 180.;
				double distance_temp = spatial_field::geodist(point_coor, node_coor);

				if(distance_temp < distance_min){
					distance_min = distance_temp;
					min_point_ID = point_ID;
				}
			}

			// Series admittance
			std::complex <double> y_series(1., 0.);
			y_series /= distance_min;
			y_series /= Power_network_inform.tech_parameters.z_conn_series;
			y_series *= z_base_high;
			y_series *= Power_network_inform.tech_parameters.line_density_conn * DSO_point_num / DSO_node_num;

			// Shunt admittance
			std::complex <double> y_shunt(1., 0.);
			y_shunt *= distance_min;
			y_shunt *= Power_network_inform.tech_parameters.y_conn_shunt;
			y_shunt *= z_base_high;
			y_shunt *= Power_network_inform.tech_parameters.line_density_conn * DSO_point_num / DSO_node_num;

//			// Update nodal admittance matrix
//			Y_n_dense(node_ID, node_num + min_point_ID) = -y_series;
//			Y_n_dense(node_num + min_point_ID, node_ID) = -y_series;
//			Y_n_dense(node_ID, node_ID) += y_series + .5 * y_shunt;
//			Y_n_dense(node_num + min_point_ID, node_num + min_point_ID) += y_series + .5 * y_shunt;

			// Triplet for series impedence
			Y_n_trip.push_back(Eigen::TripletXcd(node_ID, node_num + min_point_ID, -y_series));
			Y_n_trip.push_back(Eigen::TripletXcd(node_num + min_point_ID, node_ID, -y_series));

			// Update diagonal terms
			Y_n_Diag(node_num + min_point_ID) += y_series + .5 * y_shunt;
			Y_n_Diag(node_ID) += y_series + .5 * y_shunt;
		}
	}

	// Triplet for diagonal terms
	for(int node_iter = 0; node_iter < node_num + point_num; ++ node_iter){
		Y_n_trip.push_back(Eigen::TripletXcd(node_iter, node_iter, Y_n_Diag(node_iter)));
	}

	// Store the nodal admittance matrix and
	Power_network_inform.power_flow.nodal_admittance = Eigen::SparseMatrix <std::complex <double>> (node_num + point_num, node_num + point_num);
	Power_network_inform.power_flow.nodal_admittance.setFromTriplets(Y_n_trip.begin(), Y_n_trip.end());

	// -------------------------------------------------------------------------------
	// Determine type of buses
	// -------------------------------------------------------------------------------
	// Assume only HV power suppliers can give reactive power
	Power_network_inform.power_flow.PQ_bus.reserve(node_num);
	Power_network_inform.power_flow.PU_bus.reserve(node_num);
	Power_network_inform.power_flow.ref_bus.reserve(node_num);
	Eigen::VectorXi node_type = Eigen::VectorXi::Zero(node_num);
	int node_ref_num;
	for(int edge_iter = 0; edge_iter < Power_network_inform.cbt.entry_bz.size(); ++ edge_iter){
		if(Power_network_inform.cbt.entry_node_num(edge_iter) == 0){
			continue;
		}
		node_ref_num = (int) Power_network_inform.cbt.entry_nodes(edge_iter, 0);
		break;
	}
	std::cout << node_ref_num << "\n";
	node_type(node_ref_num) = 2;

	int hydro_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		if(node_type(node_ID) != 0){
			continue;
		}

		node_type(node_ID) = 1;
	}

	int wind_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		if(node_type(node_ID) != 0){
			continue;
		}

		node_type(node_ID) = 1;
	}

	int pump_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV.size();
	for(int agent_iter = 0; agent_iter < pump_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		if(node_type(node_ID) != 0){
			continue;
		}

		node_type(node_ID) = 1;
	}

	for(int node_iter = 1; node_iter < node_num; ++ node_iter){
		switch (node_type(node_iter)){
			case 0:
				Power_network_inform.power_flow.PQ_bus.push_back(node_iter);
			case 1:
				Power_network_inform.power_flow.PU_bus.push_back(node_iter);
			case 2:
				Power_network_inform.power_flow.ref_bus.push_back(node_iter);
		}
	}

	// -------------------------------------------------------------------------------
	// Set the solver for iterative process
	// -------------------------------------------------------------------------------
	// Entries from original nodal admittance matrix, reordered and delete reference bus
	std::vector<Eigen::TripletXcd> Permut_Y_n_trip;
	Permut_Y_n_trip.reserve(node_num + point_num);
	for(int node_iter = 0; node_iter < Power_network_inform.power_flow.PU_bus.size(); ++ node_iter){
		int row_ID = node_iter;
		int bus_ID = Power_network_inform.power_flow.PU_bus[node_iter];
		Permut_Y_n_trip.push_back(Eigen::TripletXcd(row_ID, bus_ID, 1.));
	}
	for(int node_iter = 0; node_iter < Power_network_inform.power_flow.PQ_bus.size(); ++ node_iter){
		int row_ID = Power_network_inform.power_flow.PU_bus.size() + node_iter;
		int bus_ID = Power_network_inform.power_flow.PQ_bus[node_iter];
		Permut_Y_n_trip.push_back(Eigen::TripletXcd(row_ID, bus_ID, 1.));
	}
	for(int point_iter = 0; point_iter < point_num; ++ point_iter){
		int row_ID = Power_network_inform.power_flow.PU_bus.size() + Power_network_inform.power_flow.PQ_bus.size() + point_iter;
		int bus_ID = node_num + point_iter;
		Permut_Y_n_trip.push_back(Eigen::TripletXcd(row_ID, bus_ID, 1.));
	}
	Eigen::SparseMatrix <std::complex <double>> Permut_Y_n(Permut_Y_n_trip.size(), node_num + point_num);
	Permut_Y_n.setFromTriplets(Permut_Y_n_trip.begin(), Permut_Y_n_trip.end());
	Eigen::SparseMatrix <std::complex <double>> Y_n_small = Permut_Y_n * Power_network_inform.power_flow.nodal_admittance * Permut_Y_n.transpose();

	// Matrix for solver
	std::vector<Eigen::TripletXcd> Mat_trip;
	Mat_trip.reserve(2 * Y_n_small.nonZeros() + 3 * Power_network_inform.power_flow.PU_bus.size() + 2 * (node_num + point_num));

	// Entries from reduced nodal admittance matrix
	for(int col_iter = 0; col_iter < Y_n_small.outerSize(); ++ col_iter){
		for(Eigen::SparseMatrix<std::complex <double>>::InnerIterator inner_iter(Y_n_small, col_iter); inner_iter; ++ inner_iter){
			std::complex <double> y_conj = inner_iter.value();
			y_conj = std::conj(y_conj);

			Mat_trip.push_back(Eigen::TripletXcd(inner_iter.row() - 1, inner_iter.col(), inner_iter.value()));
			Mat_trip.push_back(Eigen::TripletXcd(2 * (node_num + point_num) + inner_iter.row() - 1, 2 * (node_num + point_num) + inner_iter.col(), y_conj));
		}
	}

//	// Entries from original nodal admittance matrix, reordered
//	for(int row_iter = 0; row_iter < Power_network_inform.power_flow.PU_bus.size(); ++ row_iter){
//		int row_ID = Power_network_inform.power_flow.PU_bus[row_iter];
//
//		// Columns of PU-bus
//		for(int col_iter = 0; col_iter < Power_network_inform.power_flow.PU_bus.size(); ++ col_iter){
//			int col_ID = Power_network_inform.power_flow.PU_bus[col_iter];
//			Mat_trip_const(Eigen::TripletXcd(row_iter, col_iter, Y_n_dense(row_ID, col_ID)));
//		}
//	}


//	// Find reduced nodal admittance matrix
//	Eigen::SparseMatrix <std::complex <double>> Y_n_small = Power_network_inform.power_flow.nodal_admittance.bottomRows(node_num + point_num - 1);
//	Y_n_small = Y_n_small.rightCols(node_num + point_num - 1);
//	Power_network_inform.power_flow.solver_reg.compute(Y_n_small);
//	Power_network_inform.power_flow.solver_hat.compute(Y_n_small.conjugate());

	// -------------------------------------------------------------------------------
	// Initialize
	// -------------------------------------------------------------------------------
	int Time = configuration::parameters::Time();

	Power_network_inform.power_flow.P_edge = Eigen::MatrixXcd::Zero(Time, edge_trans_num);
	Power_network_inform.power_flow.Q_edge = Eigen::MatrixXcd::Zero(Time, edge_trans_num);
	Power_network_inform.power_flow.P_node = Eigen::MatrixXcd::Zero(Time, node_num + point_num);
	Power_network_inform.power_flow.Q_node = Eigen::MatrixXcd::Zero(Time, node_num + point_num);
	Power_network_inform.power_flow.voltage_abs = Eigen::MatrixXcd::Zero(Time, node_num + point_num);
	Power_network_inform.power_flow.voltage_arg = Eigen::MatrixXcd::Zero(Time, node_num + point_num);
}

void power_network::HELM_Node_Update(int tick, network_inform &Power_network_inform, power_market::market_whole_inform &Power_market_inform){
	int node_num = Power_network_inform.nodes.bidding_zone.size();
	int point_num = Power_network_inform.points.bidding_zone.size();

	// Update from cross-border flows
	int edge_num = Power_market_inform.agent_profiles.cross_border.size();
	for(int edge_iter = 0; edge_iter < edge_num; ++ edge_iter){
		int node_num = Power_market_inform.agent_profiles.cross_border[edge_iter].node_num;
		int entry_bz_ID = Power_market_inform.agent_profiles.cross_border[edge_iter].entry_bz_ID;

		if(node_num == 0){
			continue;
		}
		for(int node_iter = 0; node_iter < node_num; ++ node_iter){
			int node_ID = Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].node_ID;
			double real_power;
			real_power = Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results.actual_supply;
			real_power -= Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results.actual_demand;
			Power_network_inform.power_flow.P_node(tick, node_ID) += real_power;
		}
	}

	// Update from end-users
	{
		int sample_num = agent::end_user::parameters::sample_num();
		double apparent_power_base = pow(1. - pow(agent::end_user::parameters::power_factor(), 2.), .5);
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				double real_power;
				real_power = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_supply;
				Power_network_inform.power_flow.P_node(tick, node_num + point_iter) += real_power;
				real_power = -Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_demand;
				Power_network_inform.power_flow.P_node(tick, node_num + point_iter) += real_power;
				Power_network_inform.power_flow.Q_node(tick, node_num + point_iter) += real_power * apparent_power_base;
			}
		}
	}

	// Update from industrial demand
	{
		int industrial_HV_num = Power_market_inform.agent_profiles.industrial.HV.size();
		double apparent_power_base = pow(1. - pow(agent::end_user::parameters::power_factor(), 2.), .5);
		for(int agent_iter = 0; agent_iter < industrial_HV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.industrial.HV[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);

			double real_power = -Power_market_inform.agent_profiles.industrial.HV[agent_iter].results.actual_demand;
			Power_network_inform.power_flow.P_node(tick, node_num + point_ID) += real_power;
			Power_network_inform.power_flow.Q_node(tick, node_ID) += real_power * apparent_power_base;
		}
	}

	// Update from power suppliers
	int hydro_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		double real_power = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		Power_network_inform.power_flow.P_node(tick, node_ID) += real_power;
	}

	int hydro_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].point_ID;

		double real_power = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		Power_network_inform.power_flow.P_node(tick, node_num + point_ID) += real_power;
	}

	int wind_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		double real_power = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		Power_network_inform.power_flow.P_node(tick, node_ID) += real_power;
	}

	int wind_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].point_ID;

		double real_power = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		Power_network_inform.power_flow.P_node(tick, node_num + point_ID) += real_power;
	}

	int pump_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV.size();
	for(int agent_iter = 0; agent_iter < pump_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		double real_power = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		real_power -= Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].results.actual_demand;
		Power_network_inform.power_flow.P_node(tick, node_ID) += real_power;
	}

	int pump_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
	for(int agent_iter = 0; agent_iter < pump_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;

		double real_power = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		real_power -= Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].results.actual_demand;
		Power_network_inform.power_flow.P_node(tick, node_num + point_ID) += real_power;
	}
}

void power_network::HELM_Solve(int tick, network_inform &Power_network_inform){
//	int node_num = Power_network_inform.nodes.bidding_zone.size();
//	int point_num = Power_network_inform.points.bidding_zone.size();
//
//	// -------------------------------------------------------------------------------
//	// Initialization of power series coefficients
//	// -------------------------------------------------------------------------------
//	int power_terms = 100;
//	Eigen::MatrixXcd V_up_reg = Eigen::MatrixXcd::Zero(node_num + point_num, power_terms);
//	Eigen::MatrixXcd V_up_hat = Eigen::MatrixXcd::Zero(node_num + point_num, power_terms);
//	Eigen::MatrixXcd V_down_reg = Eigen::MatrixXcd::Zero(node_num + point_num, power_terms);
//	Eigen::MatrixXcd V_down_hat = Eigen::MatrixXcd::Zero(node_num + point_num, power_terms);
//
//	// -------------------------------------------------------------------------------
//	// Trivial solution (V = 1 everywhere)
//	// -------------------------------------------------------------------------------
//	V_up_reg.col(0) = Eigen::VectorXcd::Ones(node_num + point_num);
//	V_up_hat.col(0) = Eigen::VectorXcd::Ones(node_num + point_num);
//	V_down_reg.col(0) = Eigen::VectorXcd::Ones(node_num + point_num);
//	V_down_hat.col(0) = Eigen::VectorXcd::Ones(node_num + point_num);
//
//	// -------------------------------------------------------------------------------
//	// Iteratively solve the linear equations
//	// -------------------------------------------------------------------------------
//	for(int term_iter = 1; term_iter < power_terms; ++ term_iter){
//		// Solve power flow relations for V_reg and 1 / V_hat
//		{
//			Eigen::VectorXcd rhs = Power_network_inform.power_flow.power_node.row(tick).tail(node_num + point_num - 1).conjugate().transpose();
//			rhs = rhs.array() * V_down_hat.col(term_iter - 1).tail(node_num + point_num - 1).array();
//			V_up_reg.col(term_iter).tail(node_num + point_num - 1) = Power_network_inform.power_flow.solver_reg.solve(rhs);
//		}
//
//		// Solve power flow relations for V_hat and 1 / V_reg
//		{
//			Eigen::VectorXcd rhs = Power_network_inform.power_flow.power_node.row(tick).tail(node_num + point_num - 1).transpose();
//			rhs = rhs.array() * V_down_reg.col(term_iter - 1).tail(node_num + point_num - 1).array();
//			V_up_hat.col(term_iter).tail(node_num + point_num - 1) = Power_network_inform.power_flow.solver_reg.solve(rhs);
//		}
//
//		// Update reciprocal relation for V_reg and 1 / V_reg
//		{
//			Eigen::VectorXcd remnant = V_up_reg.col(term_iter).array() * V_down_reg.col(0).array();
//			for(int term_iter_2 = 1; term_iter_2 < term_iter; ++ term_iter_2){
//				remnant = remnant.array() + V_up_reg.col(term_iter - term_iter_2).array() * V_down_reg.col(term_iter_2).array();
//			}
//			V_down_reg.col(term_iter) = -remnant.array() / V_up_reg.col(0).array();
//		}
//
//		// Update reciprocal relation for V_hat and 1 / V_hat
//		{
//			Eigen::VectorXcd remnant = V_up_hat.col(term_iter).array() * V_down_hat.col(0).array();
//			for(int term_iter_2 = 1; term_iter_2 < term_iter; ++ term_iter_2){
//				remnant = remnant.array() + V_up_hat.col(term_iter - term_iter_2).array() * V_down_hat.col(term_iter_2).array();
//			}
//			V_down_hat.col(term_iter) = -remnant.array() / V_up_hat.col(0).array();
//		}
//	}
//
//	// Sanity check
//	Eigen::VectorXcd V_reg_dir = V_up_reg * Eigen::VectorXcd::Ones(node_num + point_num);
//	std::cout << "AC Power Flow:\n";
//	std::cout << "Transmission Level:\n";
//	std::cout << "Voltage:\n" << V_reg_dir.head(node_num).array().abs().minCoeff() << "\t" << V_reg_dir.head(node_num).array().abs().maxCoeff() << "\n";
//	std::cout << "Phase:\n" << V_reg_dir.head(node_num).array().arg().minCoeff() << "\t" << V_reg_dir.head(node_num).array().arg().maxCoeff() << "\n\n";
//	std::cout << "Distribution Level:\n";
//	std::cout << "Voltage:\n" << V_reg_dir.tail(point_num).array().abs().minCoeff() << "\t" << V_reg_dir.tail(point_num).array().abs().maxCoeff() << "\n";
//	std::cout << "Phase:\n" << V_reg_dir.tail(point_num).array().arg().minCoeff() << "\t" << V_reg_dir.tail(point_num).array().arg().maxCoeff() << "\n\n";
//
//	// Store results
//	Power_network_inform.power_flow.voltage.row(tick) = V_reg_dir.array();
//	Eigen::VectorXcd power_source = V_reg_dir.array() * (Power_network_inform.power_flow.nodal_admittance * V_reg_dir).array();
//	std::cout << Power_network_inform.power_flow.power_node(tick, 0) << "\t" << power_source(0) << "\n\n";
//	std::cout << Power_network_inform.power_flow.power_node(tick, 1) << "\t" << power_source(1) << "\n\n";
}
