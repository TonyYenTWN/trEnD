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

void power_network::HELM_Set(network_inform &Power_network_inform){
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

		// Triplet for series impedence
		Y_n_trip.push_back(Eigen::TripletXcd(from_ID, to_ID, -y_series));
		Y_n_trip.push_back(Eigen::TripletXcd(to_ID, from_ID, -y_series));

		// Update diagonal terms
		Y_n_Diag(from_ID) += y_series + .5 * y_shunt;
		Y_n_Diag(to_ID) += y_series + .5 * y_shunt;
	}

	// Distribution level
	for(int DSO_iter = 0; DSO_iter < DSO_num; ++ DSO_iter){
		int DSO_point_num = Power_network_inform.DSO_cluster[DSO_iter].points_ID.size();
		int DSO_node_num = Power_network_inform.DSO_cluster[DSO_iter].nodes_ID.size();

		double z_base_low = pow(Power_network_inform.tech_parameters.voltage_cutoff_distr, 2.) / Power_network_inform.tech_parameters.s_base * 3.;
		double z_base_high = pow(Power_network_inform.tech_parameters.voltage_cutoff_connection, 2.) / Power_network_inform.tech_parameters.s_base * 3.;
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
			y_series /= Power_network_inform.tech_parameters.z_distr_series;
			y_series *= z_base_high;
			y_series *= Power_network_inform.tech_parameters.line_density_connection * DSO_point_num / DSO_node_num;

			// Shunt admittance
			std::complex <double> y_shunt(1., 0.);
			y_shunt *= distance_min;
			y_shunt *= Power_network_inform.tech_parameters.y_distr_shunt;
			y_shunt *= z_base_high;
			y_shunt *= Power_network_inform.tech_parameters.line_density_connection * DSO_point_num / DSO_node_num;

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

	// -------------------------------------------------------------------------------
	// Store the nodal admittance matrix and intialize
	// -------------------------------------------------------------------------------
	int Time = configuration::parameters::Time();

	Power_network_inform.power_flow.nodal_admittance = Eigen::SparseMatrix <std::complex <double>> (node_num + point_num, node_num + point_num);
	Power_network_inform.power_flow.nodal_admittance.setFromTriplets(Y_n_trip.begin(), Y_n_trip.end());
	Power_network_inform.power_flow.power_edge = Eigen::MatrixXcd::Zero(Time, edge_trans_num);
	Power_network_inform.power_flow.power_node = Eigen::MatrixXcd::Zero(Time, node_num + point_num);
	Power_network_inform.power_flow.voltage = Eigen::MatrixXcd::Zero(Time, node_num + point_num);
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
			Power_network_inform.power_flow.power_node(tick, node_ID) += real_power;
		}
	}

	// Update from end-users
	{
		int sample_num = agent::end_user::parameters::sample_num();
		std::complex <double> apparent_power_base(1., pow(1. - pow(agent::end_user::parameters::power_factor(), 2.), .5));
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				double real_power;
				real_power = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_supply;
				Power_network_inform.power_flow.power_node(tick, node_num + point_iter) += real_power;
				real_power = -Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_demand;
				Power_network_inform.power_flow.power_node(tick, node_num + point_iter) += real_power * apparent_power_base;
			}
		}
	}

	// Update from industrial demand
	{
		int industrial_HV_num = Power_market_inform.agent_profiles.industrial.HV.size();
		std::complex <double> apparent_power_base(1., pow(1. - pow(agent::industrial::parameters::power_factor(), 2.), .5));
		for(int agent_iter = 0; agent_iter < industrial_HV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.industrial.HV[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);

			double real_power = -Power_market_inform.agent_profiles.industrial.HV[agent_iter].results.actual_demand;
			Power_network_inform.power_flow.power_node(tick, node_ID) += real_power * apparent_power_base;
		}
	}

	// Update from power suppliers
	int hydro_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		double real_power = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		Power_network_inform.power_flow.power_node(tick, node_ID) += real_power;
	}

	int hydro_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].point_ID;

		double real_power = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		Power_network_inform.power_flow.power_node(tick, node_num + point_ID) += real_power;
	}

	int wind_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		double real_power = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		Power_network_inform.power_flow.power_node(tick, node_ID) += real_power;
	}

	int wind_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].point_ID;

		double real_power = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		Power_network_inform.power_flow.power_node(tick, node_num + point_ID) += real_power;
	}

	int pump_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV.size();
	for(int agent_iter = 0; agent_iter < pump_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		double real_power = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		real_power -= Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].results.actual_demand;
		Power_network_inform.power_flow.power_node(tick, node_ID) += real_power;
	}

	int pump_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
	for(int agent_iter = 0; agent_iter < pump_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;

		double real_power = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		real_power -= Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].results.actual_demand;
		Power_network_inform.power_flow.power_node(tick, node_num + point_ID) += real_power;
	}
}

void power_network::HELM_Solve(int tick, network_inform &Power_network_inform){
	int node_num = Power_network_inform.nodes.bidding_zone.size();
	int point_num = Power_network_inform.points.bidding_zone.size();
	auto Y_n = Power_network_inform.power_flow.nodal_admittance;

	// ---------------------------------------------------------------------------------------------
	// Set the matrix for the system of iterative linear equations
	// ---------------------------------------------------------------------------------------------
	// Assume all PQ Nodes
	// Order of variables:
	// {V}, {1 / V}, {\hat V}, {1 / \hat V}
	std::vector<Eigen::TripletXcd> Mat_trip_const;
	Mat_trip_const.reserve(2 * Power_network_inform.power_flow.nodal_admittance.nonZeros() + 6 * (node_num + point_num));

	// Entries from original nodal admittance matrix
	for(int col_iter = 0; col_iter < Y_n.outerSize(); ++ col_iter){
		for(Eigen::SparseMatrix<std::complex <double>>::InnerIterator inner_iter(Y_n, col_iter); inner_iter; ++ inner_iter){
			if(inner_iter.row() == 0){
				continue;
			}

			std::complex <double> y_conj = inner_iter.value();
			y_conj = std::conj(y_conj);
			Mat_trip_const.push_back(Eigen::TripletXcd(inner_iter.row() - 1, inner_iter.col(), inner_iter.value()));
			Mat_trip_const.push_back(Eigen::TripletXcd(2 * (node_num + point_num) + inner_iter.row() - 1, 2 * (node_num + point_num) + inner_iter.col(), y_conj));
		}
	}

	// Flow conservation law
	for(int bus_iter = 0; bus_iter < node_num + point_num; ++ bus_iter){
		int row_ID;
		int col_ID;
		std::complex <double> s_bus = Power_network_inform.power_flow.power_node(tick, bus_iter);

		row_ID = node_num + point_num - 1;
		col_ID = node_num + point_num + bus_iter;
		Mat_trip_const.push_back(Eigen::TripletXcd(row_ID, col_ID, s_bus));

		row_ID += 2 * (node_num + point_num);
		col_ID += 2 * (node_num + point_num);
		Mat_trip_const.push_back(Eigen::TripletXcd(row_ID, col_ID, std::conj(s_bus)));
	}

	// -------------------------------------------------------------------------------
	// Initialization of power series coefficients
	// -------------------------------------------------------------------------------
	int power_terms = 20;
	Eigen::MatrixXcd V_up_reg = Eigen::MatrixXcd::Zero(node_num + point_num, power_terms);
	Eigen::MatrixXcd V_up_hat = Eigen::MatrixXcd::Zero(node_num + point_num, power_terms);
	Eigen::MatrixXcd V_down_reg = Eigen::MatrixXcd::Zero(node_num + point_num, power_terms);
	Eigen::MatrixXcd V_down_hat = Eigen::MatrixXcd::Zero(node_num + point_num, power_terms);

	// -------------------------------------------------------------------------------
	// Trivial solution (V = 1 everywhere)
	// -------------------------------------------------------------------------------
	V_up_reg.col(0) = Eigen::VectorXcd::Ones(node_num + point_num);
	V_up_hat.col(0) = Eigen::VectorXcd::Ones(node_num + point_num);
	V_down_reg.col(0) = Eigen::VectorXcd::Ones(node_num + point_num);
	V_down_hat.col(0) = Eigen::VectorXcd::Ones(node_num + point_num);

	// -------------------------------------------------------------------------------
	// Solve linear system for each iteration
	// -------------------------------------------------------------------------------
	int term_iter = 1;
	std::vector<Eigen::TripletXcd> Mat_trip_temp;
	Mat_trip_temp = Mat_trip_const;
	Mat_trip_temp.reserve(2 * Power_network_inform.power_flow.nodal_admittance.nonZeros() + 6 * (node_num + point_num));
	// Update reciprocal relation
	for(int bus_iter = 0; bus_iter < node_num + point_num; ++ bus_iter){
		int row_ID;
		int col_ID;

		row_ID = node_num + point_num + bus_iter;
		col_ID = bus_iter;
		Mat_trip_temp.push_back(Eigen::TripletXcd(row_ID, col_ID, V_down_reg(bus_iter, term_iter - 1)));
		col_ID += node_num + point_num;
		Mat_trip_temp.push_back(Eigen::TripletXcd(row_ID, col_ID, V_up_reg(bus_iter, term_iter - 1)));

		row_ID += 2 * (node_num + point_num);
		col_ID += node_num + point_num;
		Mat_trip_temp.push_back(Eigen::TripletXcd(row_ID, col_ID, V_down_hat(bus_iter, term_iter - 1)));
		col_ID += node_num + point_num;
		Mat_trip_temp.push_back(Eigen::TripletXcd(row_ID, col_ID, V_up_hat(bus_iter, term_iter - 1)));
	}

	// Update rhs of the equation
	Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(4* (node_num + point_num));
	rhs.head(node_num + point_num - 1) = Power_network_inform.power_flow.power_node.row(tick).tail(node_num + point_num - 1).conjugate().transpose();
	rhs.head(node_num + point_num - 1) = rhs.head(node_num + point_num - 1).array() / V_down_hat.col(term_iter - 1).tail(node_num + point_num - 1).array();
	rhs.segment(2 * (node_num + point_num), node_num + point_num - 1) = Power_network_inform.power_flow.power_node.row(tick).tail(node_num + point_num - 1).transpose();
	rhs.segment(2 * (node_num + point_num), node_num + point_num - 1) = rhs.segment(2 * (node_num + point_num), node_num + point_num - 1).array() / V_down_reg.col(term_iter - 1).tail(node_num + point_num - 1).array();
}
