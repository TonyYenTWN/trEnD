// Source file for power flow analysis
#include "src/spatial_field/geostat.h"
#include "power_network.h"
//#include "src/power_market/power_market.h"

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
// {1}^T [Y] {V}(s)  = 0.
// {1}^T [Conj(Y)] {\hat V}(s) = 0.

power_network::power_flow power_network::HELM_Set(int system_type, Eigen::VectorXi node_type, network_inform &Power_network_inform){
	int node_num = Power_network_inform.nodes.bidding_zone.size();
	int point_num = Power_network_inform.points.bidding_zone.size();
	int edge_trans_num = Power_network_inform.edges.distance.size();
	int edge_distr_num = 0;
	int edge_conn_num = 0;
	int DSO_num = Power_network_inform.DSO_cluster.size();
	for(int DSO_iter = 0; DSO_iter < DSO_num; ++ DSO_iter){
		edge_distr_num += Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() * (Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() - 1) / 2;
		edge_conn_num += Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() * Power_network_inform.DSO_cluster[DSO_iter].nodes_ID.size();
	}
	double pi = boost::math::constants::pi<double>();

	power_flow power_flow;

	// -------------------------------------------------------------------------------
	// Set the nodal admittance matrix
	// -------------------------------------------------------------------------------
	power_flow.nodal_admittance = Eigen::SparseMatrix <std::complex <double>> (node_num + point_num, node_num + point_num);
	std::vector<Eigen::TripletXcd> Y_n_trip;
	Y_n_trip.reserve(node_num + point_num + 2 * (edge_trans_num + edge_conn_num + edge_distr_num));
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
			for(int col_iter = 0; col_iter < DSO_point_num ; ++ col_iter){
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
			for(int col_iter = 0; col_iter < DSO_point_num ; ++ col_iter){
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
			int min_point_ID;
			double distance_min = std::numeric_limits<double>::infinity();

			for(int col_iter = 0; col_iter < DSO_point_num ; ++ col_iter){
				int point_ID = Power_network_inform.DSO_cluster[DSO_iter].points_ID[col_iter];
				int node_ID = Power_network_inform.DSO_cluster[DSO_iter].nodes_ID[row_iter];
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
			Y_n_trip.push_back(Eigen::TripletXcd(node_ID, node_num + point_ID, -y_series));
			Y_n_trip.push_back(Eigen::TripletXcd(node_num + point_ID, node_ID, -y_series));

			// Update diagonal terms
			Y_n_Diag(node_num + point_ID) += y_series + .5 * y_shunt;
			Y_n_Diag(node_ID) += y_series + .5 * y_shunt;
		}
	}

	// Triplet for diagonal terms
	for(int node_iter = 0; node_iter < node_num + point_num; ++ node_iter){
		Y_n_trip.push_back(Eigen::TripletXcd(node_iter, node_iter, Y_n_Diag(node_iter)));
	}

	return power_flow;
}

