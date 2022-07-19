// Source file for dispatch filtering of DSOs in Norway
// Note that in Norwegian terminology, HV distribution network (> 33kV) is called regional network.
#include <iostream>
//#include <chrono>
//#include "../basic/LP_gpa.h"
#include "src/spatial_field/geostat.h"
#include "power_market.h"

void power_market::DSO_Markets_Set(markets_inform &DSO_Markets, power_network::network_inform &Power_network_inform, int Time){
	double pi = boost::math::constants::pi<double>();
	DSO_Markets.clear();
	DSO_Markets = std::vector <market_inform> (Power_network_inform.DSO_cluster.size());

	// Initialize markets at each clustered DSO
	for(int DSO_iter = 0; DSO_iter < DSO_Markets.size(); ++ DSO_iter){
		// Input parameters of DSO market
		DSO_Markets[DSO_iter].num_zone = Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() + Power_network_inform.DSO_cluster[DSO_iter].nodes_ID.size();
		DSO_Markets[DSO_iter].time_intervals = Time;
		DSO_Markets[DSO_iter].set_bidded_price();

		// Set node admittance matrix
		// Use fractional Laplacian here :
		// y(x_1, x_2) = N / 4 / pi / z / |x_1 - x_2|^(2 + d), where
		// y(x_1, x_2): admittance between x_1 and x_2
		// N: mean line density at the neighborhood of x_1 containing also x_2
		// z: per length impedence of power lines
		// d: fractional dimension of the power line distribution
		double z_base_low = pow(Power_network_inform.tech_parameters.voltage_cutoff_distr, 2.) / Power_network_inform.tech_parameters.s_base / 3.;
		double z_base_high = pow(66., 2.) / Power_network_inform.tech_parameters.s_base / 3.;
		double partition_func = 0.;
		DSO_Markets[DSO_iter].network.num_vertice = DSO_Markets[DSO_iter].num_zone;
		Eigen::MatrixXd admittance = Eigen::MatrixXd::Ones(DSO_Markets[DSO_iter].network.num_vertice, DSO_Markets[DSO_iter].network.num_vertice);
		Eigen::MatrixXd distance = Eigen::MatrixXd::Zero(DSO_Markets[DSO_iter].network.num_vertice, DSO_Markets[DSO_iter].network.num_vertice);

//		std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
//		std::cout << DSO_iter << "\n";
//		std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
		for(int row_iter = 0; row_iter < DSO_Markets[DSO_iter].network.num_vertice - 1; ++ row_iter){
			for(int col_iter = row_iter + 1; col_iter < DSO_Markets[DSO_iter].network.num_vertice; ++ col_iter){
				if(row_iter < Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() && col_iter < Power_network_inform.DSO_cluster[DSO_iter].points_ID.size()){
					int point_ID_1 = Power_network_inform.DSO_cluster[DSO_iter].points_ID[row_iter];
					int point_ID_2 = Power_network_inform.DSO_cluster[DSO_iter].points_ID[col_iter];
					distance(row_iter, col_iter) = Power_network_inform.points.distance(point_ID_1, point_ID_2);
					admittance(row_iter, col_iter) = 1. / pow(distance(row_iter, col_iter) * 1E-3, 1. + Power_network_inform.tech_parameters.fraction_dim_distr);
					partition_func += admittance(row_iter, col_iter);
				}
				else{
					if(row_iter < Power_network_inform.DSO_cluster[DSO_iter].points_ID.size()){
						int point_ID = Power_network_inform.DSO_cluster[DSO_iter].points_ID[row_iter];
						int node_ID = Power_network_inform.DSO_cluster[DSO_iter].nodes_ID[col_iter - Power_network_inform.DSO_cluster[DSO_iter].points_ID.size()];
						Eigen::Vector2d point_coor = Eigen::Vector2d(Power_network_inform.points.lon(point_ID), Power_network_inform.points.lat(point_ID));
						Eigen::Vector2d node_coor = Eigen::Vector2d(Power_network_inform.nodes.lon(node_ID), Power_network_inform.nodes.lat(node_ID));
						point_coor *= pi / 180.;
						node_coor *= pi / 180.;
						distance(row_iter, col_iter) = spatial_field::geodist(point_coor, node_coor);
					}
					else{
						int node_ID_1 = Power_network_inform.DSO_cluster[DSO_iter].nodes_ID[row_iter - Power_network_inform.DSO_cluster[DSO_iter].points_ID.size()];
						int node_ID_2 = Power_network_inform.DSO_cluster[DSO_iter].nodes_ID[col_iter - Power_network_inform.DSO_cluster[DSO_iter].points_ID.size()];
						Eigen::Vector2d node_coor_1 = Eigen::Vector2d(Power_network_inform.nodes.lon(node_ID_1), Power_network_inform.nodes.lat(node_ID_1));
						Eigen::Vector2d node_coor_2 = Eigen::Vector2d(Power_network_inform.nodes.lon(node_ID_2), Power_network_inform.nodes.lat(node_ID_2));
						node_coor_1 *= pi / 180.;
						node_coor_2 *= pi / 180.;
						distance(row_iter, col_iter) = spatial_field::geodist(node_coor_1, node_coor_2);
					}
				}
			}
		}
		admittance.topLeftCorner(Power_network_inform.DSO_cluster[DSO_iter].points_ID.size(), Power_network_inform.DSO_cluster[DSO_iter].points_ID.size()) *= Power_network_inform.tech_parameters.line_density_distr * Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() * z_base_low;
		admittance.topLeftCorner(Power_network_inform.DSO_cluster[DSO_iter].points_ID.size(), Power_network_inform.DSO_cluster[DSO_iter].points_ID.size()) /= partition_func * Power_network_inform.tech_parameters.z_distr_series.imag();
		admittance.rightCols(Power_network_inform.DSO_cluster[DSO_iter].nodes_ID.size()) *= z_base_high / Power_network_inform.tech_parameters.z_distr_series.imag();
		int demo = std::min(15, DSO_Markets[DSO_iter].network.num_vertice);
//		std::cout << admittance.topLeftCorner(demo, demo) << "\n\n";
//		std::cout << distance.topLeftCorner(demo, demo) << "\n\n";
		admittance = admittance.array() / distance.array();
//		std::cout << admittance.topLeftCorner(demo, demo) << "\n\n";
//		std::cout << "\n\n";

		// Set compact incidence matrix and edge admittance matrix
		double tol = 1.;
		DSO_Markets[DSO_iter].network.incidence.reserve(DSO_Markets[DSO_iter].network.num_vertice * DSO_Markets[DSO_iter].network.num_vertice);
		DSO_Markets[DSO_iter].network.admittance.reserve(DSO_Markets[DSO_iter].network.num_vertice * DSO_Markets[DSO_iter].network.num_vertice);

		for(int row_iter = 0; row_iter < DSO_Markets[DSO_iter].network.num_vertice - 1; ++ row_iter){
			for(int col_iter = row_iter + 1; col_iter < DSO_Markets[DSO_iter].network.num_vertice; ++ col_iter){
				if(row_iter < Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() && col_iter < Power_network_inform.DSO_cluster[DSO_iter].points_ID.size()){
					if(admittance(row_iter, col_iter) > tol){
						DSO_Markets[DSO_iter].network.incidence.push_back(Eigen::Vector2i(row_iter, col_iter));
						DSO_Markets[DSO_iter].network.admittance.push_back(admittance(row_iter , col_iter));
						// power limit?
					}
				}
				else{
					DSO_Markets[DSO_iter].network.incidence.push_back(Eigen::Vector2i(row_iter, col_iter));
					DSO_Markets[DSO_iter].network.admittance.push_back(admittance(row_iter , col_iter));
					// power limit?
				}
			}
		}

//		DSO_Markets[DSO_iter].network.num_edges = DSO_Markets[DSO_iter].num_zone * (DSO_Markets[DSO_iter].num_zone - 1) / 2;
//		if(DSO_Markets[DSO_iter].num_zone > 1){
//			DSO_Markets[DSO_iter].network.incidence_matrix = Eigen::MatrixXi(DSO_Markets[DSO_iter].network.num_edges, 2);
//			DSO_Markets[DSO_iter].network.admittance_vector = Eigen::VectorXd(DSO_Markets[DSO_iter].network.num_edges);
//			int edge_iter = 0;
//			for(int row_iter = 0; row_iter < DSO_Markets[DSO_iter].network.num_vertice - 1; ++ row_iter){
//				for(int col_iter = row_iter + 1; col_iter < DSO_Markets[DSO_iter].network.num_vertice; ++ col_iter){
//					DSO_Markets[DSO_iter].network.incidence_matrix(edge_iter, 0) = row_iter;
//					DSO_Markets[DSO_iter].network.incidence_matrix(edge_iter, 1) = col_iter;
//					DSO_Markets[DSO_iter].network.admittance_vector(edge_iter) = Power_network_inform.tech_parameters.line_density_distr / 4. / pi;
//					DSO_Markets[DSO_iter].network.admittance_vector(edge_iter) /= Power_network_inform.tech_parameters.z_distr_series.imag() * 1E3;
//					double distance;
//
//					if(row_iter < Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() && col_iter < Power_network_inform.DSO_cluster[DSO_iter].points_ID.size()){
//						int point_ID_1 = Power_network_inform.DSO_cluster[DSO_iter].points_ID[row_iter];
//						int point_ID_2 = Power_network_inform.DSO_cluster[DSO_iter].points_ID[col_iter];
//						distance = Power_network_inform.points.distance(point_ID_1, point_ID_2);
//					}
//					else{
//						if(row_iter < Power_network_inform.DSO_cluster[DSO_iter].points_ID.size()){
//							int point_ID = Power_network_inform.DSO_cluster[DSO_iter].points_ID[row_iter];
//							int node_ID = Power_network_inform.DSO_cluster[DSO_iter].nodes_ID[col_iter - Power_network_inform.DSO_cluster[DSO_iter].points_ID.size()];
//							Eigen::Vector2d point_coor = Eigen::Vector2d(Power_network_inform.points.lon(point_ID), Power_network_inform.points.lat(point_ID));
//							Eigen::Vector2d node_coor = Eigen::Vector2d(Power_network_inform.nodes.lon(point_ID), Power_network_inform.nodes.lat(point_ID));
//							point_coor *= pi / 180.;
//							node_coor *= pi / 180.;
//							distance = spatial_field::geodist(point_coor, node_coor);
//						}
//						else{
//							int node_ID_1 = Power_network_inform.DSO_cluster[DSO_iter].nodes_ID[row_iter - Power_network_inform.DSO_cluster[DSO_iter].points_ID.size()];
//							int node_ID_2 = Power_network_inform.DSO_cluster[DSO_iter].nodes_ID[col_iter - Power_network_inform.DSO_cluster[DSO_iter].points_ID.size()];
//							Eigen::Vector2d node_coor_1 = Eigen::Vector2d(Power_network_inform.points.lon(node_ID_1), Power_network_inform.points.lat(node_ID_1));
//							Eigen::Vector2d node_coor_2 = Eigen::Vector2d(Power_network_inform.nodes.lon(node_ID_2), Power_network_inform.nodes.lat(node_ID_2));
//							node_coor_1 *= pi / 180.;
//							node_coor_2 *= pi / 180.;
//							distance = spatial_field::geodist(node_coor_1, node_coor_2);
//						}
//					}
//					DSO_Markets[DSO_iter].network.admittance_vector(edge_iter) /= pow(distance * 1E-3, 2. + Power_network_inform.tech_parameters.fraction_dim_distr);
//					edge_iter += 1;
//				}
//			}
//		}

		// Set voltage and power constraints at each edge
//		DSO_Markets[DSO_iter].network.voltage_constraint = Eigen::MatrixXd::Ones(DSO_Markets[DSO_iter].network.num_vertice, 2);
//		DSO_Markets[DSO_iter].network.voltage_constraint.col(0) *= -pi / 12;
//		DSO_Markets[DSO_iter].network.voltage_constraint.col(1) *= pi / 12;
//		DSO_Markets[DSO_iter].network.power_constraint = Eigen::MatrixXd::Ones(DSO_Markets[DSO_iter].network.num_edges, 2);
//		DSO_Markets[DSO_iter].network.power_constraint.col(0) *= -50.;
//		DSO_Markets[DSO_iter].network.power_constraint.col(1) *= 50.;

		// Initialization of process variables
		power_market::Market_Initialization(DSO_Markets[DSO_iter]);

		// Initialization of output variables
		DSO_Markets[DSO_iter].confirmed_supply = Eigen::MatrixXd::Zero(Time, DSO_Markets[DSO_iter].num_zone);
		DSO_Markets[DSO_iter].confirmed_demand = Eigen::MatrixXd::Zero(Time, DSO_Markets[DSO_iter].num_zone);
		DSO_Markets[DSO_iter].confirmed_price = Eigen::MatrixXd(Time, DSO_Markets[DSO_iter].num_zone);
		DSO_Markets[DSO_iter].network.confirmed_power = Eigen::MatrixXd(Time, DSO_Markets[DSO_iter].network.num_edges);
		DSO_Markets[DSO_iter].network.confirmed_voltage = Eigen::MatrixXd(Time, DSO_Markets[DSO_iter].network.num_vertice);
	}
}

//int main(){
//	network_inform Power_network_inform;
//	power_network_input_process(Power_network_inform);
//
//	int Time = 8760;
//	std::string fin_name_moc = "input/merit_order_curve_q_assimilated_2021.csv";
//	std::string fin_name_demand = "input/residual_load_default_forecast_2021.csv";
//	market_inform International_Market;
//	International_Market_Set(International_Market, Time, fin_name_moc, fin_name_demand);
//
//	market_inform TSO_Market;
//	LP_object TSO_Problem;
//	TSO_Market_Set(TSO_Market, Power_network_inform, 1);
//	Flow_Based_Market_LP_Set(TSO_Market, TSO_Problem);
//
//	std::string fin_point_demand = "../area_deamand/processed/nominal_mean_demand_field_10km_annual_mean.csv";
//	DSO_Markets DSO_Markets;
//	DSO_Markets_Set(DSO_Markets, Power_network_inform, 1);
//	Submitted_bid_calculation(0, DSO_Markets, TSO_Market, International_Market, Power_network_inform, fin_point_demand);
//}
