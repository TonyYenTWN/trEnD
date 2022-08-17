// Source file for power network data input and process
#include <iostream>
#include "src/spatial_field/geostat.h"
#include "power_network.h"

// In-file functions
namespace {
	void cbt_data_input(power_network::network_inform &Power_network_inform, std::string fin_cbt, std::string fin_entry){
		// Read cross-border transmission data
		auto fin_cbt_dim = basic::get_file_dim(fin_cbt);
		auto fin_entry_dim = basic::get_file_dim(fin_entry);

		// Store bidding zone names
		Power_network_inform.cbt.bz_names = basic::get_col_name(fin_cbt, fin_cbt_dim[1]);

		// Store entry points
		Power_network_inform.cbt.entry_nodes = basic::read_file(fin_entry_dim[0], fin_entry_dim[1], fin_entry);
		Power_network_inform.cbt.entry_nodes = Power_network_inform.cbt.entry_nodes.rightCols(fin_entry_dim[1] - 1);
		Power_network_inform.cbt.entry_nodes = Power_network_inform.cbt.entry_nodes.array().max(0) - 1.;

		// Find number of entry nodes for each bidding zones
		Power_network_inform.cbt.entry_node_num = Eigen::VectorXi::Zero(fin_entry_dim[0]);
		for(int zone_iter = 0; zone_iter < fin_entry_dim[0]; ++ zone_iter){
			for(int node_iter = 0; node_iter < fin_entry_dim[1] - 1; ++ node_iter){
				if(Power_network_inform.cbt.entry_nodes(zone_iter, node_iter) < 0.){
					break;
				}
				Power_network_inform.cbt.entry_node_num(zone_iter) += 1;
			}
		}

		// Store cbt constraints
		Power_network_inform.cbt.flow_constraint = basic::read_file(fin_cbt_dim[0], fin_cbt_dim[1], fin_cbt);
	}

	// Must read transmission data before points (DSO cluster initialize here)
	void tranmission_data_input(power_network::network_inform &Power_network_inform, Eigen::MatrixXd bz_inform, std::string fin_node, std::string fin_edge, std::string fin_edge_simp){
		// Read power network data
		auto fin_node_dim = basic::get_file_dim(fin_node);
		auto fin_edge_dim = basic::get_file_dim(fin_edge);
		auto fin_edge_simp_dim = basic::get_file_dim(fin_edge_simp);
		auto node_inform = basic::read_file(fin_node_dim[0], fin_node_dim[1], fin_node);
		auto edge_inform = basic::read_file(fin_edge_dim[0], fin_edge_dim[1], fin_edge);
		auto edge_simp_inform = basic::read_file(fin_edge_simp_dim[0], fin_edge_simp_dim[1], fin_edge_simp);

		// Initialize node ID for DSO-Clusters
		int cluster_num = int(node_inform.col(1).maxCoeff());
		Power_network_inform.DSO_cluster.clear();
		Power_network_inform.DSO_cluster = std::vector <power_network::DSO_cluster> (cluster_num);
		for(int cluster_iter = 0; cluster_iter < cluster_num; ++ cluster_iter){
			Power_network_inform.DSO_cluster[cluster_iter].nodes_ID.clear();
			Power_network_inform.DSO_cluster[cluster_iter].nodes_ID.reserve(fin_node_dim[0]);
		}

		// Organize node cata
		Power_network_inform.nodes.bidding_zone = Eigen::VectorXi(fin_node_dim[0]);
		Power_network_inform.nodes.cluster = Eigen::VectorXi(fin_node_dim[0]);
		Power_network_inform.nodes.voltage_base = Eigen::VectorXi(fin_node_dim[0]);
		Power_network_inform.nodes.in_cluster_ID = Eigen::VectorXi(fin_node_dim[0]);
		for(int node_iter = 0; node_iter < fin_node_dim[0]; ++ node_iter){
			Power_network_inform.nodes.bidding_zone(node_iter) = bz_inform(int(node_inform(node_iter, 0)) - 1, 1) - 1;
			Power_network_inform.nodes.cluster(node_iter) = int(node_inform(node_iter, 1)) - 1;
			Power_network_inform.nodes.voltage_base(node_iter) = int(node_inform(node_iter, 2));
			Power_network_inform.DSO_cluster[Power_network_inform.nodes.cluster(node_iter)].nodes_ID.push_back(node_iter);
			Power_network_inform.nodes.in_cluster_ID(node_iter) = Power_network_inform.DSO_cluster[Power_network_inform.nodes.cluster(node_iter)].nodes_ID.size() - 1;
		}
		Power_network_inform.nodes.x = node_inform.col(node_inform.cols() - 4);
		Power_network_inform.nodes.y = node_inform.col(node_inform.cols() - 3);
		Power_network_inform.nodes.lon = node_inform.col(node_inform.cols() - 2);
		Power_network_inform.nodes.lat = node_inform.col(node_inform.cols() - 1);

		// Organize edge data
		Power_network_inform.edges.from = Eigen::VectorXi(fin_edge_dim[0]);
		Power_network_inform.edges.to = Eigen::VectorXi(fin_edge_dim[0]);
		Power_network_inform.edges.voltage_base = Eigen::VectorXi(fin_edge_dim[0]);
		for(int edge_iter = 0; edge_iter < fin_edge_dim[0]; ++ edge_iter){
			Power_network_inform.edges.from(edge_iter) = int(edge_inform(edge_iter, 0)) - 1;
			Power_network_inform.edges.to(edge_iter) = int(edge_inform(edge_iter, 1)) - 1;
			Power_network_inform.edges.voltage_base(edge_iter) = int(edge_inform(edge_iter, 4));
		}
		Power_network_inform.edges.distance = edge_inform.col(5);
	}

	void points_data_input(power_network::network_inform &Power_network_inform, Eigen::MatrixXd bz_inform, std::string fin_point, std::string fin_point_matrix){
		// Read point data
		auto fin_point_dim = basic::get_file_dim(fin_point);
		auto fin_point_matrix_dim = basic::get_file_dim(fin_point_matrix);
		auto point_inform = basic::read_file(fin_point_dim[0], fin_point_dim[1], fin_point);
		auto point_matrix = basic::read_file(fin_point_matrix_dim[0], fin_point_matrix_dim[1], fin_point_matrix);

		// Initialize point ID for DSO-Clusters
		for(int cluster_iter = 0; cluster_iter < Power_network_inform.DSO_cluster.size(); ++ cluster_iter){
			Power_network_inform.DSO_cluster[cluster_iter].points_ID.clear();
			Power_network_inform.DSO_cluster[cluster_iter].points_ID.reserve(fin_point_dim[0]);
		}

		// Organize point data
		Power_network_inform.points.bidding_zone = Eigen::VectorXi(fin_point_dim[0]);
		Power_network_inform.points.node = Eigen::VectorXi(fin_point_dim[0]);
		Power_network_inform.points.in_cluster_ID = Eigen::VectorXi(fin_point_dim[0]);
		for(int point_iter = 0; point_iter < fin_point_dim[0]; ++ point_iter){
			int node_ID = int(point_inform(point_iter, 1)) - 1;
			int bz_ID = bz_inform(int(point_inform(point_iter, 0)) - 1, 1) - 1;
			Power_network_inform.points.node(point_iter) = node_ID;
			Power_network_inform.points.bidding_zone(point_iter) = bz_ID;
			Power_network_inform.DSO_cluster[Power_network_inform.nodes.cluster(Power_network_inform.points.node(point_iter))].points_ID.push_back(point_iter);
			Power_network_inform.points.in_cluster_ID(point_iter) = Power_network_inform.DSO_cluster[Power_network_inform.nodes.cluster(Power_network_inform.points.node(point_iter))].points_ID.size() - 1;
		}
		Power_network_inform.points.population_density = point_inform.col(2);
		Power_network_inform.points.x = point_inform.col(point_inform.cols() - 4);
		Power_network_inform.points.y = point_inform.col(point_inform.cols() - 3);
		Power_network_inform.points.lon = point_inform.col(point_inform.cols() - 2);
		Power_network_inform.points.lat = point_inform.col(point_inform.cols() - 1);

		// Read coordinate grid data
		Power_network_inform.points.coordinate_grid = Eigen::MatrixXi(fin_point_matrix_dim[0], fin_point_matrix_dim[1]);
		for(int x_coor_iter = 0; x_coor_iter < fin_point_matrix_dim[0]; ++ x_coor_iter){
			for(int y_coor_iter = 0; y_coor_iter < fin_point_matrix_dim[1]; ++ y_coor_iter){
				Power_network_inform.points.coordinate_grid(x_coor_iter, y_coor_iter) = std::max(int(point_matrix(x_coor_iter, y_coor_iter)), 0) - 1;
			}
		}

		// Calculate distance matrix
		power_network::point_distance_cov(Power_network_inform.points, 10.);
		//std::cout << Power_network_inform.points.distance.bottomRightCorner(5, 5) << "\n\n";
	}

	void plant_data_input(power_network::network_inform &Power_network_inform, std::string fin_hydro, std::string fin_wind){
		// Read plant data
		auto fin_hydro_dim = basic::get_file_dim(fin_hydro);
		auto fin_wind_dim = basic::get_file_dim(fin_wind);
		auto hydro_inform = basic::read_file(fin_hydro_dim[0], fin_hydro_dim[1], fin_hydro);
		auto wind_inform = basic::read_file(fin_wind_dim[0], fin_wind_dim[1], fin_wind);

		// Organize hydro power plant data
		Power_network_inform.plants.hydro.node = Eigen::VectorXi(fin_hydro_dim[0]);
		Power_network_inform.plants.hydro.type = Eigen::VectorXi(fin_hydro_dim[0]);
		for(int hydro_iter = 0; hydro_iter < fin_hydro_dim[0]; ++ hydro_iter){
			Power_network_inform.plants.hydro.node(hydro_iter) = int(hydro_inform(hydro_iter, 0)) - 1;
			Power_network_inform.plants.hydro.type(hydro_iter) = int(hydro_inform(hydro_iter, 1)) - 1;
		}
		Power_network_inform.plants.hydro.cap = hydro_inform.col(3).array().abs();
		Power_network_inform.plants.hydro.x = hydro_inform.col(hydro_inform.cols() - 4);
		Power_network_inform.plants.hydro.y = hydro_inform.col(hydro_inform.cols() - 3);
		Power_network_inform.plants.hydro.lon = hydro_inform.col(hydro_inform.cols() - 2);
		Power_network_inform.plants.hydro.lat = hydro_inform.col(hydro_inform.cols() - 1);

		// Organize wind power plant data
		Power_network_inform.plants.wind.node = Eigen::VectorXi(fin_wind_dim[0]);
		for(int wind_iter = 0; wind_iter < fin_wind_dim[0]; ++ wind_iter){
			Power_network_inform.plants.wind.node(wind_iter) = int(wind_inform(wind_iter, 0)) - 1;
		}
		Power_network_inform.plants.wind.cap = wind_inform.col(3);
		Power_network_inform.plants.wind.x = wind_inform.col(wind_inform.cols() - 4);
		Power_network_inform.plants.wind.y = wind_inform.col(wind_inform.cols() - 3);
		Power_network_inform.plants.wind.lon = wind_inform.col(wind_inform.cols() - 2);
		Power_network_inform.plants.wind.lat = wind_inform.col(wind_inform.cols() - 1);

		//std::cout << Power_network_inform.plants.wind.node.tail(10).transpose() << "\n";
		//std::cout << Power_network_inform.plants.wind.lat.tail(10).transpose() << "\n";
	}
}

// Find the distance and covariance between points
void power_network::point_distance_cov(points &point, double lambda){
	double pi = boost::math::constants::pi<double>();
	double tol = 1E-8;
	Eigen::MatrixXd covariance = Eigen::MatrixXd::Ones(point.lon.size(), point.lon.size());
	point.distance = Eigen::MatrixXd::Zero(point.lon.size(), point.lon.size());

	#pragma omp parallel
	{
		#pragma omp for
		for(int row_iter = 0; row_iter < point.lon.size() - 1; ++ row_iter){
			for(int col_iter = row_iter + 1; col_iter < point.lon.size(); ++ col_iter){
				point.distance(row_iter, col_iter) = spatial_field::geodist(Eigen::Vector2d(point.lon(row_iter) * pi / 180., point.lat(row_iter) * pi / 180.), Eigen::Vector2d(point.lon(col_iter) * pi / 180., point.lat(col_iter) * pi / 180.));
				point.distance(col_iter, row_iter) = point.distance(row_iter, col_iter);
				covariance(row_iter, col_iter) = exp(-pow(point.distance(row_iter, col_iter), 1.) / point.grid_length / lambda);
				covariance(col_iter, row_iter) = covariance(row_iter, col_iter);
			}
		}
	}

	point.covariance = covariance.sparseView(tol);
}

void power_network::power_network_input_process(network_inform &Power_network_inform, std::string parent_dir){
	auto fin_bz = parent_dir + "DSO_Bidding_Zone.csv";
	auto fin_cbt = parent_dir + "cbt_constraint.csv";
	auto fin_entry = parent_dir + "cbt_entry_nodes.csv";
	auto fin_node = parent_dir + "transmission_nodes.csv";
	auto fin_edge = parent_dir + "transmission_edges.csv";
	auto fin_edge_simp = parent_dir + "transmission_edges_pu_simp.csv";
	auto fin_point = parent_dir + "point_info.csv";
	auto fin_point_matrix = parent_dir + "point_matrix.csv";
	auto fin_hydro = parent_dir + "hydro_plants.csv";
	auto fin_wind = parent_dir + "wind_plants.csv";

	auto fin_bz_dim = basic::get_file_dim(fin_bz);
	auto bz_inform = basic::read_file(fin_bz_dim[0], fin_bz_dim[1], fin_bz);

	cbt_data_input(Power_network_inform, fin_cbt, fin_entry);
	tranmission_data_input(Power_network_inform, bz_inform, fin_node, fin_edge, fin_edge_simp);
	points_data_input(Power_network_inform, bz_inform, fin_point, fin_point_matrix);
	plant_data_input(Power_network_inform, fin_hydro, fin_wind);

	// Set power line density for the distribution network
	Power_network_inform.set_line_density();

	// Construct base voltage and base impedance level maps
	Power_network_inform.tech_parameters.set_level_maps(Power_network_inform.nodes);
}

//int main(){
//	network_inform Power_network_inform;
//	power_network_input_process(Power_network_inform);
//};
