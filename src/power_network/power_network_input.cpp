// Source file for power network data input and process
#include <iostream>
#include "src/spatial_field/geostat.h"
#include "power_network.h"

// Must read transmission data before points (DSO cluster initialize here)
void tranmission_data_input(network_inform &Power_network_inform, Eigen::MatrixXd bz_inform, std::string fin_node, std::string fin_edge_orig, std::string fin_edge_simp){
	// Read power network data
	auto fin_node_dim = basic::get_file_dim(fin_node);
	auto fin_edge_orig_dim = basic::get_file_dim(fin_edge_orig);
	auto fin_edge_simp_dim = basic::get_file_dim(fin_edge_simp);
	auto node_inform = basic::read_file(fin_node_dim[0], fin_node_dim[1], fin_node);
	auto edge_orig_inform = basic::read_file(fin_edge_orig_dim[0], fin_edge_orig_dim[1], fin_edge_orig);
	auto edge_simp_inform = basic::read_file(fin_edge_simp_dim[0], fin_edge_simp_dim[1], fin_edge_simp);

	// Initialize node ID for DSO-Clusters
	int cluster_num = int(node_inform.col(1).maxCoeff());
	Power_network_inform.DSO_cluster.clear();
	Power_network_inform.DSO_cluster = std::vector <DSO_cluster> (cluster_num);
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
//	for(int node_iter = 0; node_iter < Power_network_inform.DSO_cluster[0].nodes_ID.size(); ++ node_iter){
//		std::cout << Power_network_inform.DSO_cluster[0].nodes_ID[node_iter] << " ";
//	}
//	std::cout << "\n\n";

	// Organize original edge data
	Power_network_inform.edges_orig.from = Eigen::VectorXi(fin_edge_orig_dim[0]);
	Power_network_inform.edges_orig.to = Eigen::VectorXi(fin_edge_orig_dim[0]);
	Power_network_inform.edges_orig.voltage_base = Eigen::VectorXi(fin_edge_orig_dim[0]);
	for(int edge_iter = 0; edge_iter < fin_edge_orig_dim[0]; ++ edge_iter){
		Power_network_inform.edges_orig.from(edge_iter) = int(edge_orig_inform(edge_iter, 0)) - 1;
		Power_network_inform.edges_orig.to(edge_iter) = int(edge_orig_inform(edge_iter, 1)) - 1;
		Power_network_inform.edges_orig.voltage_base(edge_iter) = int(edge_orig_inform(edge_iter, 4));
	}
	Power_network_inform.edges_orig.distance = edge_orig_inform.col(5);

	// Organize simplified egde data
	Power_network_inform.edges_simp.from = Eigen::VectorXi(fin_edge_simp_dim[0]);
	Power_network_inform.edges_simp.to = Eigen::VectorXi(fin_edge_simp_dim[0]);
	for(int edge_iter = 0; edge_iter < fin_edge_simp_dim[0]; ++ edge_iter){
		Power_network_inform.edges_simp.from(edge_iter) = int(edge_simp_inform(edge_iter, 0)) - 1;
		Power_network_inform.edges_simp.to(edge_iter) = int(edge_simp_inform(edge_iter, 1)) - 1;
	}
	Power_network_inform.edges_simp.conductance = edge_simp_inform.col(2);
}

void points_data_input(network_inform &Power_network_inform, Eigen::MatrixXd bz_inform, std::string fin_point, std::string fin_point_matrix){
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
		Power_network_inform.points.bidding_zone(point_iter) = bz_inform(int(point_inform(point_iter, 0)) - 1, 1) - 1;
		Power_network_inform.points.node(point_iter) = int(point_inform(point_iter, 1)) - 1;
		Power_network_inform.DSO_cluster[Power_network_inform.nodes.cluster(Power_network_inform.points.node(point_iter))].points_ID.push_back(point_iter);
		Power_network_inform.points.in_cluster_ID(point_iter) = Power_network_inform.DSO_cluster[Power_network_inform.nodes.cluster(Power_network_inform.points.node(point_iter))].points_ID.size() - 1;
	}
	Power_network_inform.points.population_density = point_inform.col(2);
	Power_network_inform.points.x = point_inform.col(point_inform.cols() - 4);
	Power_network_inform.points.y = point_inform.col(point_inform.cols() - 3);
	Power_network_inform.points.lon = point_inform.col(point_inform.cols() - 2);
	Power_network_inform.points.lat = point_inform.col(point_inform.cols() - 1);
//	for(int node_iter = 0; node_iter < Power_network_inform.DSO_cluster[0].nodes_ID.size(); ++ node_iter){
//		std::cout << Power_network_inform.DSO_cluster[0].nodes_ID[node_iter] << " ";
//	}
//	std::cout << "\n\n";
//	for(int point_iter = 0; point_iter < Power_network_inform.DSO_cluster[0].points_ID.size(); ++ point_iter){
//		std::cout << Power_network_inform.DSO_cluster[0].points_ID[point_iter] << " ";
//	}
//	std::cout << "\n\n";
//	std::cout << Power_network_inform.DSO_cluster.size() << "\n\n";

	// Read coordinate grid data
	Power_network_inform.points.coordinate_grid = Eigen::MatrixXi(fin_point_matrix_dim[0], fin_point_matrix_dim[1]);
	for(int x_coor_iter = 0; x_coor_iter < fin_point_matrix_dim[0]; ++ x_coor_iter){
		for(int y_coor_iter = 0; y_coor_iter < fin_point_matrix_dim[1]; ++ y_coor_iter){
			Power_network_inform.points.coordinate_grid(x_coor_iter, y_coor_iter) = std::max(int(point_matrix(x_coor_iter, y_coor_iter)), 0) - 1;
		}
	}

	// Calculate distance matrix
	geostat::point_distance_cov(Power_network_inform.points, 10.);
	//std::cout << Power_network_inform.points.distance.bottomRightCorner(5, 5) << "\n\n";
}

void plant_data_input(network_inform &Power_network_inform, std::string fin_hydro, std::string fin_wind){
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

void power_network_input_process(network_inform &Power_network_inform, std::string parent_dir){
	auto fin_bz = parent_dir + "DSO_Bidding_Zone.csv";
	auto fin_node = parent_dir + "transmission_nodes.csv";
	auto fin_edge_orig = parent_dir + "transmission_edges.csv";
	auto fin_edge_simp = parent_dir + "transmission_edges_pu_simp.csv";
	auto fin_point = parent_dir + "point_info.csv";
	auto fin_point_matrix = parent_dir + "point_matrix.csv";
	auto fin_hydro = parent_dir + "hydro_plants.csv";
	auto fin_wind = parent_dir + "wind_plants.csv";

	auto fin_bz_dim = basic::get_file_dim(fin_bz);
	auto bz_inform = basic::read_file(fin_bz_dim[0], fin_bz_dim[1], fin_bz);

	tranmission_data_input(Power_network_inform, bz_inform, fin_node, fin_edge_orig, fin_edge_simp);
	points_data_input(Power_network_inform, bz_inform, fin_point, fin_point_matrix);
	plant_data_input(Power_network_inform, fin_hydro, fin_wind);
}

//int main(){
//	network_inform Power_network_inform;
//	power_network_input_process(Power_network_inform);
//};
