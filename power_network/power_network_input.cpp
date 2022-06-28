// Source file for power network data input and process
#include <iostream>
#include "power_network.h"

void tranmission_data_input(network_inform &Power_network_inform, std::string fin_bz, std::string fin_node, std::string fin_edge_orig, std::string fin_edge_simp, std::string fin_pu){
	// Read power network data
	auto fin_bz_dim = get_file_dim(fin_bz);
	auto fin_node_dim = get_file_dim(fin_node);
	auto fin_edge_orig_dim = get_file_dim(fin_edge_orig);
	auto fin_edge_simp_dim = get_file_dim(fin_edge_simp);
	auto fin_pu_dim = get_file_dim(fin_pu);
	auto bz_inform = read_file(fin_bz_dim[0], fin_bz_dim[1], fin_bz);
	auto node_inform = read_file(fin_node_dim[0], fin_node_dim[1], fin_node);
	auto edge_orig_inform = read_file(fin_edge_orig_dim[0], fin_edge_orig_dim[1], fin_edge_orig);
	auto edge_simp_inform = read_file(fin_edge_simp_dim[0], fin_edge_simp_dim[1], fin_edge_simp);
	auto pu_inform = read_file(fin_pu_dim[0], fin_pu_dim[1], fin_pu);
	//std::cout << node_inform.topRows(10) << "\n\n";
	
	// Organize node cata
	Power_network_inform.nodes.bidding_zone = Eigen::VectorXi(fin_node_dim[0]);
	Power_network_inform.nodes.cluster = Eigen::VectorXi(fin_node_dim[0]);
	Power_network_inform.nodes.voltage_base = Eigen::VectorXi(fin_node_dim[0]);
	for(int node_iter = 0; node_iter < fin_node_dim[0]; ++ node_iter){
		Power_network_inform.nodes.bidding_zone(node_iter) = bz_inform(int(node_inform(node_iter, 0)) - 1, 1) - 1;
		Power_network_inform.nodes.cluster(node_iter) = int(node_inform(node_iter, 1)) - 1;
		Power_network_inform.nodes.voltage_base(node_iter) = int(node_inform(node_iter, 2));
	}
	Power_network_inform.nodes.x = node_inform.col(3);
	Power_network_inform.nodes.y = node_inform.col(4);
	Power_network_inform.nodes.lon = node_inform.col(5);
	Power_network_inform.nodes.lat = node_inform.col(6);
	
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
	
	std::cout << Power_network_inform.edges_orig.from.head(10).transpose() << "\n";
	std::cout << Power_network_inform.edges_orig.to.head(10).transpose() << "\n";
	std::cout << Power_network_inform.edges_orig.voltage_base.head(10).transpose() << "\n";
	std::cout << Power_network_inform.edges_orig.distance.head(10).transpose() << "\n";
}

int main(){
	auto fin_bz = "input/DSO_Bidding_Zone.csv";
	auto fin_node = "input/transmission_nodes.csv";
	auto fin_edge_orig = "input/transmission_edges.csv";
	auto fin_edge_simp = "input/transmission_edges_pu_simp.csv";
	auto fin_pu = "input/transmission_pu.csv";
	network_inform Power_network_inform;
	tranmission_data_input(Power_network_inform, fin_bz, fin_node, fin_edge_orig, fin_edge_simp, fin_pu);
	
}