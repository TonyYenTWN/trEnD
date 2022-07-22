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
		double z_base_high = pow(Power_network_inform.tech_parameters.voltage_cutoff_connection, 2.) / Power_network_inform.tech_parameters.s_base / 3.;
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
		admittance = admittance.array() / distance.array();

		// Set compact incidence matrix and edge admittance matrix
		double tol = 1.;
		DSO_Markets[DSO_iter].network.incidence.reserve(DSO_Markets[DSO_iter].network.num_vertice * DSO_Markets[DSO_iter].network.num_vertice);
		DSO_Markets[DSO_iter].network.admittance.reserve(DSO_Markets[DSO_iter].network.num_vertice * DSO_Markets[DSO_iter].network.num_vertice);
		std::vector <double> power_limit;
		power_limit.reserve(DSO_Markets[DSO_iter].network.num_vertice * DSO_Markets[DSO_iter].network.num_vertice);

		for(int row_iter = 0; row_iter < DSO_Markets[DSO_iter].network.num_vertice - 1; ++ row_iter){
			for(int col_iter = row_iter + 1; col_iter < DSO_Markets[DSO_iter].network.num_vertice; ++ col_iter){
				if(row_iter < Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() && col_iter < Power_network_inform.DSO_cluster[DSO_iter].points_ID.size()){
					if(admittance(row_iter, col_iter) > tol){
						DSO_Markets[DSO_iter].network.incidence.push_back(Eigen::Vector2i(row_iter, col_iter));
						DSO_Markets[DSO_iter].network.admittance.push_back(admittance(row_iter , col_iter));
						power_limit.push_back(Power_network_inform.tech_parameters.voltage_cutoff_distr);
					}
				}
				else{
					DSO_Markets[DSO_iter].network.incidence.push_back(Eigen::Vector2i(row_iter, col_iter));
					DSO_Markets[DSO_iter].network.admittance.push_back(admittance(row_iter , col_iter));
					power_limit.push_back(Power_network_inform.tech_parameters.voltage_cutoff_connection);
				}
			}
		}
		DSO_Markets[DSO_iter].network.num_edges = DSO_Markets[DSO_iter].network.incidence.size();

		// Set voltage and power constraints at each edge
		DSO_Markets[DSO_iter].network.voltage_constraint = Eigen::MatrixXd::Ones(DSO_Markets[DSO_iter].network.num_vertice, 2);
		DSO_Markets[DSO_iter].network.voltage_constraint.col(0) *= -Power_network_inform.tech_parameters.theta_limit;
		DSO_Markets[DSO_iter].network.voltage_constraint.col(1) *= Power_network_inform.tech_parameters.theta_limit;
		DSO_Markets[DSO_iter].network.power_constraint = Eigen::MatrixXd (DSO_Markets[DSO_iter].network.num_edges, 2);
		DSO_Markets[DSO_iter].network.power_constraint.col(1) = Eigen::Map <Eigen::VectorXd> (power_limit.data(), power_limit.size());
		DSO_Markets[DSO_iter].network.power_constraint.col(1) /= Power_network_inform.tech_parameters.s_base;
		DSO_Markets[DSO_iter].network.power_constraint.col(0) = -DSO_Markets[DSO_iter].network.power_constraint.col(1);

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

void power_market::Source_Node_Set(market_inform &DSO_Market, power_network::DSO_cluster &DSO_cluster){
	for(int node_iter = 0; node_iter < DSO_cluster.nodes_ID.size(); ++ node_iter){
		DSO_Market.submitted_demand.col(DSO_cluster.points_ID.size() + node_iter)	= Eigen::VectorXd::Zero(DSO_Market.price_intervals + 2);
		DSO_Market.submitted_supply(0, DSO_cluster.points_ID.size() + node_iter) =  std::numeric_limits<double>::infinity();
	}
}

void power_market::Sink_Node_Set(market_inform &DSO_Market, power_network::DSO_cluster &DSO_cluster){
	for(int node_iter = 0; node_iter < DSO_cluster.nodes_ID.size(); ++ node_iter){
		DSO_Market.submitted_supply.col(DSO_cluster.points_ID.size() + node_iter) = Eigen::VectorXd::Zero(DSO_Market.price_intervals + 2);
		DSO_Market.submitted_demand(DSO_Market.price_intervals + 1, DSO_cluster.points_ID.size() + node_iter) =  std::numeric_limits<double>::infinity();
	}
}

void power_market::DSO_Market_Results_Get(market_inform &Market, alglib::minlpstate &Problem, power_network::DSO_cluster &DSO_cluster, bool supply){
	alglib::real_1d_array sol;
	alglib::minlpreport rep;
	alglib::minlpresults(Problem, sol, rep);
	Eigen::VectorXd sol_vec = Eigen::Map <Eigen::VectorXd> (sol.getcontent(), sol.length());

	if(supply){
		// Store filtered supply bids at spatial points
		for(int point_iter = 0; point_iter < DSO_cluster.points_ID.size(); ++ point_iter){
			int row_start = 2 * Market.network.num_vertice + point_iter * (Market.price_intervals + 2);
			Market.filtered_supply.col(point_iter) = sol_vec.segment(row_start, Market.price_intervals + 2).array().max(0);
		}
	}
	else{
		// Store filtered demand bids at spatial points
		for(int point_iter = 0; point_iter < DSO_cluster.points_ID.size(); ++ point_iter){
			int row_start = 2 * Market.network.num_vertice + point_iter * (Market.price_intervals + 2);
			Market.filtered_demand.col(point_iter) = -(sol_vec.segment(row_start, Market.price_intervals + 2).array().min(0));
		}
	}

//	std::cout << sol_vec.segment(Market.network.num_vertice, DSO_cluster.points_ID.size()).minCoeff() << " " << sol_vec.segment(Market.network.num_vertice, DSO_cluster.points_ID.size()).maxCoeff() << " " << .5 * sol_vec.segment(Market.network.num_vertice, Market.network.num_vertice).array().abs().sum() << "\n";
//	std::cout << sol_vec.head(Market.network.num_vertice).minCoeff() << " " << sol_vec.head(Market.network.num_vertice).maxCoeff()  << "\n";
//	std::cout << sol_vec.tail(Market.network.num_edges).minCoeff() << " " << sol_vec.tail(Market.network.num_edges).maxCoeff() << "\n\n";
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
