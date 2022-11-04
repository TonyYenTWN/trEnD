// Source file for dispatch filtering of DSOs in Norway
// Note that in Norwegian terminology, HV distribution network (> 33kV) is called regional network.
#include <iostream>
#include "src/spatial_field/geostat.h"
#include "power_market.h"

void power_market::DSO_Markets_Set(markets_inform &DSO_Markets, power_network::network_inform &Power_network_inform, int Time){
	double pi = boost::math::constants::pi<double>();
	DSO_Markets.clear();
	DSO_Markets = std::vector <market_inform> (Power_network_inform.DSO_cluster.size());

	// Initialize markets at each clustered DSO
	for(int DSO_iter = 0; DSO_iter < DSO_Markets.size(); ++ DSO_iter){
		// Input parameters of DSO market
		int DSO_point_num = Power_network_inform.DSO_cluster[DSO_iter].points_ID.size();
		int DSO_node_num = Power_network_inform.DSO_cluster[DSO_iter].nodes_ID.size();
		if(DSO_point_num == 0){
			DSO_Markets[DSO_iter].num_zone = DSO_node_num;
			DSO_Markets[DSO_iter].network.num_edges = 0;
			//std::cout << DSO_iter << ":\t" << DSO_Markets[DSO_iter].num_zone << "\t" << DSO_Markets[DSO_iter].network.num_edges << "\n";
			continue;
		}

		DSO_Markets[DSO_iter].num_zone = DSO_point_num + DSO_node_num;
		DSO_Markets[DSO_iter].time_intervals = Time;
		parameters::bidded_price(DSO_Markets[DSO_iter].bidded_price_map);

		// Set compact incidence matrix and edge admittance matrix
		double tol = 1.;
		double power_limit_connection = 1.;
		double power_limit_distr = .5;
		DSO_Markets[DSO_iter].network.incidence.reserve(DSO_Markets[DSO_iter].network.num_vertice * DSO_Markets[DSO_iter].network.num_vertice);
		DSO_Markets[DSO_iter].network.admittance.reserve(DSO_Markets[DSO_iter].network.num_vertice * DSO_Markets[DSO_iter].network.num_vertice);
		std::vector <double> power_limit;
		power_limit.reserve(DSO_Markets[DSO_iter].network.num_vertice * DSO_Markets[DSO_iter].network.num_vertice);

		// Calculate edge admittance
		// Use fractional Laplacian here :
		// y(x_1, x_2) = N / 4 / pi / z / |x_1 - x_2|^(2 + d), where
		// y(x_1, x_2): admittance between x_1 and x_2
		// N: mean line density at the neighborhood of x_1 containing also x_2
		// z: per length impedence of power lines
		// d: fractional dimension of the power line distribution
		double z_base_low = pow(Power_network_inform.tech_parameters.voltage_cutoff_distr, 2.) / Power_network_inform.tech_parameters.s_base * 3.;
		double z_base_high = pow(Power_network_inform.tech_parameters.voltage_cutoff_connection, 2.) / Power_network_inform.tech_parameters.s_base * 3.;
		double partition_func = 0.;
		DSO_Markets[DSO_iter].network.num_vertice = DSO_Markets[DSO_iter].num_zone;
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
				// Series admittance
				double y_series = 1.;
				y_series /= distance(row_iter, col_iter);
				y_series /= Power_network_inform.tech_parameters.z_distr_series.imag();
				y_series *= z_base_low;
				y_series *= num_line(row_iter, col_iter);

				DSO_Markets[DSO_iter].network.incidence.push_back(Eigen::Vector2i(row_iter, col_iter));
				DSO_Markets[DSO_iter].network.admittance.push_back(y_series);
				power_limit.push_back(power_limit_distr * Power_network_inform.tech_parameters.voltage_cutoff_distr * num_line(row_iter, col_iter));
			}
		}

		// Connection between nodes and points
		for(int row_iter = 0; row_iter < DSO_node_num; ++ row_iter){
			int node_ID = Power_network_inform.DSO_cluster[DSO_iter].nodes_ID[row_iter];
			int min_point_ID;
			double distance_min = std::numeric_limits<double>::infinity();
			double line_density_connection = Power_network_inform.tech_parameters.line_density_connection * DSO_point_num / DSO_node_num;

			for(int col_iter = 0; col_iter < DSO_point_num ; ++ col_iter){
				int point_ID = Power_network_inform.DSO_cluster[DSO_iter].points_ID[col_iter];
				Eigen::Vector2d point_coor = Eigen::Vector2d(Power_network_inform.points.lon(point_ID), Power_network_inform.points.lat(point_ID));
				Eigen::Vector2d node_coor = Eigen::Vector2d(Power_network_inform.nodes.lon(node_ID), Power_network_inform.nodes.lat(node_ID));
				point_coor *= pi / 180.;
				node_coor *= pi / 180.;
				double distance_temp = spatial_field::geodist(point_coor, node_coor);

				if(distance_temp < distance_min){
					distance_min = distance_temp;
					min_point_ID = col_iter;
				}
			}

			// Series admittance
			double y_series = 1.;
			y_series /= distance_min;
			y_series /= Power_network_inform.tech_parameters.z_distr_series.imag();
			y_series *= z_base_high;
			y_series *= line_density_connection;

			DSO_Markets[DSO_iter].network.incidence.push_back(Eigen::Vector2i(min_point_ID, DSO_point_num + row_iter));
			DSO_Markets[DSO_iter].network.admittance.push_back(y_series);
			power_limit.push_back(power_limit_connection * line_density_connection);
		}
		DSO_Markets[DSO_iter].network.num_edges = DSO_Markets[DSO_iter].network.incidence.size();
		//std::cout << DSO_iter << ":\t" << DSO_Markets[DSO_iter].num_zone << "\t" << DSO_Markets[DSO_iter].network.num_edges << "\n";

//		Eigen::MatrixXd admittance = Eigen::MatrixXd::Ones(DSO_Markets[DSO_iter].network.num_vertice, DSO_Markets[DSO_iter].network.num_vertice);
//		Eigen::MatrixXd num_line = Eigen::MatrixXd::Ones(DSO_Markets[DSO_iter].network.num_vertice, DSO_Markets[DSO_iter].network.num_vertice);
//		Eigen::MatrixXd distance = Eigen::MatrixXd::Zero(DSO_Markets[DSO_iter].network.num_vertice, DSO_Markets[DSO_iter].network.num_vertice);
//
////		std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
////		std::cout << DSO_iter << "\n";
////		std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
//		for(int row_iter = 0; row_iter < DSO_Markets[DSO_iter].network.num_vertice - 1; ++ row_iter){
//			for(int col_iter = row_iter + 1; col_iter < DSO_Markets[DSO_iter].network.num_vertice; ++ col_iter){
//				if(row_iter < point_num && col_iter < point_num){
//					int point_ID_1 = Power_network_inform.DSO_cluster[DSO_iter].points_ID[row_iter];
//					int point_ID_2 = Power_network_inform.DSO_cluster[DSO_iter].points_ID[col_iter];
//					distance(row_iter, col_iter) = Power_network_inform.points.distance(point_ID_1, point_ID_2);
//					admittance(row_iter, col_iter) = 1. / pow(distance(row_iter, col_iter) * 1E-3, 1. + Power_network_inform.tech_parameters.fraction_dim_distr);
//					partition_func += admittance(row_iter, col_iter);
//				}
//				else{
//					if(row_iter < point_num){
//						int point_ID = Power_network_inform.DSO_cluster[DSO_iter].points_ID[row_iter];
//						int node_ID = Power_network_inform.DSO_cluster[DSO_iter].nodes_ID[col_iter - Power_network_inform.DSO_cluster[DSO_iter].points_ID.size()];
//						Eigen::Vector2d point_coor = Eigen::Vector2d(Power_network_inform.points.lon(point_ID), Power_network_inform.points.lat(point_ID));
//						Eigen::Vector2d node_coor = Eigen::Vector2d(Power_network_inform.nodes.lon(node_ID), Power_network_inform.nodes.lat(node_ID));
//						point_coor *= pi / 180.;
//						node_coor *= pi / 180.;
//						distance(row_iter, col_iter) = spatial_field::geodist(point_coor, node_coor);
//					}
//					else{
//						int node_ID_1 = Power_network_inform.DSO_cluster[DSO_iter].nodes_ID[row_iter - Power_network_inform.DSO_cluster[DSO_iter].points_ID.size()];
//						int node_ID_2 = Power_network_inform.DSO_cluster[DSO_iter].nodes_ID[col_iter - Power_network_inform.DSO_cluster[DSO_iter].points_ID.size()];
//						Eigen::Vector2d node_coor_1 = Eigen::Vector2d(Power_network_inform.nodes.lon(node_ID_1), Power_network_inform.nodes.lat(node_ID_1));
//						Eigen::Vector2d node_coor_2 = Eigen::Vector2d(Power_network_inform.nodes.lon(node_ID_2), Power_network_inform.nodes.lat(node_ID_2));
//						node_coor_1 *= pi / 180.;
//						node_coor_2 *= pi / 180.;
//						distance(row_iter, col_iter) = spatial_field::geodist(node_coor_1, node_coor_2);
//					}
//				}
//			}
//		}
//		admittance.topLeftCorner(point_num, point_num) *= Power_network_inform.tech_parameters.line_density_distr * point_num;
//		num_line.topLeftCorner(point_num, point_num) = admittance.topLeftCorner(point_num, point_num) / partition_func;
//		admittance.topLeftCorner(point_num, point_num) /= partition_func * Power_network_inform.tech_parameters.z_distr_series.imag() / z_base_low;
//		admittance.rightCols(Power_network_inform.DSO_cluster[DSO_iter].nodes_ID.size()) *= z_base_high / Power_network_inform.tech_parameters.z_distr_series.imag();
//		num_line.rightCols(Power_network_inform.DSO_cluster[DSO_iter].nodes_ID.size()) *= Power_network_inform.tech_parameters.line_density_connection * point_num;
//		admittance = admittance.array() / distance.array();

//		// Set compact incidence matrix and edge admittance matrix
//		double tol = 1.;
//		double power_limit_connection = 1.;
//		double power_limit_distr = .5;
//		DSO_Markets[DSO_iter].network.incidence.reserve(DSO_Markets[DSO_iter].network.num_vertice * DSO_Markets[DSO_iter].network.num_vertice);
//		DSO_Markets[DSO_iter].network.admittance.reserve(DSO_Markets[DSO_iter].network.num_vertice * DSO_Markets[DSO_iter].network.num_vertice);
//		std::vector <double> power_limit;
//		power_limit.reserve(DSO_Markets[DSO_iter].network.num_vertice * DSO_Markets[DSO_iter].network.num_vertice);
//
//		for(int row_iter = 0; row_iter < DSO_Markets[DSO_iter].network.num_vertice - 1; ++ row_iter){
//			for(int col_iter = row_iter + 1; col_iter < DSO_Markets[DSO_iter].network.num_vertice; ++ col_iter){
//				if(row_iter < point_num && col_iter < point_num){
//					if(admittance(row_iter, col_iter) > tol){
//						DSO_Markets[DSO_iter].network.incidence.push_back(Eigen::Vector2i(row_iter, col_iter));
//						DSO_Markets[DSO_iter].network.admittance.push_back(admittance(row_iter , col_iter));
//						power_limit.push_back(power_limit_distr * Power_network_inform.tech_parameters.voltage_cutoff_distr * num_line(row_iter, col_iter));
//					}
//				}
//				else{
//					DSO_Markets[DSO_iter].network.incidence.push_back(Eigen::Vector2i(row_iter, col_iter));
//					DSO_Markets[DSO_iter].network.admittance.push_back(admittance(row_iter, col_iter));
//					power_limit.push_back(power_limit_connection * Power_network_inform.tech_parameters.voltage_cutoff_connection * num_line(row_iter, col_iter));
//				}
//			}
//		}
//		DSO_Markets[DSO_iter].network.num_edges = DSO_Markets[DSO_iter].network.incidence.size();

		// Set voltage and power constraints at each edge
		DSO_Markets[DSO_iter].network.voltage_constraint = Eigen::MatrixXd::Ones(DSO_Markets[DSO_iter].network.num_vertice, 2);
		DSO_Markets[DSO_iter].network.voltage_constraint.col(0) *= -Power_network_inform.tech_parameters.theta_distr_limit;
		DSO_Markets[DSO_iter].network.voltage_constraint.col(1) *= Power_network_inform.tech_parameters.theta_distr_limit;
		DSO_Markets[DSO_iter].network.power_constraint = Eigen::MatrixXd (DSO_Markets[DSO_iter].network.num_edges, 2);
		DSO_Markets[DSO_iter].network.power_constraint.col(1) = Eigen::Map <Eigen::VectorXd> (power_limit.data(), power_limit.size());
		DSO_Markets[DSO_iter].network.power_constraint.col(1) /= Power_network_inform.tech_parameters.s_base;
		DSO_Markets[DSO_iter].network.power_constraint.col(0) = -DSO_Markets[DSO_iter].network.power_constraint.col(1);

		// Initialization of process variables
		power_market::Market_Initialization(DSO_Markets[DSO_iter]);

		// Initialization of output variables
		DSO_Markets[DSO_iter].confirmed.supply = Eigen::MatrixXd::Zero(Time, DSO_Markets[DSO_iter].num_zone);
		DSO_Markets[DSO_iter].confirmed.demand = Eigen::MatrixXd::Zero(Time, DSO_Markets[DSO_iter].num_zone);
		DSO_Markets[DSO_iter].confirmed.price = Eigen::MatrixXd::Zero(Time, DSO_Markets[DSO_iter].num_zone);
		DSO_Markets[DSO_iter].confirmed.ratio_supply = Eigen::VectorXd::Zero(DSO_Markets[DSO_iter].num_zone);
		DSO_Markets[DSO_iter].confirmed.ratio_demand = Eigen::VectorXd::Zero(DSO_Markets[DSO_iter].num_zone);
	}
}

namespace{
	void Source_Node_Set(power_market::market_inform &DSO_Market, power_network::DSO_cluster &DSO_cluster){
		for(int node_iter = 0; node_iter < DSO_cluster.nodes_ID.size(); ++ node_iter){
			DSO_Market.submitted_demand.col(DSO_cluster.points_ID.size() + node_iter)	= Eigen::VectorXd::Zero(DSO_Market.price_intervals + 2);
			DSO_Market.submitted_supply(0, DSO_cluster.points_ID.size() + node_iter) =  std::numeric_limits<double>::infinity();
		}
	}

	void Sink_Node_Set(power_market::market_inform &DSO_Market, power_network::DSO_cluster &DSO_cluster){
		for(int node_iter = 0; node_iter < DSO_cluster.nodes_ID.size(); ++ node_iter){
			DSO_Market.submitted_supply.col(DSO_cluster.points_ID.size() + node_iter) = Eigen::VectorXd::Zero(DSO_Market.price_intervals + 2);
			DSO_Market.submitted_demand(DSO_Market.price_intervals + 1, DSO_cluster.points_ID.size() + node_iter) =  std::numeric_limits<double>::infinity();
		}
	}
}

void power_market::Filtered_bid_demand_calculation(int tick, market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
	for(int DSO_iter = 0; DSO_iter < Power_market_inform.DSO_Markets.size(); ++ DSO_iter){
		if(Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() == 0){
			continue;
		}

		// Initialize submit bids of DSO markets
		Market_Initialization(Power_market_inform.DSO_Markets[DSO_iter]);

		// Create source nodes
		Source_Node_Set(Power_market_inform.DSO_Markets[DSO_iter], Power_network_inform.DSO_cluster[DSO_iter]);
	}

	// Residential demand at each point
	int point_num = Power_network_inform.points.bidding_zone.size();
	int sample_num = Power_market_inform.agent_profiles.end_users[0].size();
	for(int point_iter = 0; point_iter < point_num; ++ point_iter){
		int point_ID = Power_network_inform.points.in_cluster_ID(point_iter);
		int node_ID = Power_network_inform.points.node(point_iter);
		int DSO_ID = Power_network_inform.nodes.cluster(node_ID);

		for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
			Power_market_inform.DSO_Markets[DSO_ID].submitted_demand.col(point_ID) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.filter_demand;
		}
	}

	// LV Power suppliers
	int pump_LV_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
	for(int agent_iter = 0; agent_iter < pump_LV_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);
		int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
		point_ID = Power_network_inform.points.in_cluster_ID(point_ID);

		Power_market_inform.DSO_Markets[DSO_ID].submitted_demand.col(point_ID) += Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.filter_demand;
	}

	for(int DSO_iter = 0; DSO_iter < Power_market_inform.DSO_Markets.size(); ++ DSO_iter){
		if(Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() == 0){
			continue;
		}

		Flow_Based_Market_Optimization(Power_market_inform.DSO_Markets[DSO_iter], Power_market_inform.DSO_Problems[DSO_iter]);
		DSO_Market_Results_Get(tick, Power_market_inform.DSO_Markets[DSO_iter], Power_market_inform.DSO_Problems[DSO_iter], Power_network_inform.DSO_cluster[DSO_iter], 0);
	}
}

void power_market::Filtered_bid_supply_calculation(int tick, market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
	for(int DSO_iter = 0; DSO_iter < Power_market_inform.DSO_Markets.size(); ++ DSO_iter){
		if(Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() == 0){
			continue;
		}

		// Initialize submit bids of DSO markets
		Market_Initialization(Power_market_inform.DSO_Markets[DSO_iter]);

		// Create source nodes
		Sink_Node_Set(Power_market_inform.DSO_Markets[DSO_iter], Power_network_inform.DSO_cluster[DSO_iter]);
	}

	// Residential demand at each point
	int point_num = Power_network_inform.points.bidding_zone.size();
	int sample_num = Power_market_inform.agent_profiles.end_users[0].size();
	for(int point_iter = 0; point_iter < point_num; ++ point_iter){
		int point_ID = Power_network_inform.points.in_cluster_ID(point_iter);
		int node_ID = Power_network_inform.points.node(point_iter);
		int DSO_ID = Power_network_inform.nodes.cluster(node_ID);

		for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
			Power_market_inform.DSO_Markets[DSO_ID].submitted_supply.col(point_ID) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.filter_supply;
		}
	}

	// LV Power suppliers
	int hydro_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);
		int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
		point_ID = Power_network_inform.points.in_cluster_ID(point_ID);

		Power_market_inform.DSO_Markets[DSO_ID].submitted_supply.col(point_ID) += Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids.filter_supply;
	}

	int wind_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);
		int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
		point_ID = Power_network_inform.points.in_cluster_ID(point_ID);

		Power_market_inform.DSO_Markets[DSO_ID].submitted_supply.col(point_ID) += Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.filter_supply;
	}

	int pump_LV_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
	for(int agent_iter = 0; agent_iter < pump_LV_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);
		int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
		point_ID = Power_network_inform.points.in_cluster_ID(point_ID);

		Power_market_inform.DSO_Markets[DSO_ID].submitted_supply.col(point_ID) += Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.filter_supply;
	}

	for(int DSO_iter = 0; DSO_iter < Power_market_inform.DSO_Markets.size(); ++ DSO_iter){
		if(Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() == 0){
			continue;
		}

		Flow_Based_Market_Optimization(Power_market_inform.DSO_Markets[DSO_iter], Power_market_inform.DSO_Problems[DSO_iter]);
		DSO_Market_Results_Get(tick, Power_market_inform.DSO_Markets[DSO_iter], Power_market_inform.DSO_Problems[DSO_iter], Power_network_inform.DSO_cluster[DSO_iter], 1);
	}
}

void power_market::DSO_Market_Results_Get(int tick, market_inform &Market, alglib::minlpstate &Problem, power_network::DSO_cluster &DSO_cluster, bool supply){
	alglib::real_1d_array sol;
	alglib::minlpreport rep;
	alglib::minlpresults(Problem, sol, rep);
	Eigen::VectorXd sol_vec = Eigen::Map <Eigen::VectorXd> (sol.getcontent(), sol.length());

	if(supply){
		// Store filtered supply bids at spatial points
		for(int point_iter = 0; point_iter < DSO_cluster.points_ID.size(); ++ point_iter){
			// Store power sink
			int row_start = 2 * Market.network.num_vertice + point_iter * (Market.price_intervals + 2);
			Market.confirmed.supply(tick, point_iter) = (sol_vec.segment(row_start, Market.price_intervals + 2).array().max(0)).sum();

			// Store nodal prices
			Market.confirmed.price(tick, point_iter) = Market.bidded_price_map.bidded_price(1) + rep.lagbc[row_start + 1];
			Market.confirmed.price(tick, point_iter) = int(Market.confirmed.price(tick, point_iter)) + .5;
			if(Market.confirmed.price(tick, point_iter) < Market.bidded_price_map.bidded_price(1)){
				Market.confirmed.price(tick, point_iter) = Market.bidded_price_map.bidded_price(0);
			}
			else if(Market.confirmed.price(tick, point_iter) > Market.bidded_price_map.bidded_price(Market.price_intervals)){
				Market.confirmed.price(tick, point_iter) = Market.bidded_price_map.bidded_price(Market.price_intervals + 1);
			}

			// Store ratio at nodes
			for(int price_iter = 0; price_iter < Market.price_intervals + 2; ++ price_iter){
				if(Market.bidded_price_map.bidded_price(price_iter) >= Market.confirmed.price(tick, point_iter) || price_iter == Market.price_intervals + 1){
					if(sol[row_start + price_iter] >= 0.){
						Market.confirmed.ratio_demand(point_iter) = std::min(Market.submitted_demand(price_iter, point_iter), Market.submitted_supply(price_iter, point_iter) - sol[row_start + price_iter]);
						Market.confirmed.ratio_supply(point_iter) = Market.confirmed.ratio_demand(point_iter) + sol[row_start + price_iter];
						Market.confirmed.ratio_supply(point_iter) /= Market.submitted_supply(price_iter, point_iter) + 1E-12;
					}
					else{
						Market.confirmed.ratio_supply(point_iter) = std::min(Market.submitted_supply(price_iter, point_iter), Market.submitted_demand(price_iter, point_iter) + sol[row_start + price_iter]);
						Market.confirmed.ratio_supply(point_iter) /= Market.submitted_supply(price_iter, point_iter) + 1E-12;
					}
					break;
				}
			}
		}
	}
	else{
		// Store filtered demand bids at spatial points
		for(int point_iter = 0; point_iter < DSO_cluster.points_ID.size(); ++ point_iter){
			// Store power sink
			int row_start = 2 * Market.network.num_vertice + point_iter * (Market.price_intervals + 2);
			Market.confirmed.demand(tick, point_iter) = -(sol_vec.segment(row_start, Market.price_intervals + 2).array().min(0)).sum();

			// Store nodal prices
			Market.confirmed.price(tick, point_iter) = Market.bidded_price_map.bidded_price(0) + rep.lagbc[row_start];
			Market.confirmed.price(tick, point_iter) = int(Market.confirmed.price(tick, point_iter)) + .5;
			if(Market.confirmed.price(tick, point_iter) < Market.bidded_price_map.bidded_price(1)){
				Market.confirmed.price(tick, point_iter) = Market.bidded_price_map.bidded_price(0);
			}
			else if(Market.confirmed.price(tick, point_iter) > Market.bidded_price_map.bidded_price(Market.price_intervals)){
				Market.confirmed.price(tick, point_iter) = Market.bidded_price_map.bidded_price(Market.price_intervals + 1);
			}

			// Store ratio at nodes
			for(int price_iter = 0; price_iter < Market.price_intervals + 2; ++ price_iter){
				if(Market.bidded_price_map.bidded_price(price_iter) >= Market.confirmed.price(tick, point_iter) || price_iter == Market.price_intervals + 1){
					if(sol[row_start + price_iter] >= 0.){
						Market.confirmed.ratio_demand(point_iter) = std::min(Market.submitted_demand(price_iter, point_iter), Market.submitted_supply(price_iter, point_iter) - sol[row_start + price_iter]);
						Market.confirmed.ratio_demand(point_iter) /= Market.submitted_demand(price_iter, point_iter) + 1E-12;
					}
					else{
						Market.confirmed.ratio_supply(point_iter) = std::min(Market.submitted_supply(price_iter, point_iter), Market.submitted_demand(price_iter, point_iter) + sol[row_start + price_iter]);
						Market.confirmed.ratio_demand(point_iter) = Market.confirmed.ratio_supply(point_iter) - sol[row_start + price_iter];
						Market.confirmed.ratio_demand(point_iter) /= Market.submitted_demand(price_iter, point_iter) + 1E-12;
					}
					break;
				}
			}
		}
	}

//	std::cout << sol_vec.segment(Market.network.num_vertice, DSO_cluster.points_ID.size()).minCoeff() << " " << sol_vec.segment(Market.network.num_vertice, DSO_cluster.points_ID.size()).maxCoeff() << " " << .5 * sol_vec.segment(Market.network.num_vertice, Market.network.num_vertice).array().abs().sum() << "\n";
//	std::cout << sol_vec.head(Market.network.num_vertice).minCoeff() << " " << sol_vec.head(Market.network.num_vertice).maxCoeff()  << "\n";
//	std::cout << sol_vec.tail(Market.network.num_edges).minCoeff() << " " << sol_vec.tail(Market.network.num_edges).maxCoeff() << "\n\n";
}
