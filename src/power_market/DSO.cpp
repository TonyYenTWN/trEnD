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
		DSO_Markets[DSO_iter].num_zone = Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() + Power_network_inform.DSO_cluster[DSO_iter].nodes_ID.size();
		DSO_Markets[DSO_iter].time_intervals = Time;
		parameters::bidded_price(DSO_Markets[DSO_iter].bidded_price_map);
		int point_num = Power_network_inform.DSO_cluster[DSO_iter].points_ID.size();

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
		Eigen::MatrixXd num_line = Eigen::MatrixXd::Ones(DSO_Markets[DSO_iter].network.num_vertice, DSO_Markets[DSO_iter].network.num_vertice);
		Eigen::MatrixXd distance = Eigen::MatrixXd::Zero(DSO_Markets[DSO_iter].network.num_vertice, DSO_Markets[DSO_iter].network.num_vertice);

//		std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
//		std::cout << DSO_iter << "\n";
//		std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------------\n";
		for(int row_iter = 0; row_iter < DSO_Markets[DSO_iter].network.num_vertice - 1; ++ row_iter){
			for(int col_iter = row_iter + 1; col_iter < DSO_Markets[DSO_iter].network.num_vertice; ++ col_iter){
				if(row_iter < point_num && col_iter < point_num){
					int point_ID_1 = Power_network_inform.DSO_cluster[DSO_iter].points_ID[row_iter];
					int point_ID_2 = Power_network_inform.DSO_cluster[DSO_iter].points_ID[col_iter];
					distance(row_iter, col_iter) = Power_network_inform.points.distance(point_ID_1, point_ID_2);
					admittance(row_iter, col_iter) = 1. / pow(distance(row_iter, col_iter) * 1E-3, 1. + Power_network_inform.tech_parameters.fraction_dim_distr);
					partition_func += admittance(row_iter, col_iter);
				}
				else{
					if(row_iter < point_num){
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
		admittance.topLeftCorner(point_num, point_num) *= Power_network_inform.tech_parameters.line_density_distr * point_num;
		num_line.topLeftCorner(point_num, point_num) = admittance.topLeftCorner(point_num, point_num) / partition_func;
		admittance.topLeftCorner(point_num, point_num) /= partition_func * Power_network_inform.tech_parameters.z_distr_series.imag() / z_base_low;
		admittance.rightCols(Power_network_inform.DSO_cluster[DSO_iter].nodes_ID.size()) *= z_base_high / Power_network_inform.tech_parameters.z_distr_series.imag();
		num_line.rightCols(Power_network_inform.DSO_cluster[DSO_iter].nodes_ID.size()) *= Power_network_inform.tech_parameters.line_density_connection * point_num;
		admittance = admittance.array() / distance.array();

		// Set compact incidence matrix and edge admittance matrix
		double tol = 1.;
		DSO_Markets[DSO_iter].network.incidence.reserve(DSO_Markets[DSO_iter].network.num_vertice * DSO_Markets[DSO_iter].network.num_vertice);
		DSO_Markets[DSO_iter].network.admittance.reserve(DSO_Markets[DSO_iter].network.num_vertice * DSO_Markets[DSO_iter].network.num_vertice);
		std::vector <double> power_limit;
		power_limit.reserve(DSO_Markets[DSO_iter].network.num_vertice * DSO_Markets[DSO_iter].network.num_vertice);

		for(int row_iter = 0; row_iter < DSO_Markets[DSO_iter].network.num_vertice - 1; ++ row_iter){
			for(int col_iter = row_iter + 1; col_iter < DSO_Markets[DSO_iter].network.num_vertice; ++ col_iter){
				if(row_iter < point_num && col_iter < point_num){
					if(admittance(row_iter, col_iter) > tol){
						DSO_Markets[DSO_iter].network.incidence.push_back(Eigen::Vector2i(row_iter, col_iter));
						DSO_Markets[DSO_iter].network.admittance.push_back(admittance(row_iter , col_iter));
						power_limit.push_back(Power_network_inform.tech_parameters.voltage_cutoff_distr * num_line(row_iter, col_iter));
					}
				}
				else{
					DSO_Markets[DSO_iter].network.incidence.push_back(Eigen::Vector2i(row_iter, col_iter));
					DSO_Markets[DSO_iter].network.admittance.push_back(admittance(row_iter, col_iter));
					power_limit.push_back(Power_network_inform.tech_parameters.voltage_cutoff_connection * num_line(row_iter, col_iter));
				}
			}
		}
		DSO_Markets[DSO_iter].network.num_edges = DSO_Markets[DSO_iter].network.incidence.size();

		// Set voltage and power constraints at each edge
		DSO_Markets[DSO_iter].network.voltage_constraint = Eigen::MatrixXd::Ones(DSO_Markets[DSO_iter].network.num_vertice, 2);
		DSO_Markets[DSO_iter].network.voltage_constraint.col(0) *= -Power_network_inform.tech_parameters.theta_limit;
		DSO_Markets[DSO_iter].network.voltage_constraint.col(1) *= Power_network_inform.tech_parameters.theta_limit ;
		DSO_Markets[DSO_iter].network.power_constraint = Eigen::MatrixXd (DSO_Markets[DSO_iter].network.num_edges, 2);
		DSO_Markets[DSO_iter].network.power_constraint.col(1) = Eigen::Map <Eigen::VectorXd> (power_limit.data(), power_limit.size());
		DSO_Markets[DSO_iter].network.power_constraint.col(1) /= Power_network_inform.tech_parameters.s_base;
		DSO_Markets[DSO_iter].network.power_constraint.col(0) = -DSO_Markets[DSO_iter].network.power_constraint.col(1);

		// Initialization of process variables
		power_market::Market_Initialization(DSO_Markets[DSO_iter]);

		// Initialization of output variables
		DSO_Markets[DSO_iter].confirmed_supply = Eigen::MatrixXd::Zero(Time, DSO_Markets[DSO_iter].num_zone);
		DSO_Markets[DSO_iter].confirmed_demand = Eigen::MatrixXd::Zero(Time, DSO_Markets[DSO_iter].num_zone);
		DSO_Markets[DSO_iter].confirmed_price = Eigen::MatrixXd::Zero(Time, DSO_Markets[DSO_iter].num_zone);
		DSO_Markets[DSO_iter].confirmed_ratio_supply = Eigen::VectorXd::Zero(DSO_Markets[DSO_iter].num_zone);
		DSO_Markets[DSO_iter].confirmed_ratio_demand = Eigen::VectorXd::Zero(DSO_Markets[DSO_iter].num_zone);
	}
}

//void power_market::DSO_agents_update(int tick, agent::end_user::profiles &end_user_profiles, market_inform &TSO_Market, market_inform &International_Market, power_network::network_inform &Power_network_inform){
////	// Update state of agents from the previous tick
////	for(int point_iter = 0; point_iter < end_user_profiles.size(); ++ point_iter){
////		int node_ID = Power_network_inform.points.node(point_iter);
////		int sample_num = agent::parameters::sample_num();
////		double node_price = TSO_Market.actual_price(tick - 1, node_ID);
////
////		for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
////			end_user_profiles[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand = 0.;
////			end_user_profiles[point_iter][sample_iter].operation.EV.BESS.soc_ini -= end_user_profiles[point_iter][sample_iter].operation.EV.BESS.normalized_scheduled_capacity_profile(0);
////			end_user_profiles[point_iter][sample_iter].operation.BESS.soc_ini -= end_user_profiles[point_iter][sample_iter].operation.BESS.normalized_scheduled_capacity_profile(0);
////
////			// Demand settlement
////			// Order of reduction: BESS charge, EV charge, smart appliance, inflexible demand
////			double marginal_demand = 0.;
////			marginal_demand += (node_price == end_user_profiles[point_iter][sample_iter].operation.demand_flex_price) * (end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile(0) >= 0.) * end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile(0);
////			marginal_demand += (node_price == end_user_profiles[point_iter][sample_iter].operation.demand_inflex_price) * (end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile(0) >= 0.) * end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile(0);
////			double demand_gap = (1. - TSO_Market.actual_price_ratio(tick - 1, node_ID)) * marginal_demand;
////			demand_gap += (node_price > end_user_profiles[point_iter][sample_iter].operation.demand_flex_price) * (end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile(0) >= 0.) * end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile(0);
////			demand_gap += (node_price > end_user_profiles[point_iter][sample_iter].operation.demand_inflex_price) * (end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile(0) >= 0.) * end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile(0);
////
////			if(demand_gap > 0.){
////				// Reduce BESS charge
////				double reduction = std::min(demand_gap, end_user_profiles[point_iter][sample_iter].operation.BESS.capacity_scale - end_user_profiles[point_iter][sample_iter].operation.BESS.normalized_scheduled_capacity_profile(0));
////				end_user_profiles[point_iter][sample_iter].operation.BESS.soc_ini -= reduction;
////				end_user_profiles[point_iter][sample_iter].operation.BESS.normalized_scheduled_capacity_profile(0) += reduction;
////				demand_gap -= reduction;
////
////				if(demand_gap > 0.){
////					// Reduce EV charge
////					reduction = std::min(demand_gap, end_user_profiles[point_iter][sample_iter].operation.EV.BESS.capacity_scale - end_user_profiles[point_iter][sample_iter].operation.EV.BESS.normalized_scheduled_capacity_profile(0));
////					end_user_profiles[point_iter][sample_iter].operation.EV.BESS.soc_ini -= reduction;
////					end_user_profiles[point_iter][sample_iter].operation.EV.BESS.normalized_scheduled_capacity_profile(0) += reduction;
////					demand_gap -= reduction;
////
////					if(demand_gap > 0.){
////						// Reduce smart appliance
////						reduction = std::min(demand_gap, end_user_profiles[point_iter][sample_iter].operation.smart_appliance.normalized_scheduled_profile(0));
////						end_user_profiles[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand += reduction;
////						demand_gap -= reduction;
////					}
////				}
////			}
////
////			// Supply settlement
////			// Order of reduction: BESS discharge, EV discharge, inflexible supply
////			double marginal_supply = 0.;
////			marginal_supply += (node_price == end_user_profiles[point_iter][sample_iter].operation.supply_flex_price) * (end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile(0) < 0.) * end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile(0);
////			marginal_supply += (node_price == end_user_profiles[point_iter][sample_iter].operation.supply_inflex_price) * (end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile(0) < 0.) * end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile(0);
////			double supply_gap = (1. - TSO_Market.actual_price_ratio(tick - 1, node_ID)) * marginal_supply;
////			supply_gap += (node_price > end_user_profiles[point_iter][sample_iter].operation.supply_flex_price) * (end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile(0) < 0.) * end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile(0);
////			supply_gap += (node_price > end_user_profiles[point_iter][sample_iter].operation.supply_inflex_price) * (end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile(0) < 0.) * end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile(0);
////
////			if(supply_gap > 0.){
////				// Reduce BESS charge
////				double reduction = std::min(supply_gap, end_user_profiles[point_iter][sample_iter].operation.BESS.capacity_scale + end_user_profiles[point_iter][sample_iter].operation.BESS.normalized_scheduled_capacity_profile(0));
////				end_user_profiles[point_iter][sample_iter].operation.BESS.soc_ini += reduction;
////				end_user_profiles[point_iter][sample_iter].operation.BESS.normalized_scheduled_capacity_profile(0) -= reduction;
////				supply_gap -= reduction;
////
////				if(supply_gap > 0.){
////					// Reduce EV charge
////					reduction = std::min(supply_gap, end_user_profiles[point_iter][sample_iter].operation.EV.BESS.capacity_scale + end_user_profiles[point_iter][sample_iter].operation.EV.BESS.normalized_scheduled_capacity_profile(0));
////					end_user_profiles[point_iter][sample_iter].operation.EV.BESS.soc_ini += reduction;
////					end_user_profiles[point_iter][sample_iter].operation.EV.BESS.normalized_scheduled_capacity_profile(0) -= reduction;
////					supply_gap -= reduction;
////				}
////			}
////		}
////	}
//
////	// Renew the expected price and demand profiles
////	int sample_num = end_user_profiles[0].size();
////	int foresight_time = agent::parameters::foresight_time();
////	double residential_ratio = agent::parameters::residential_ratio();
////
////	// Update of forecast price profile
////	auto expected_price_sorted = power_market::International_Market_Price_Sorted(tick, International_Market);
////
////	// Update of forecast demand profile and operation strategies
////	Eigen::VectorXd weight(sample_num);
////	weight = Eigen::VectorXd::Constant(sample_num, 1. / sample_num);
////	for(int point_iter = 0; point_iter < end_user_profiles.size(); ++ point_iter){
////		int bz_ID = Power_network_inform.points.bidding_zone(point_iter);
////		for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
////			end_user_profiles[point_iter][sample_iter].operation.weight = weight[sample_iter];
////
////			// Normalized default demand profile in the foresight timeframe
////			end_user_profiles[point_iter][sample_iter].operation.normalized_default_demand_profile.head(foresight_time - 1) = end_user_profiles[point_iter][sample_iter].operation.normalized_default_demand_profile.tail(foresight_time - 1);
////			end_user_profiles[point_iter][sample_iter].operation.normalized_default_demand_profile(foresight_time - 1) = residential_ratio * Power_network_inform.points.nominal_mean_demand_field(point_iter, tick + foresight_time - 1);
////
////			// PV
////			end_user_profiles[point_iter][sample_iter].operation.normalized_default_PV_profile.head(foresight_time - 1) = end_user_profiles[point_iter][sample_iter].operation.normalized_default_PV_profile.tail(foresight_time - 1);
////			end_user_profiles[point_iter][sample_iter].operation.normalized_default_PV_profile(foresight_time - 1) = Power_network_inform.points.solar_cf(point_iter, tick + foresight_time - 1);
////			end_user_profiles[point_iter][sample_iter].operation.normalized_default_PV_profile(foresight_time - 1) *= end_user_profiles[point_iter][sample_iter].operation.PV_scale;
////
////			// Smart appliance
////			end_user_profiles[point_iter][sample_iter].operation.smart_appliance.scale = .2;
////			end_user_profiles[point_iter][sample_iter].operation.smart_appliance.flexibility_factor = .5;
////			agent::end_user::smart_appliance_schedule(expected_price_sorted[bz_ID], end_user_profiles[point_iter][sample_iter].operation.normalized_default_demand_profile, end_user_profiles[point_iter][sample_iter].operation.smart_appliance);
////
////			// EV
////			end_user_profiles[point_iter][sample_iter].operation.EV.usage_default_period = Eigen::VectorXi::Zero(foresight_time);
////			end_user_profiles[point_iter][sample_iter].operation.EV.house_default_period = Eigen::VectorXi::Ones(foresight_time);
////			agent::end_user::EV_schedule(foresight_time, expected_price_sorted[bz_ID], end_user_profiles[point_iter][sample_iter].operation.EV);
////
////			// BESS
////			agent::end_user::storage_schedule_LP_optimize(foresight_time, expected_price_sorted[bz_ID], end_user_profiles[point_iter][sample_iter].operation.BESS);
////
////			// Update schedule profile and prices
////			end_user_profiles[point_iter][sample_iter].operation.supply_inflex_price = International_Market.price_range_inflex(0);
////			end_user_profiles[point_iter][sample_iter].operation.supply_flex_price = International_Market.price_range_inflex(0);
////			end_user_profiles[point_iter][sample_iter].operation.demand_inflex_price = International_Market.price_range_inflex(1);
////			end_user_profiles[point_iter][sample_iter].operation.demand_flex_price = International_Market.price_range_inflex(1);
////			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile = end_user_profiles[point_iter][sample_iter].operation.normalized_default_demand_profile;
////			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile *= (1. - end_user_profiles[point_iter][sample_iter].operation.smart_appliance.scale);
////			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile -= end_user_profiles[point_iter][sample_iter].operation.normalized_default_PV_profile;
////			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile = end_user_profiles[point_iter][sample_iter].operation.smart_appliance.normalized_scheduled_profile;
////			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile -= end_user_profiles[point_iter][sample_iter].operation.EV.BESS.normalized_scheduled_capacity_profile;
////			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile -= end_user_profiles[point_iter][sample_iter].operation.BESS.normalized_scheduled_capacity_profile;
////		}
////	}
//}

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
			Power_market_inform.DSO_Markets[DSO_ID].submitted_demand.col(point_ID) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand;
		}
	}

	// LV Power suppliers
	int pump_LV_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
	for(int agent_iter = 0; agent_iter < pump_LV_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);
		int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
		point_ID = Power_network_inform.points.in_cluster_ID(point_ID);

		Power_market_inform.DSO_Markets[DSO_ID].submitted_demand.col(point_ID) += Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_demand;
	}

	for(int DSO_iter = 0; DSO_iter < Power_market_inform.DSO_Markets.size(); ++ DSO_iter){
		Flow_Based_Market_Optimization(Power_market_inform.DSO_Markets[DSO_iter], Power_market_inform.DSO_Problems[DSO_iter]);
		DSO_Market_Results_Get(tick, Power_market_inform.DSO_Markets[DSO_iter], Power_market_inform.DSO_Problems[DSO_iter], Power_network_inform.DSO_cluster[DSO_iter], 0);
	}
}

void power_market::Filtered_bid_supply_calculation(int tick, market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
	for(int DSO_iter = 0; DSO_iter < Power_market_inform.DSO_Markets.size(); ++ DSO_iter){
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
			Power_market_inform.DSO_Markets[DSO_ID].submitted_supply.col(point_ID) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply;
		}
	}

	// LV Power suppliers
	int hydro_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);
		int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
		point_ID = Power_network_inform.points.in_cluster_ID(point_ID);

		Power_market_inform.DSO_Markets[DSO_ID].submitted_supply.col(point_ID) += Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids.redispatch_supply;
	}

	int wind_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);
		int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
		point_ID = Power_network_inform.points.in_cluster_ID(point_ID);

		Power_market_inform.DSO_Markets[DSO_ID].submitted_supply.col(point_ID) += Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.redispatch_supply;
	}

	int pump_LV_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
	for(int agent_iter = 0; agent_iter < pump_LV_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);
		int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
		point_ID = Power_network_inform.points.in_cluster_ID(point_ID);

		Power_market_inform.DSO_Markets[DSO_ID].submitted_supply.col(point_ID) += Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_supply;
	}

	for(int DSO_iter = 0; DSO_iter < Power_market_inform.DSO_Markets.size(); ++ DSO_iter){
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
			Market.confirmed_supply(tick, point_iter) = (sol_vec.segment(row_start, Market.price_intervals + 2).array().max(0)).sum();

			// Store nodal prices
			Market.confirmed_price(tick, point_iter) = Market.bidded_price_map.bidded_price(1) + rep.lagbc[row_start + 1];
			Market.confirmed_price(tick, point_iter) = int(Market.confirmed_price(tick, point_iter)) + .5;
			if(Market.confirmed_price(tick, point_iter) < Market.bidded_price_map.bidded_price(1)){
				Market.confirmed_price(tick, point_iter) = Market.bidded_price_map.bidded_price(0);
			}
			else if(Market.confirmed_price(tick, point_iter) > Market.bidded_price_map.bidded_price(Market.price_intervals)){
				Market.confirmed_price(tick, point_iter) = Market.bidded_price_map.bidded_price(Market.price_intervals + 1);
			}

			// Store ratio at nodes
			for(int price_iter = 0; price_iter < Market.price_intervals + 2; ++ price_iter){
				if(Market.bidded_price_map.bidded_price(price_iter) >= Market.confirmed_price(tick, point_iter) || price_iter == Market.price_intervals + 1){
					if(sol[row_start + price_iter] >= 0.){
						Market.confirmed_ratio_demand(point_iter) = std::min(Market.submitted_demand(price_iter, point_iter), Market.submitted_supply(price_iter, point_iter) - sol[row_start + price_iter]);
						Market.confirmed_ratio_supply(point_iter) = Market.confirmed_ratio_demand(point_iter) + sol[row_start + price_iter];
						Market.confirmed_ratio_demand(point_iter) /= Market.submitted_demand(price_iter, point_iter) + 1E-12;
					}
					else{
						Market.confirmed_ratio_supply(point_iter) = std::min(Market.submitted_supply(price_iter, point_iter), Market.submitted_demand(price_iter, point_iter) + sol[row_start + price_iter]);
						Market.confirmed_ratio_supply(point_iter) /= Market.submitted_supply(price_iter, point_iter) + 1E-12;
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
			Market.confirmed_demand(tick, point_iter) = -(sol_vec.segment(row_start, Market.price_intervals + 2).array().min(0)).sum();

			// Store nodal prices
			Market.confirmed_price(tick, point_iter) = Market.bidded_price_map.bidded_price(1) + rep.lagbc[row_start + 1];
			Market.confirmed_price(tick, point_iter) = int(Market.confirmed_price(tick, point_iter)) + .5;
			if(Market.confirmed_price(tick, point_iter) < Market.bidded_price_map.bidded_price(1)){
				Market.confirmed_price(tick, point_iter) = Market.bidded_price_map.bidded_price(0);
			}
			else if(Market.confirmed_price(tick, point_iter) > Market.bidded_price_map.bidded_price(Market.price_intervals)){
				Market.confirmed_price(tick, point_iter) = Market.bidded_price_map.bidded_price(Market.price_intervals + 1);
			}

			// Store ratio at nodes
			for(int price_iter = 0; price_iter < Market.price_intervals + 2; ++ price_iter){
				if(Market.bidded_price_map.bidded_price(price_iter) >= Market.confirmed_price(tick, point_iter) || price_iter == Market.price_intervals + 1){
					if(sol[row_start + price_iter] >= 0.){
						Market.confirmed_ratio_demand(point_iter) = std::min(Market.submitted_demand(price_iter, point_iter), Market.submitted_supply(price_iter, point_iter) - sol[row_start + price_iter]);
						Market.confirmed_ratio_demand(point_iter) /= Market.submitted_demand(price_iter, point_iter) + 1E-12;
					}
					else{
						Market.confirmed_ratio_supply(point_iter) = std::min(Market.submitted_supply(price_iter, point_iter), Market.submitted_demand(price_iter, point_iter) + sol[row_start + price_iter]);
						Market.confirmed_ratio_demand(point_iter) = Market.confirmed_ratio_supply(point_iter) - sol[row_start + price_iter];
						Market.confirmed_ratio_demand(point_iter) /= Market.submitted_demand(price_iter, point_iter) + 1E-12;
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
