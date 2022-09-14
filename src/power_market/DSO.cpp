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
					DSO_Markets[DSO_iter].network.admittance.push_back(admittance(row_iter, col_iter));
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
		DSO_Markets[DSO_iter].filtered_price_supply = Eigen::MatrixXd::Zero(Time, DSO_Markets[DSO_iter].num_zone);
		DSO_Markets[DSO_iter].filtered_price_demand = Eigen::MatrixXd::Zero(Time, DSO_Markets[DSO_iter].num_zone);
		DSO_Markets[DSO_iter].filtered_ratio_supply = Eigen::MatrixXd::Zero(Time, DSO_Markets[DSO_iter].num_zone);
		DSO_Markets[DSO_iter].filtered_ratio_demand = Eigen::MatrixXd::Zero(Time, DSO_Markets[DSO_iter].num_zone);
//		DSO_Markets[DSO_iter].confirmed_supply = Eigen::MatrixXd::Zero(Time, DSO_Markets[DSO_iter].num_zone);
//		DSO_Markets[DSO_iter].confirmed_demand = Eigen::MatrixXd::Zero(Time, DSO_Markets[DSO_iter].num_zone);
//		DSO_Markets[DSO_iter].confirmed_price = Eigen::MatrixXd(Time, DSO_Markets[DSO_iter].num_zone);
//		DSO_Markets[DSO_iter].network.confirmed_power = Eigen::MatrixXd(Time, DSO_Markets[DSO_iter].network.num_edges);
//		DSO_Markets[DSO_iter].network.confirmed_voltage = Eigen::MatrixXd(Time, DSO_Markets[DSO_iter].network.num_vertice);
	}
}

agent::end_user::profiles power_market::DSO_agents_set(market_inform &International_Market, power_network::network_inform &Power_network_inform){
	int foresight_time = agent::parameters::foresight_time();
	double residential_ratio = agent::parameters::residential_ratio();

	agent::end_user::profiles end_user_profiles(Power_network_inform.points.bidding_zone.rows());
	int sample_num = agent::parameters::sample_num();
	for(int point_iter = 0; point_iter < end_user_profiles.size(); ++ point_iter){
		end_user_profiles[point_iter] = std::vector <agent::end_user::profile> (sample_num);
	}

	// Initialization of forecast price profile
	auto expected_price_sorted = power_market::International_Market_Price_Sorted(0, International_Market);

	// Initialization of forecast demand profile and operation strategies
	Eigen::VectorXd weight(sample_num);
	weight = Eigen::VectorXd::Constant(sample_num, 1. / sample_num);
	auto Problem = agent::end_user::storage_schedule_LP_mold(foresight_time);
	for(int point_iter = 0; point_iter < end_user_profiles.size(); ++ point_iter){
		int bz_ID = Power_network_inform.points.bidding_zone(point_iter);

		for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
			end_user_profiles[point_iter][sample_iter].operation.weight = weight[sample_iter];

			// Normalized default demand profile in the foresight timeframe
			end_user_profiles[point_iter][sample_iter].operation.normalized_default_demand_profile = residential_ratio * (Power_network_inform.points.nominal_mean_demand_field.row(point_iter)).head(foresight_time);

			// PV
			end_user_profiles[point_iter][sample_iter].operation.PV_scale = .01;
			end_user_profiles[point_iter][sample_iter].operation.normalized_default_PV_profile = Eigen::VectorXd(foresight_time);
			for(int tick = 0; tick < foresight_time; ++ tick){
				end_user_profiles[point_iter][sample_iter].operation.normalized_default_PV_profile(tick) = Power_network_inform.points.solar_cf(point_iter, tick);
				end_user_profiles[point_iter][sample_iter].operation.normalized_default_PV_profile(tick) *= end_user_profiles[point_iter][sample_iter].operation.PV_scale;
			}

			// Smart appliance
			end_user_profiles[point_iter][sample_iter].operation.smart_appliance.scale = 0.;
			end_user_profiles[point_iter][sample_iter].operation.smart_appliance.flexibility_factor = .5;
			agent::end_user::smart_appliance_schedule(expected_price_sorted[bz_ID], end_user_profiles[point_iter][sample_iter].operation.normalized_default_demand_profile, end_user_profiles[point_iter][sample_iter].operation.smart_appliance);

			// EV
			end_user_profiles[point_iter][sample_iter].operation.EV.BESS.soc_ini = .5 * end_user_profiles[point_iter][sample_iter].operation.EV.BESS.energy_scale;
			end_user_profiles[point_iter][sample_iter].operation.EV.usage_default_period = Eigen::VectorXi::Zero(foresight_time);
			end_user_profiles[point_iter][sample_iter].operation.EV.house_default_period = Eigen::VectorXi::Ones(foresight_time);
			agent::end_user::EV_schedule(foresight_time, expected_price_sorted[bz_ID], end_user_profiles[point_iter][sample_iter].operation.EV);

			// BESS
			end_user_profiles[point_iter][sample_iter].operation.BESS.soc_ini = .5 * end_user_profiles[point_iter][sample_iter].operation.BESS.energy_scale;
			end_user_profiles[point_iter][sample_iter].operation.BESS.Problem = Problem;
			agent::end_user::storage_schedule_LP_optimize(foresight_time, expected_price_sorted[bz_ID], end_user_profiles[point_iter][sample_iter].operation.BESS);

			// Update schedule profile and prices
			end_user_profiles[point_iter][sample_iter].operation.supply_inflex_price_ID = 0;
			end_user_profiles[point_iter][sample_iter].operation.supply_flex_price_ID = 0;
			end_user_profiles[point_iter][sample_iter].operation.demand_inflex_price_ID = International_Market.price_intervals + 1;
			end_user_profiles[point_iter][sample_iter].operation.demand_flex_price_ID = International_Market.price_intervals + 1;
			end_user_profiles[point_iter][sample_iter].operation.supply_inflex_price = International_Market.bidded_price(end_user_profiles[point_iter][sample_iter].operation.supply_inflex_price_ID);
			end_user_profiles[point_iter][sample_iter].operation.supply_flex_price = International_Market.bidded_price(end_user_profiles[point_iter][sample_iter].operation.supply_flex_price_ID);
			end_user_profiles[point_iter][sample_iter].operation.demand_inflex_price = International_Market.bidded_price(end_user_profiles[point_iter][sample_iter].operation.demand_inflex_price_ID);
			end_user_profiles[point_iter][sample_iter].operation.demand_flex_price = International_Market.bidded_price(end_user_profiles[point_iter][sample_iter].operation.demand_flex_price_ID);
			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile = end_user_profiles[point_iter][sample_iter].operation.normalized_default_demand_profile;
			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile *= (1. - end_user_profiles[point_iter][sample_iter].operation.smart_appliance.scale);
			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile -= end_user_profiles[point_iter][sample_iter].operation.normalized_default_PV_profile;
			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile = end_user_profiles[point_iter][sample_iter].operation.smart_appliance.normalized_scheduled_profile;
			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile -= end_user_profiles[point_iter][sample_iter].operation.EV.BESS.normalized_scheduled_capacity_profile;
			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile -= end_user_profiles[point_iter][sample_iter].operation.BESS.normalized_scheduled_capacity_profile;
			//std::cout << end_user_profiles[point_iter][sample_iter].operation.normalized_default_demand_profile.transpose() << "\n";
			//std::cout << end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile.transpose() << "\n\n";
		}
	}

	return end_user_profiles;
}

void power_market::DSO_agents_update(int tick, agent::end_user::profiles &end_user_profiles, market_inform &TSO_Market, market_inform &International_Market, power_network::network_inform &Power_network_inform){
	// Update state of agents from the previous tick
	for(int point_iter = 0; point_iter < end_user_profiles.size(); ++ point_iter){
		int node_ID = Power_network_inform.points.node(point_iter);
		int sample_num = agent::parameters::sample_num();
		double node_price = TSO_Market.actual_price(tick - 1, node_ID);

		for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
			end_user_profiles[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand = 0.;
			end_user_profiles[point_iter][sample_iter].operation.EV.BESS.soc_ini -= end_user_profiles[point_iter][sample_iter].operation.EV.BESS.normalized_scheduled_capacity_profile(0);
			end_user_profiles[point_iter][sample_iter].operation.BESS.soc_ini -= end_user_profiles[point_iter][sample_iter].operation.BESS.normalized_scheduled_capacity_profile(0);

			// Demand settlement
			// Order of reduction: BESS charge, EV charge, smart appliance, inflexible demand
			double marginal_demand = 0.;
			marginal_demand += (node_price == end_user_profiles[point_iter][sample_iter].operation.demand_flex_price) * (end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile(0) >= 0.) * end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile(0);
			marginal_demand += (node_price == end_user_profiles[point_iter][sample_iter].operation.demand_inflex_price) * (end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile(0) >= 0.) * end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile(0);
			double demand_gap = (1. - TSO_Market.actual_price_ratio(tick - 1, node_ID)) * marginal_demand;
			demand_gap += (node_price > end_user_profiles[point_iter][sample_iter].operation.demand_flex_price) * (end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile(0) >= 0.) * end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile(0);
			demand_gap += (node_price > end_user_profiles[point_iter][sample_iter].operation.demand_inflex_price) * (end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile(0) >= 0.) * end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile(0);

			if(demand_gap > 0.){
				// Reduce BESS charge
				double reduction = std::min(demand_gap, end_user_profiles[point_iter][sample_iter].operation.BESS.capacity_scale - end_user_profiles[point_iter][sample_iter].operation.BESS.normalized_scheduled_capacity_profile(0));
				end_user_profiles[point_iter][sample_iter].operation.BESS.soc_ini -= reduction;
				end_user_profiles[point_iter][sample_iter].operation.BESS.normalized_scheduled_capacity_profile(0) += reduction;
				demand_gap -= reduction;

				if(demand_gap > 0.){
					// Reduce EV charge
					reduction = std::min(demand_gap, end_user_profiles[point_iter][sample_iter].operation.EV.BESS.capacity_scale - end_user_profiles[point_iter][sample_iter].operation.EV.BESS.normalized_scheduled_capacity_profile(0));
					end_user_profiles[point_iter][sample_iter].operation.EV.BESS.soc_ini -= reduction;
					end_user_profiles[point_iter][sample_iter].operation.EV.BESS.normalized_scheduled_capacity_profile(0) += reduction;
					demand_gap -= reduction;

					if(demand_gap > 0.){
						// Reduce smart appliance
						reduction = std::min(demand_gap, end_user_profiles[point_iter][sample_iter].operation.smart_appliance.normalized_scheduled_profile(0));
						end_user_profiles[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand += reduction;
						demand_gap -= reduction;
					}
				}
			}

			// Supply settlement
			// Order of reduction: BESS discharge, EV discharge, inflexible supply
			double marginal_supply = 0.;
			marginal_supply += (node_price == end_user_profiles[point_iter][sample_iter].operation.supply_flex_price) * (end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile(0) < 0.) * end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile(0);
			marginal_supply += (node_price == end_user_profiles[point_iter][sample_iter].operation.supply_inflex_price) * (end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile(0) < 0.) * end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile(0);
			double supply_gap = (1. - TSO_Market.actual_price_ratio(tick - 1, node_ID)) * marginal_supply;
			supply_gap += (node_price > end_user_profiles[point_iter][sample_iter].operation.supply_flex_price) * (end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile(0) < 0.) * end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile(0);
			supply_gap += (node_price > end_user_profiles[point_iter][sample_iter].operation.supply_inflex_price) * (end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile(0) < 0.) * end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile(0);

			if(supply_gap > 0.){
				// Reduce BESS charge
				double reduction = std::min(supply_gap, end_user_profiles[point_iter][sample_iter].operation.BESS.capacity_scale + end_user_profiles[point_iter][sample_iter].operation.BESS.normalized_scheduled_capacity_profile(0));
				end_user_profiles[point_iter][sample_iter].operation.BESS.soc_ini += reduction;
				end_user_profiles[point_iter][sample_iter].operation.BESS.normalized_scheduled_capacity_profile(0) -= reduction;
				supply_gap -= reduction;

				if(supply_gap > 0.){
					// Reduce EV charge
					reduction = std::min(supply_gap, end_user_profiles[point_iter][sample_iter].operation.EV.BESS.capacity_scale + end_user_profiles[point_iter][sample_iter].operation.EV.BESS.normalized_scheduled_capacity_profile(0));
					end_user_profiles[point_iter][sample_iter].operation.EV.BESS.soc_ini += reduction;
					end_user_profiles[point_iter][sample_iter].operation.EV.BESS.normalized_scheduled_capacity_profile(0) -= reduction;
					supply_gap -= reduction;
				}
			}
		}
	}

	// Renew the expected price and demand profiles
	int sample_num = end_user_profiles[0].size();
	int foresight_time = agent::parameters::foresight_time();
	double residential_ratio = agent::parameters::residential_ratio();

	// Update of forecast price profile
	auto expected_price_sorted = power_market::International_Market_Price_Sorted(tick, International_Market);

	// Update of forecast demand profile and operation strategies
	Eigen::VectorXd weight(sample_num);
	weight = Eigen::VectorXd::Constant(sample_num, 1. / sample_num);
	for(int point_iter = 0; point_iter < end_user_profiles.size(); ++ point_iter){
		int bz_ID = Power_network_inform.points.bidding_zone(point_iter);
		for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
			end_user_profiles[point_iter][sample_iter].operation.weight = weight[sample_iter];

			// Normalized default demand profile in the foresight timeframe
			end_user_profiles[point_iter][sample_iter].operation.normalized_default_demand_profile.head(foresight_time - 1) = end_user_profiles[point_iter][sample_iter].operation.normalized_default_demand_profile.tail(foresight_time - 1);
			end_user_profiles[point_iter][sample_iter].operation.normalized_default_demand_profile(foresight_time - 1) = residential_ratio * Power_network_inform.points.nominal_mean_demand_field(point_iter, tick + foresight_time - 1);

			// PV
			end_user_profiles[point_iter][sample_iter].operation.normalized_default_PV_profile.head(foresight_time - 1) = end_user_profiles[point_iter][sample_iter].operation.normalized_default_PV_profile.tail(foresight_time - 1);
			end_user_profiles[point_iter][sample_iter].operation.normalized_default_PV_profile(foresight_time - 1) = Power_network_inform.points.solar_cf(point_iter, tick + foresight_time - 1);
			end_user_profiles[point_iter][sample_iter].operation.normalized_default_PV_profile(foresight_time - 1) *= end_user_profiles[point_iter][sample_iter].operation.PV_scale;

			// Smart appliance
			end_user_profiles[point_iter][sample_iter].operation.smart_appliance.scale = .2;
			end_user_profiles[point_iter][sample_iter].operation.smart_appliance.flexibility_factor = .5;
			agent::end_user::smart_appliance_schedule(expected_price_sorted[bz_ID], end_user_profiles[point_iter][sample_iter].operation.normalized_default_demand_profile, end_user_profiles[point_iter][sample_iter].operation.smart_appliance);

			// EV
			end_user_profiles[point_iter][sample_iter].operation.EV.usage_default_period = Eigen::VectorXi::Zero(foresight_time);
			end_user_profiles[point_iter][sample_iter].operation.EV.house_default_period = Eigen::VectorXi::Ones(foresight_time);
			agent::end_user::EV_schedule(foresight_time, expected_price_sorted[bz_ID], end_user_profiles[point_iter][sample_iter].operation.EV);

			// BESS
			agent::end_user::storage_schedule_LP_optimize(foresight_time, expected_price_sorted[bz_ID], end_user_profiles[point_iter][sample_iter].operation.BESS);

			// Update schedule profile and prices
			end_user_profiles[point_iter][sample_iter].operation.supply_inflex_price = International_Market.price_range_inflex(0);
			end_user_profiles[point_iter][sample_iter].operation.supply_flex_price = International_Market.price_range_inflex(0);
			end_user_profiles[point_iter][sample_iter].operation.demand_inflex_price = International_Market.price_range_inflex(1);
			end_user_profiles[point_iter][sample_iter].operation.demand_flex_price = International_Market.price_range_inflex(1);
			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile = end_user_profiles[point_iter][sample_iter].operation.normalized_default_demand_profile;
			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile *= (1. - end_user_profiles[point_iter][sample_iter].operation.smart_appliance.scale);
			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_inflex_profile -= end_user_profiles[point_iter][sample_iter].operation.normalized_default_PV_profile;
			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile = end_user_profiles[point_iter][sample_iter].operation.smart_appliance.normalized_scheduled_profile;
			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile -= end_user_profiles[point_iter][sample_iter].operation.EV.BESS.normalized_scheduled_capacity_profile;
			end_user_profiles[point_iter][sample_iter].operation.normalized_scheduled_residual_demand_flex_profile -= end_user_profiles[point_iter][sample_iter].operation.BESS.normalized_scheduled_capacity_profile;
		}
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

void power_market::DSO_Market_Results_Get(int tick, market_inform &Market, alglib::minlpstate &Problem, power_network::DSO_cluster &DSO_cluster, bool supply){
	alglib::real_1d_array sol;
	alglib::minlpreport rep;
	alglib::minlpresults(Problem, sol, rep);
	Eigen::VectorXd sol_vec = Eigen::Map <Eigen::VectorXd> (sol.getcontent(), sol.length());

	if(supply){
		// Store filtered supply bids at spatial points
		for(int point_iter = 0; point_iter < DSO_cluster.points_ID.size(); ++ point_iter){
			int row_start = 2 * Market.network.num_vertice + point_iter * (Market.price_intervals + 2);
			Market.filtered_supply.col(point_iter) = sol_vec.segment(row_start, Market.price_intervals + 2).array().max(0);
			Market.filtered_price_supply(tick, point_iter) = Market.bidded_price(0) + rep.lagbc[row_start];
			Market.filtered_price_supply(tick, point_iter) = std::min(Market.filtered_price_supply(tick, point_iter), Market.price_range_inflex(1));
			Market.filtered_price_supply(tick, point_iter) = std::max(Market.filtered_price_supply(tick, point_iter), Market.price_range_inflex(0));

//			for(int price_iter = 0; price_iter < Market.price_intervals + 2; ++ price_iter ){
//				std::cout << rep.lagbc[row_start + price_iter] << "\t";
//			}
//			std::cout << "\n\n";
			//std::cout << Market.filtered_supply.col(point_iter).transpose() << "\n\n";
			//std::cout << (Market.submitted_supply.col(point_iter) - Market.filtered_supply.col(point_iter)).transpose() << "\n\n";
		}
	}
	else{
		// Store filtered demand bids at spatial points
		for(int point_iter = 0; point_iter < DSO_cluster.points_ID.size(); ++ point_iter){
			int row_start = 2 * Market.network.num_vertice + point_iter * (Market.price_intervals + 2);
			Market.filtered_demand.col(point_iter) = -(sol_vec.segment(row_start, Market.price_intervals + 2).array().min(0));
			Market.filtered_price_demand(tick, point_iter) = Market.bidded_price(0) + rep.lagbc[row_start];
			Market.filtered_price_demand(tick, point_iter) = std::min(Market.filtered_price_demand(tick, point_iter), Market.price_range_inflex(1));
			Market.filtered_price_demand(tick, point_iter) = std::max(Market.filtered_price_demand(tick, point_iter), Market.price_range_inflex(0));

			//std::cout << Market.submitted_demand.col(point_iter).transpose() << "\n\n";
			//std::cout << Market.filtered_demand.col(point_iter).transpose() << "\n\n";
			//std::cout << (Market.submitted_demand.col(point_iter) - Market.filtered_demand.col(point_iter)).transpose() << "\n\n";
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
