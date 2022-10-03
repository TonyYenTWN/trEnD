// Source file for operational strategies of different agents
#include "agent_func.h"

namespace{
	void agent_results_set(agent::results &results){
		int Time = power_market::parameters::Time();

		results.confirmed_supply = Eigen::VectorXd(Time);
		results.confirmed_demand = Eigen::VectorXd(Time);
		results.actual_supply = Eigen::VectorXd(Time);
		results.actual_demand = Eigen::VectorXd(Time);
	}

	void agent_bids_initialization(agent::bids &bids){
		int price_interval = power_market::parameters::price_interval();

		bids.submitted_supply_inflex = Eigen::VectorXd::Zero(price_interval + 2);
		bids.submitted_demand_inflex = Eigen::VectorXd::Zero(price_interval + 2);
		bids.submitted_supply_flex = Eigen::VectorXd::Zero(price_interval + 2);
		bids.submitted_demand_flex = Eigen::VectorXd::Zero(price_interval + 2);
		bids.redispatch_supply = Eigen::VectorXd::Zero(price_interval + 2);
		bids.redispatch_demand = Eigen::VectorXd::Zero(price_interval + 2);
		bids.filtered_supply = Eigen::VectorXd::Zero(price_interval + 2);
		bids.filtered_demand = Eigen::VectorXd::Zero(price_interval + 2);
		bids.accepted_supply = Eigen::VectorXd::Zero(price_interval + 2);
		bids.accepted_demand = Eigen::VectorXd::Zero(price_interval + 2);
		bids.balancing_supply = Eigen::VectorXd::Zero(price_interval + 2);
		bids.balancing_demand = Eigen::VectorXd::Zero(price_interval + 2);
	}

	void agent_submitted_bids_scale(double scale, agent::bids &bids){
		bids.submitted_supply_inflex *= scale;
		bids.submitted_demand_inflex *= scale;
		bids.submitted_supply_flex *= scale;
		bids.submitted_demand_flex *= scale;
	}

	agent::aggregator::profiles aggregator_set(power_market::market_inform &International_Market, power_network::network_inform &Power_network_inform){
		int foresight_time = agent::aggregator::parameters::foresight_time();
		int point_num = Power_network_inform.points.bidding_zone.size();

		agent::aggregator::profiles aggregator_profiles(point_num);
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int bz_ID = Power_network_inform.points.bidding_zone(point_iter);

			aggregator_profiles[point_iter].point_ID = point_iter;
			aggregator_profiles[point_iter].price_expected_profile = International_Market.confirmed_price.col(bz_ID).head(foresight_time);
			aggregator_profiles[point_iter].price_demand_profile = International_Market.confirmed_price.col(bz_ID).head(foresight_time);
			aggregator_profiles[point_iter].price_supply_profile = International_Market.confirmed_price.col(bz_ID).head(foresight_time);
		}

		return aggregator_profiles;
	}

	agent::end_user::profiles end_user_set(power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int foresight_time = agent::end_user::parameters::foresight_time();
		int point_num = Power_network_inform.points.bidding_zone.size();
		int load_shift_time = agent::end_user::parameters::load_shift_time();
		int price_interval = power_market::parameters::price_interval();
		double residential_ratio = agent::parameters::residential_ratio();

		agent::end_user::profiles end_user_profiles(point_num);
		int sample_num = agent::end_user::parameters::sample_num();
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			end_user_profiles[point_iter] = std::vector <agent::end_user::profile> (sample_num);
		}

		// Initialization of forecast demand profile and operation strategies
		Eigen::VectorXd weight(sample_num);
		weight = Eigen::VectorXd::Constant(sample_num, 1. / sample_num);
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				// Initialization of investment parameters
				end_user_profiles[point_iter][sample_iter].investment.decision.dynamic_tariff = 1;
				end_user_profiles[point_iter][sample_iter].investment.decision.smart_appliance = 1;
				end_user_profiles[point_iter][sample_iter].investment.decision.PV = 1;
				end_user_profiles[point_iter][sample_iter].investment.decision.BESS = 1;
				end_user_profiles[point_iter][sample_iter].investment.decision.EV_self_charging = 1;
				end_user_profiles[point_iter][sample_iter].investment.decision.reverse_flow = 1;
				end_user_profiles[point_iter][sample_iter].investment.decision.redispatch = 1;
				end_user_profiles[point_iter][sample_iter].investment.decision.control_reserve = 1;

				// Initialization of operational parameters
				end_user_profiles[point_iter][sample_iter].operation.foresight_time = foresight_time;
				end_user_profiles[point_iter][sample_iter].operation.weight = weight(sample_iter);
				int load_shift_time_temp = std::min(load_shift_time, foresight_time / 2);
				end_user_profiles[point_iter][sample_iter].operation.smart_appliance.shift_time = load_shift_time_temp;
				end_user_profiles[point_iter][sample_iter].operation.BESS.soc = end_user_profiles[point_iter][sample_iter].operation.BESS.energy_scale / 2;
				end_user_profiles[point_iter][sample_iter].operation.BESS.soc *= end_user_profiles[point_iter][sample_iter].investment.decision.BESS;
				end_user_profiles[point_iter][sample_iter].operation.EV.BESS.soc = end_user_profiles[point_iter][sample_iter].operation.EV.BESS.energy_scale / 2;
				end_user_profiles[point_iter][sample_iter].operation.EV.BESS.soc *= end_user_profiles[point_iter][sample_iter].investment.decision.EV_self_charging;

				// Initialization of input profiles
				end_user_profiles[point_iter][sample_iter].operation.EV.house_default_period = Eigen::VectorXi::Ones(foresight_time);
				end_user_profiles[point_iter][sample_iter].operation.EV.default_demand_profile = Eigen::VectorXd::Zero(foresight_time);
				end_user_profiles[point_iter][sample_iter].operation.EV.default_demand_profile *= end_user_profiles[point_iter][sample_iter].investment.decision.EV_self_charging;
				end_user_profiles[point_iter][sample_iter].operation.default_demand_profile = Power_network_inform.points.nominal_mean_demand_field.row(point_iter).head(foresight_time);
				//end_user_profiles[point_iter][sample_iter].operation.default_demand_profile -= end_user_profiles[point_iter][sample_iter].operation.EV.default_demand_profile;
				end_user_profiles[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand = Eigen::VectorXd::Zero(foresight_time + load_shift_time_temp);
				for(int tick = load_shift_time_temp; tick < foresight_time + load_shift_time_temp; ++ tick){
					end_user_profiles[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand(tick) = Power_network_inform.points.nominal_mean_demand_field(point_iter, tick - load_shift_time_temp);
					end_user_profiles[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand(tick) *= end_user_profiles[point_iter][sample_iter].investment.decision.smart_appliance * end_user_profiles[point_iter][sample_iter].operation.smart_appliance.scale;
				}
				end_user_profiles[point_iter][sample_iter].operation.default_demand_profile *= 1. - end_user_profiles[point_iter][sample_iter].investment.decision.smart_appliance * end_user_profiles[point_iter][sample_iter].operation.smart_appliance.scale;
				end_user_profiles[point_iter][sample_iter].operation.default_PV_profile = Eigen::VectorXd::Zero(foresight_time);
				end_user_profiles[point_iter][sample_iter].operation.price_demand_profile = Power_market_inform.agent_profiles.aggregators[point_iter].price_demand_profile;
				end_user_profiles[point_iter][sample_iter].operation.price_supply_profile = Power_market_inform.agent_profiles.aggregators[point_iter].price_supply_profile;

				// Set the LP problem
				agent::end_user::end_user_LP_set(end_user_profiles[point_iter][sample_iter]);

				// Set bids and results information
				agent_bids_initialization(end_user_profiles[point_iter][sample_iter].operation.bids);
				//agent_results_set(end_user_profiles[point_iter][sample_iter].operation.results);

				// Optimization and update process variables
				agent::end_user::end_user_LP_optimize(0, end_user_profiles[point_iter][sample_iter]);

				// Scale the bids correctly
				double scale = end_user_profiles[point_iter][sample_iter].operation.weight;
				scale *= agent::parameters::residential_ratio();
				scale *= Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
				agent_submitted_bids_scale(scale, end_user_profiles[point_iter][sample_iter].operation.bids);
			}
		}

		return end_user_profiles;
	}

	agent::industrial::profiles industrial_set(power_network::network_inform &Power_network_inform){
		int point_num = Power_network_inform.points.bidding_zone.size();
		int price_interval = power_market::parameters::price_interval();

		agent::industrial::profiles industrial_profiles;
		industrial_profiles.HV.reserve(point_num);
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int point_ID = point_iter;
			agent::industrial::profile profile_temp;
			profile_temp.point_ID = point_ID;

			double bid_inflex_industrial = Power_network_inform.points.nominal_mean_demand_field(point_iter, 0);
			bid_inflex_industrial *= Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
			bid_inflex_industrial *= 1. - agent::parameters::residential_ratio();
			double bid_flex_industrial = bid_inflex_industrial;
			bid_inflex_industrial *= 1. - agent::industrial::flexible_ratio();
			bid_flex_industrial *= agent::industrial::flexible_ratio();
			bid_flex_industrial /= price_interval;

			// Set bids information
			agent_bids_initialization(profile_temp.bids);
			//agent_results_set(profile_temp.results);
			profile_temp.bids.submitted_demand_flex(price_interval + 1) = bid_inflex_industrial;
			profile_temp.bids.submitted_demand_flex.segment(1, price_interval) = Eigen::VectorXd::Constant(price_interval, bid_flex_industrial);
			profile_temp.bids.redispatch_demand = profile_temp.bids.submitted_demand_flex;
			industrial_profiles.HV.push_back(profile_temp);
		}

		return industrial_profiles;
	}

	agent::power_supplier::profiles power_supplier_set(power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int point_num = Power_network_inform.points.bidding_zone.size();
		int price_interval = power_market::parameters::price_interval();
		int hydro_num = Power_network_inform.plants.hydro.cap.size();
		int wind_num = Power_network_inform.plants.wind.cap.size();
		double cutoff_power = agent::power_supplier::parameters::cutoff_power();

		agent::power_supplier::profiles power_supplier_profiles;
		power_supplier_profiles.hydro.HV_hybrid.reserve(hydro_num);
		power_supplier_profiles.hydro.HV_plant.reserve(hydro_num);
		power_supplier_profiles.hydro.LV_hybrid.reserve(hydro_num);
		power_supplier_profiles.hydro.LV_plant.reserve(hydro_num);
		power_supplier_profiles.pump_storage.HV.reserve(hydro_num);
		power_supplier_profiles.pump_storage.LV.reserve(hydro_num);
		for(int agent_iter = 0; agent_iter < hydro_num; ++ agent_iter){
			int x_ID = int((Power_network_inform.plants.hydro.x(agent_iter) - Power_network_inform.points.x.minCoeff()) / Power_network_inform.points.grid_length + .5);
			int y_ID = int((Power_network_inform.plants.hydro.y(agent_iter) - Power_network_inform.points.y.minCoeff()) / Power_network_inform.points.grid_length + .5);
			int point_ID = Power_network_inform.points.coordinate_grid(x_ID, y_ID);
			if(point_ID == -1){
				continue;
			}
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);
			Eigen::VectorXd bid_vec = Power_market_inform.International_Market.merit_order_curve.col(bz_ID);
			bid_vec *= Power_network_inform.plants.hydro.cap(agent_iter);
			bid_vec /= (Power_market_inform.International_Market.merit_order_curve.col(bz_ID).sum());

			if(Power_network_inform.plants.hydro.type(agent_iter) < 4){
				agent::power_supplier::plant_profile profile_temp;
				profile_temp.point_ID = point_ID;
				profile_temp.cap = Power_network_inform.plants.hydro.cap(agent_iter);

				// Set bids information
				agent_bids_initialization(profile_temp.bids);
				//agent_results_set(profile_temp.results);
				profile_temp.bids.submitted_supply_flex = bid_vec;
				profile_temp.bids.redispatch_supply = bid_vec;

				// High voltage power plants connect directly to transmission network
				if(profile_temp.cap >= cutoff_power){
					power_supplier_profiles.hydro.HV_plant.push_back(profile_temp);
				}
				// Low voltage power plants feed into distribution network
				else{
					power_supplier_profiles.hydro.LV_plant.push_back(profile_temp);
				}
			}
			else if(Power_network_inform.plants.hydro.type(agent_iter) == 5){
				agent::power_supplier::storage_profile profile_temp;
				profile_temp.point_ID = point_ID;
				profile_temp.cap = Power_network_inform.plants.hydro.cap(agent_iter);

				// Set bids information
				agent_bids_initialization(profile_temp.bids);
				//agent_results_set(profile_temp.results);
				profile_temp.bids.submitted_supply_flex = bid_vec;
				profile_temp.bids.redispatch_supply = bid_vec;

				// High voltage power plants connect directly to transmission network
				if(profile_temp.cap >= cutoff_power){
					power_supplier_profiles.pump_storage.HV.push_back(profile_temp);
				}
				// Low voltage power plants feed into distribution network
				else{
					power_supplier_profiles.pump_storage.LV.push_back(profile_temp);
				}
			}
		}

		power_supplier_profiles.wind.HV_hybrid.reserve(wind_num);
		power_supplier_profiles.wind.HV_plant.reserve(wind_num);
		power_supplier_profiles.wind.LV_hybrid.reserve(wind_num);
		power_supplier_profiles.wind.LV_plant.reserve(wind_num);
		for(int agent_iter = 0; agent_iter < wind_num; ++ agent_iter){
			int x_ID = int((Power_network_inform.plants.wind.x(agent_iter) - Power_network_inform.points.x.minCoeff()) / Power_network_inform.points.grid_length + .5);
			int y_ID = int((Power_network_inform.plants.wind.y(agent_iter) - Power_network_inform.points.y.minCoeff()) / Power_network_inform.points.grid_length + .5);
			int point_ID = Power_network_inform.points.coordinate_grid(x_ID, y_ID);
			if(point_ID == -1){
				continue;
			}
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);
			int price_supply_flex_ID = 0;
			double bid_quan = Power_network_inform.points.wind_on_cf(point_ID, 0) * Power_network_inform.plants.wind.cap(agent_iter);

			agent::power_supplier::plant_profile profile_temp;
			profile_temp.point_ID = point_ID;
			profile_temp.cap = Power_network_inform.plants.wind.cap(agent_iter);

			// Set bids information
			agent_bids_initialization(profile_temp.bids);
			//agent_results_set(profile_temp.results);
			profile_temp.bids.submitted_supply_flex(price_supply_flex_ID) = bid_quan;
			profile_temp.bids.redispatch_supply = profile_temp.bids.submitted_supply_flex;

			// High voltage power plants connect directly to transmission network
			if(profile_temp.cap >= cutoff_power){
				power_supplier_profiles.hydro.HV_plant.push_back(profile_temp);
			}
			// Low voltage power plants feed into distribution network
			else{
				power_supplier_profiles.hydro.LV_plant.push_back(profile_temp);
			}
		}

		return power_supplier_profiles;
	}

	void end_user_redispatch_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int point_num = Power_network_inform.points.bidding_zone.size();
		int sample_num = agent::end_user::parameters::sample_num();
		int price_interval = power_market::parameters::price_interval();

		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int bz_ID = Power_network_inform.points.bidding_zone(point_iter);
			double marginal_price = Power_market_inform.International_Market.confirmed_price(tick, bz_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				// Inflexible bids have first / last priority for dispatch
				if(marginal_price_ID > 0){
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(0) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_inflex.head(marginal_price_ID).sum();
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(0) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_inflex.head(marginal_price_ID).sum();
				}
				if(marginal_price_ID < price_interval + 1){
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(price_interval + 1) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_inflex.tail(price_interval + 1 - marginal_price_ID).sum();
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(price_interval + 1) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_inflex.tail(price_interval + 1 - marginal_price_ID).sum();
				}
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(0) += (1. - Power_market_inform.International_Market.confirmed_ratio_demand(bz_ID)) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_inflex(marginal_price_ID);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(0) += Power_market_inform.International_Market.confirmed_ratio_supply(bz_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_inflex(marginal_price_ID);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(price_interval + 1) += Power_market_inform.International_Market.confirmed_ratio_demand(bz_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_inflex(marginal_price_ID);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(price_interval + 1) += (1. - Power_market_inform.International_Market.confirmed_ratio_supply(bz_ID)) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_inflex(marginal_price_ID);

				// Flexible bids have redispatcch priority in between
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_flex;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_flex;
			}
		}
	}

	void end_user_balancing_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, bool DSO_filter_flag){
		int point_num = Power_network_inform.points.bidding_zone.size();
		int sample_num = agent::end_user::parameters::sample_num();
		int price_interval = power_market::parameters::price_interval();

		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int node_ID = Power_network_inform.points.node(point_iter);
			double marginal_price = Power_market_inform.TSO_Market.confirmed_price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				if(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].investment.decision.control_reserve){
					if(DSO_filter_flag){

					}
					else{
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand;
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply;
					}
				}
				else{
					if(DSO_filter_flag){

					}
					else{
						if(marginal_price_ID > 0){
							Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand(0) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand.head(marginal_price_ID).sum();
							Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply(0) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply.head(marginal_price_ID).sum();
						}
						if(marginal_price_ID < price_interval + 1){
							Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand(price_interval + 1) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand.tail(price_interval + 1 - marginal_price_ID).sum();
							Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply(price_interval + 1) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply.tail(price_interval + 1 - marginal_price_ID).sum();
						}
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand(0) += (1. - Power_market_inform.TSO_Market.confirmed_ratio_demand(node_ID)) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(marginal_price_ID);
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply(0) += Power_market_inform.TSO_Market.confirmed_ratio_supply(node_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(marginal_price_ID);
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply(price_interval + 1) += Power_market_inform.TSO_Market.confirmed_ratio_demand(node_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(marginal_price_ID);
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply(price_interval + 1) += (1. - Power_market_inform.TSO_Market.confirmed_ratio_supply(node_ID)) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(marginal_price_ID);
					}
				}
			}
		}
	}
}

void agent::agents_set(power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
	Power_market_inform.agent_profiles.aggregators = aggregator_set(Power_market_inform.International_Market, Power_network_inform);
	Power_market_inform.agent_profiles.end_users = end_user_set(Power_market_inform, Power_network_inform);
	Power_market_inform.agent_profiles.industrial = industrial_set(Power_network_inform);
	Power_market_inform.agent_profiles.power_supplier = power_supplier_set(Power_market_inform, Power_network_inform);
}

void agent::agents_redispatch_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
	end_user_redispatch_update(tick, Power_market_inform, Power_network_inform);
}

void agent::agents_balancing_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, bool DSO_filter_flag){
	end_user_balancing_update(tick, Power_market_inform, Power_network_inform, DSO_filter_flag);
}
