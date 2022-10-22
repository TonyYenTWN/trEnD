// Source file for operational strategies of different agents
#include "agent_func.h"

namespace{
	void agent_settlement_set(agent::settlement &settlement){
		settlement.cost.balancing = 0.;
		settlement.cost.EOM = 0.;
		settlement.cost.redispatch = 0.;
		settlement.price.balancing = 0.;
		settlement.price.EOM = 0.;
		settlement.price.redispatch = 0.;
		settlement.utility.balancing = 0.;
		settlement.utility.EOM = 0.;
		settlement.utility.redispatch = 0.;
		settlement.volume_demand.balancing = 0.;
		settlement.volume_demand.EOM = 0.;
		settlement.volume_demand.redispatch = 0.;
		settlement.volume_supply.balancing = 0.;
		settlement.volume_supply.EOM = 0.;
		settlement.volume_supply.redispatch = 0.;
	}

	void agent_results_set(agent::results &results){
		results.cleared_supply = 0.;
		results.cleared_demand = 0.;
		results.confirmed_supply = 0.;
		results.confirmed_demand = 0.;
		results.actual_supply = 0.;
		results.actual_demand = 0.;
	}

	void agent_bids_initialization(agent::bids &bids){
		int price_interval = power_market::parameters::price_interval();

		bids.submitted_supply_inflex = Eigen::VectorXd::Zero(price_interval + 2);
		bids.submitted_demand_inflex = Eigen::VectorXd::Zero(price_interval + 2);
		bids.submitted_supply_flex = Eigen::VectorXd::Zero(price_interval + 2);
		bids.submitted_demand_flex = Eigen::VectorXd::Zero(price_interval + 2);
		bids.redispatch_supply = Eigen::VectorXd::Zero(price_interval + 2);
		bids.redispatch_demand = Eigen::VectorXd::Zero(price_interval + 2);
		bids.filter_supply = Eigen::VectorXd::Zero(price_interval + 2);
		bids.filter_demand = Eigen::VectorXd::Zero(price_interval + 2);
		bids.balancing_supply = Eigen::VectorXd::Zero(price_interval + 2);
		bids.balancing_demand = Eigen::VectorXd::Zero(price_interval + 2);
	}

	void agent_submitted_bids_scale(double scale, agent::bids &bids){
		bids.submitted_supply_inflex *= scale;
		bids.submitted_demand_inflex *= scale;
		bids.submitted_supply_flex *= scale;
		bids.submitted_demand_flex *= scale;
	}

	void agent_scheduled_results_calculation(int bz_ID, int node_ID, int marginal_price_ID, int original_price_ID, power_market::market_whole_inform &Power_market_inform, agent::bids &bids, agent::results &results){
		int price_interval = power_market::parameters::price_interval();

		if(marginal_price_ID > 0){
			results.confirmed_supply += bids.redispatch_supply.head(marginal_price_ID).sum();
		}
		if(marginal_price_ID < price_interval + 1){
			results.confirmed_demand += bids.redispatch_demand.tail(price_interval + 1 - marginal_price_ID).sum();
		}
		results.confirmed_supply += Power_market_inform.TSO_Market.confirmed.ratio_supply(node_ID) * bids.redispatch_supply(marginal_price_ID);
		results.confirmed_demand += Power_market_inform.TSO_Market.confirmed.ratio_demand(node_ID) * bids.redispatch_demand(marginal_price_ID);

		if(original_price_ID > 0){
			results.cleared_supply += bids.filter_supply.head(original_price_ID).sum();
		}
		if(original_price_ID < price_interval + 1){
			results.cleared_demand += bids.filter_demand.tail(price_interval + 1 - original_price_ID).sum();
		}
		results.cleared_supply += Power_market_inform.International_Market.confirmed.ratio_supply(bz_ID)  * bids.filter_supply(original_price_ID);
		results.cleared_demand += Power_market_inform.International_Market.confirmed.ratio_demand(bz_ID) * bids.filter_demand(original_price_ID);
	}

	void agent_redispatch_settlement_calculation(int tick, int node_ID, double original_price, power_market::market_whole_inform &Power_market_inform, agent::bids &bids, agent::results &results, agent::settlement &settlement){
		int price_interval = power_market::parameters::price_interval();
		double redispatch_price_max = power_market::parameters::redispatch_price_max();

		// Settlement of redispatch
		// Supply side
		double cleared_supply_gap = results.cleared_supply;
		double confirmed_supply_gap = results.confirmed_supply;
		double min_supply_gap = std::min(cleared_supply_gap, confirmed_supply_gap);
		double max_supply_gap = std::max(cleared_supply_gap, confirmed_supply_gap);
		bool reduced_flag_supply = (confirmed_supply_gap == min_supply_gap);
		double margin_quan_supply;
		int margin_ID_supply;
		for(int price_iter = 0; price_iter < price_interval; ++ price_iter){
			margin_quan_supply = bids.filter_supply(price_iter);

			if(min_supply_gap > margin_quan_supply){
				min_supply_gap -= margin_quan_supply;
				max_supply_gap -= margin_quan_supply;
			}
			else{
				max_supply_gap -= min_supply_gap;
				margin_quan_supply -= min_supply_gap;
				margin_ID_supply = price_iter ;
				break;
			}
		}
		for(int price_iter = margin_ID_supply; price_iter < price_interval; ++ price_iter){
			double redispatch_price = abs(original_price - Power_market_inform.price_map.bidded_price(price_iter));
			if(reduced_flag_supply){
				redispatch_price = std::min(redispatch_price, redispatch_price_max);
			}

			if(price_iter > margin_ID_supply){
				margin_quan_supply = bids.filter_supply(price_iter);
			}

			if(max_supply_gap > margin_quan_supply){
				max_supply_gap -= margin_quan_supply;
				Power_market_inform.TSO_Market.redispatch.supply_down(tick, node_ID) += reduced_flag_supply * margin_quan_supply;
				Power_market_inform.TSO_Market.redispatch.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * margin_quan_supply;
				settlement.volume_supply.redispatch +=  (1 - 2 * reduced_flag_supply) * margin_quan_supply;
				Power_market_inform.TSO_Market.redispatch.cost_supply(tick, node_ID) += redispatch_price * margin_quan_supply;
				settlement.utility.redispatch += redispatch_price * margin_quan_supply;
			}
			else{
				Power_market_inform.TSO_Market.redispatch.supply_down(tick, node_ID) += reduced_flag_supply * max_supply_gap;
				Power_market_inform.TSO_Market.redispatch.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * max_supply_gap;
				settlement.volume_supply.redispatch +=  (1 - 2 * reduced_flag_supply) * max_supply_gap;
				Power_market_inform.TSO_Market.redispatch.cost_supply(tick, node_ID) += redispatch_price * max_supply_gap;
				settlement.utility.redispatch += redispatch_price * max_supply_gap;
				break;
			}
		}

		// Demand side
		double cleared_demand_gap = results.cleared_demand;
		double confirmed_demand_gap = results.confirmed_demand;
		double min_demand_gap = std::min(cleared_demand_gap, confirmed_demand_gap);
		double max_demand_gap = std::max(cleared_demand_gap, confirmed_demand_gap);
		bool reduced_flag_demand = (confirmed_demand_gap == min_demand_gap);
		double margin_quan_demand;
		int margin_ID_demand;
		for(int price_iter = price_interval + 1; price_iter >= 0; -- price_iter){
			margin_quan_demand = bids.filter_demand(price_iter);

			if(min_demand_gap > margin_quan_demand){
				min_demand_gap -= margin_quan_demand;
				max_demand_gap -= margin_quan_demand;
			}
			else{
				max_demand_gap -= min_demand_gap;
				margin_quan_demand -= min_demand_gap;
				margin_ID_demand = price_iter ;
				break;
			}
		}
		for(int price_iter = margin_ID_demand; price_iter >= 0; -- price_iter){
			double redispatch_price = abs(original_price - Power_market_inform.price_map.bidded_price(price_iter));
			if(reduced_flag_demand){
				redispatch_price = std::min(redispatch_price, redispatch_price_max);
			}

			if(price_iter > margin_ID_demand){
				margin_quan_demand = bids.filter_demand(price_iter);
			}

			if(max_demand_gap > margin_quan_demand){
				max_demand_gap -= margin_quan_demand;
				Power_market_inform.TSO_Market.redispatch.demand_down(tick, node_ID) += reduced_flag_demand * margin_quan_demand;
				Power_market_inform.TSO_Market.redispatch.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * margin_quan_demand;
				settlement.volume_demand.redispatch +=  (1 - 2 * reduced_flag_demand) * margin_quan_demand;
				Power_market_inform.TSO_Market.redispatch.cost_demand(tick, node_ID) += redispatch_price * margin_quan_demand;
				settlement.utility.redispatch += redispatch_price * margin_quan_demand;
			}
			else{
				Power_market_inform.TSO_Market.redispatch.demand_down(tick, node_ID) += reduced_flag_demand * max_demand_gap;
				Power_market_inform.TSO_Market.redispatch.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * max_demand_gap;
				settlement.volume_demand.redispatch +=  (1 - 2 * reduced_flag_demand) * max_demand_gap;
				Power_market_inform.TSO_Market.redispatch.cost_demand(tick, node_ID) += redispatch_price * max_demand_gap;
				settlement.utility.redispatch += redispatch_price * max_demand_gap;
				break;
			}
		}
	}

	agent::aggregator::profiles aggregator_set(power_market::market_inform &International_Market, power_network::network_inform &Power_network_inform){
		int foresight_time = agent::aggregator::parameters::foresight_time();
		int point_num = Power_network_inform.points.bidding_zone.size();

		agent::aggregator::profiles aggregator_profiles(point_num);
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int bz_ID = Power_network_inform.points.bidding_zone(point_iter);

			aggregator_profiles[point_iter].point_ID = point_iter;
			aggregator_profiles[point_iter].price_expected_profile = International_Market.confirmed.price.col(bz_ID).head(foresight_time);
			aggregator_profiles[point_iter].price_demand_profile = International_Market.confirmed.price.col(bz_ID).head(foresight_time);
			aggregator_profiles[point_iter].price_supply_profile = International_Market.confirmed.price.col(bz_ID).head(foresight_time);
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
				end_user_profiles[point_iter][sample_iter].operation.PV_scale = .01;
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
				end_user_profiles[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand = Eigen::VectorXd::Zero(foresight_time + load_shift_time_temp);
				for(int tick = load_shift_time_temp; tick < foresight_time + load_shift_time_temp; ++ tick){
					end_user_profiles[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand(tick) = Power_network_inform.points.nominal_mean_demand_field(point_iter, tick - load_shift_time_temp);
					end_user_profiles[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand(tick) *= end_user_profiles[point_iter][sample_iter].investment.decision.smart_appliance * end_user_profiles[point_iter][sample_iter].operation.smart_appliance.scale;
				}
				end_user_profiles[point_iter][sample_iter].operation.default_demand_profile *= 1. - end_user_profiles[point_iter][sample_iter].investment.decision.smart_appliance * end_user_profiles[point_iter][sample_iter].operation.smart_appliance.scale;
				end_user_profiles[point_iter][sample_iter].operation.default_PV_profile = Power_network_inform.points.solar_cf.row(point_iter).head(foresight_time);
				end_user_profiles[point_iter][sample_iter].operation.default_PV_profile *= end_user_profiles[point_iter][sample_iter].operation.PV_scale;
				end_user_profiles[point_iter][sample_iter].operation.price_demand_profile = Power_market_inform.agent_profiles.aggregators[point_iter].price_demand_profile;
				end_user_profiles[point_iter][sample_iter].operation.price_supply_profile = Power_market_inform.agent_profiles.aggregators[point_iter].price_supply_profile;

				// Set the LP problem
				agent::end_user::end_user_LP_set(end_user_profiles[point_iter][sample_iter]);

				// Set bids and results information
				agent_bids_initialization(end_user_profiles[point_iter][sample_iter].operation.bids);
				agent_results_set(end_user_profiles[point_iter][sample_iter].operation.results);
				agent_settlement_set(end_user_profiles[point_iter][sample_iter].operation.settlement);

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
			agent_results_set(profile_temp.results);
			agent_settlement_set(profile_temp.settlement);
			profile_temp.bids.submitted_demand_flex(price_interval + 1) = bid_inflex_industrial;
			profile_temp.bids.submitted_demand_flex.segment(1, price_interval) = Eigen::VectorXd::Constant(price_interval, bid_flex_industrial);
			profile_temp.bids.redispatch_demand = profile_temp.bids.submitted_demand_flex;
			profile_temp.bids.filter_demand = profile_temp.bids.submitted_demand_flex;
			profile_temp.bids.balancing_demand = profile_temp.bids.submitted_demand_flex;
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
			bid_vec /= bid_vec.sum();
			bid_vec *= Power_network_inform.plants.hydro.cap(agent_iter);

			if(Power_network_inform.plants.hydro.type(agent_iter) < 4){
				agent::power_supplier::plant_profile profile_temp;
				profile_temp.point_ID = point_ID;
				profile_temp.cap = Power_network_inform.plants.hydro.cap(agent_iter);

				// Set bids information
				agent_bids_initialization(profile_temp.bids);
				agent_results_set(profile_temp.results);
				agent_settlement_set(profile_temp.settlement);
				profile_temp.bids.submitted_supply_flex = bid_vec;
				profile_temp.bids.redispatch_supply = bid_vec;
				profile_temp.bids.filter_supply = bid_vec;
				profile_temp.bids.balancing_supply = bid_vec;

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
				agent_results_set(profile_temp.results);
				agent_settlement_set(profile_temp.settlement);
				profile_temp.bids.submitted_supply_flex = bid_vec;
				profile_temp.bids.redispatch_supply = bid_vec;
				profile_temp.bids.filter_supply = bid_vec;
				profile_temp.bids.balancing_supply = bid_vec;

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
			agent_results_set(profile_temp.results);
			agent_settlement_set(profile_temp.settlement);
			profile_temp.bids.submitted_supply_flex(price_supply_flex_ID) = bid_quan;
			profile_temp.bids.redispatch_supply = profile_temp.bids.submitted_supply_flex;
			profile_temp.bids.filter_supply = profile_temp.bids.submitted_supply_flex;
			profile_temp.bids.balancing_supply = profile_temp.bids.submitted_supply_flex;

			// High voltage power plants connect directly to transmission network
			if(profile_temp.cap >= cutoff_power){
				power_supplier_profiles.wind.HV_plant.push_back(profile_temp);
			}
			// Low voltage power plants feed into distribution network
			else{
				power_supplier_profiles.wind.LV_plant.push_back(profile_temp);
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
			double marginal_price = Power_market_inform.International_Market.confirmed.price(tick, bz_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				// Inflexible bids have first / last priority for dispatch
				if(marginal_price_ID > 0){
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_supply += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_inflex.head(marginal_price_ID).sum();
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_supply += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_flex.head(marginal_price_ID).sum();
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(0) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_inflex.head(marginal_price_ID).sum();
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(0) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_inflex.head(marginal_price_ID).sum();
				}
				if(marginal_price_ID < price_interval + 1){
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_demand += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_inflex.tail(price_interval + 1 - marginal_price_ID).sum();
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_demand += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_flex.tail(price_interval + 1 - marginal_price_ID).sum();
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(price_interval + 1) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_inflex.tail(price_interval + 1 - marginal_price_ID).sum();
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(price_interval + 1) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_inflex.tail(price_interval + 1 - marginal_price_ID).sum();
				}
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_demand += Power_market_inform.International_Market.confirmed.ratio_demand(bz_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_inflex(marginal_price_ID);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_demand += Power_market_inform.International_Market.confirmed.ratio_demand(bz_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_flex(marginal_price_ID);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_supply += Power_market_inform.International_Market.confirmed.ratio_supply(bz_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_inflex(marginal_price_ID);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_supply += Power_market_inform.International_Market.confirmed.ratio_supply(bz_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_flex(marginal_price_ID);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(0) += (1. - Power_market_inform.International_Market.confirmed.ratio_demand(bz_ID)) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_inflex(marginal_price_ID);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(0) += Power_market_inform.International_Market.confirmed.ratio_supply(bz_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_inflex(marginal_price_ID);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(price_interval + 1) += Power_market_inform.International_Market.confirmed.ratio_demand(bz_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_inflex(marginal_price_ID);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(price_interval + 1) += (1. - Power_market_inform.International_Market.confirmed.ratio_supply(bz_ID)) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_inflex(marginal_price_ID);

				// Flexible bids have redispatcch priority in between
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_flex;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_flex;

				// Update filter bids
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.filter_demand = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.filter_supply = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply;
			}
		}
	}

	void end_user_filter_demand_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int point_num = Power_network_inform.points.bidding_zone.size();
		int sample_num = agent::end_user::parameters::sample_num();
		int price_interval = power_market::parameters::price_interval();

		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int point_ID = Power_network_inform.points.in_cluster_ID(point_iter);
			int node_ID = Power_network_inform.points.node(point_iter);
			int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			double marginal_price = Power_market_inform.DSO_Markets[DSO_ID].confirmed.price(tick, point_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				//Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.filter_demand = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand;
				if(marginal_price_ID > 0){
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand.head(marginal_price_ID) *= 0.;
				}
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(marginal_price_ID) *= Power_market_inform.DSO_Markets[DSO_ID].confirmed.ratio_demand(point_ID);
			}
		}
	}

	void power_supplier_filter_demand_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int price_interval = power_market::parameters::price_interval();

		int pump_LV_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
		for(int agent_iter = 0; agent_iter < pump_LV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			point_ID = Power_network_inform.points.in_cluster_ID(point_ID);
			double marginal_price = Power_market_inform.DSO_Markets[DSO_ID].confirmed.price(tick, point_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			//Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.filter_demand = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_demand;
			if(marginal_price_ID > 0){
				Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_demand.head(marginal_price_ID) *= 0.;
			}
			Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_demand(marginal_price_ID) *= Power_market_inform.DSO_Markets[DSO_ID].confirmed.ratio_demand(point_ID);
			Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.balancing_demand = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_demand;
		}
	}

	void end_user_filter_supply_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int point_num = Power_network_inform.points.bidding_zone.size();
		int sample_num = agent::end_user::parameters::sample_num();
		int price_interval = power_market::parameters::price_interval();

		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int point_ID = Power_network_inform.points.in_cluster_ID(point_iter);
			int node_ID = Power_network_inform.points.node(point_iter);
			int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			double marginal_price = Power_market_inform.DSO_Markets[DSO_ID].confirmed.price(tick, point_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				//Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.filter_supply = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply;
				if(marginal_price_ID < price_interval + 1){
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply.tail(price_interval + 1 - marginal_price_ID) *= 0.;
				}
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(marginal_price_ID) *= Power_market_inform.DSO_Markets[DSO_ID].confirmed.ratio_supply(point_ID);
			}
		}
	}

	void power_supplier_filter_supply_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int price_interval = power_market::parameters::price_interval();

		int hydro_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant.size();
		for(int agent_iter = 0; agent_iter < hydro_LV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			point_ID = Power_network_inform.points.in_cluster_ID(point_ID);
			double marginal_price = Power_market_inform.DSO_Markets[DSO_ID].confirmed.price(tick, point_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			//Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids.filter_supply = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids.redispatch_supply;
			if(marginal_price_ID < price_interval + 1){
				Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids.redispatch_supply.tail(price_interval + 1 - marginal_price_ID) *= 0.;
			}
			Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids.redispatch_supply(marginal_price_ID) *= Power_market_inform.DSO_Markets[DSO_ID].confirmed.ratio_supply(point_ID);
			Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids.balancing_supply = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids.redispatch_supply;
		}

		int wind_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant.size();
		for(int agent_iter = 0; agent_iter < wind_LV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			point_ID = Power_network_inform.points.in_cluster_ID(point_ID);
			double marginal_price = Power_market_inform.DSO_Markets[DSO_ID].confirmed.price(tick, point_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			//Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.filter_supply = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.redispatch_supply;
			if(marginal_price_ID < price_interval + 1){
				Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.redispatch_supply.tail(price_interval + 1 - marginal_price_ID) *= 0.;
			}
			Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.redispatch_supply(marginal_price_ID) *= Power_market_inform.DSO_Markets[DSO_ID].confirmed.ratio_supply(point_ID);
			Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.balancing_supply = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.redispatch_supply;
		}

		int pump_LV_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
		for(int agent_iter = 0; agent_iter < pump_LV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			point_ID = Power_network_inform.points.in_cluster_ID(point_ID);
			double marginal_price = Power_market_inform.DSO_Markets[DSO_ID].confirmed.price(tick, point_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			//Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.filter_supply = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_supply;
			if(marginal_price_ID < price_interval + 1){
				Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_supply.tail(price_interval + 1 - marginal_price_ID) *= 0.;
			}
			Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_supply(marginal_price_ID) *= Power_market_inform.DSO_Markets[DSO_ID].confirmed.ratio_supply(point_ID);
			Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.balancing_supply = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_supply;
		}
	}

	void end_user_balancing_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int point_num = Power_network_inform.points.bidding_zone.size();
		int sample_num = agent::end_user::parameters::sample_num();
		int price_interval = power_market::parameters::price_interval();
		double redispatch_price_max = power_market::parameters::redispatch_price_max();

		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int node_ID = Power_network_inform.points.node(point_iter);
			int bz_ID = Power_network_inform.points.bidding_zone(point_iter);
			double marginal_price = Power_market_inform.TSO_Market.confirmed.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];
			double original_price = Power_market_inform.International_Market.confirmed.price(tick, bz_ID);
			int original_price_ID = Power_market_inform.price_map.price_ID[original_price];

			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				if(marginal_price_ID > 0){
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.confirmed_supply += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply.head(marginal_price_ID).sum();
					if(!Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].investment.decision.control_reserve){
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand(0) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand.head(marginal_price_ID).sum();
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply(0) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply.head(marginal_price_ID).sum();
					}
				}
				if(marginal_price_ID < price_interval + 1){
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.confirmed_demand += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand.tail(price_interval + 1 - marginal_price_ID).sum();
					if(!Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].investment.decision.control_reserve){
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand(price_interval + 1) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand.tail(price_interval + 1 - marginal_price_ID).sum();
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply(price_interval + 1) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply.tail(price_interval + 1 - marginal_price_ID).sum();
					}
				}
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.confirmed_supply += Power_market_inform.TSO_Market.confirmed.ratio_supply(node_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(marginal_price_ID);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.confirmed_demand += Power_market_inform.TSO_Market.confirmed.ratio_demand(node_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(marginal_price_ID);
				if(!Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].investment.decision.control_reserve){
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand(0) += (1. - Power_market_inform.TSO_Market.confirmed.ratio_demand(node_ID)) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(marginal_price_ID);
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply(0) += Power_market_inform.TSO_Market.confirmed.ratio_supply(node_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(marginal_price_ID);
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand(price_interval + 1) += Power_market_inform.TSO_Market.confirmed.ratio_demand(node_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(marginal_price_ID);
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply(price_interval + 1) += (1. - Power_market_inform.TSO_Market.confirmed.ratio_supply(node_ID)) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(marginal_price_ID);
				}

				if(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].investment.decision.control_reserve){
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand;
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply;
				}

				// Imbalance accounting
				double imbalance = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.confirmed_demand;
				imbalance *= Power_network_inform.points.imbalance_field(point_iter, tick);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand(price_interval + 1) += (imbalance >= 0.) * imbalance;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply(0) -= (imbalance < 0.) * imbalance;

				// Settlement in EOM
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement.volume_supply.EOM += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.confirmed_supply;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement.volume_demand.EOM += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.confirmed_demand;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement.price.EOM += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.confirmed_demand * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.price_demand_profile(0);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement.utility.EOM += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.confirmed_supply * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.price_supply_profile(0);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement.utility.EOM += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.confirmed_demand * Power_market_inform.price_map.bidded_price(price_interval + 1);

				// Settlement of redispatch
				agent_redispatch_settlement_calculation(tick, node_ID, original_price, Power_market_inform, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement);
			}
		}
	}

	void industrial_balancing_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int price_interval = power_market::parameters::price_interval();
		double redispatch_price_max = power_market::parameters::redispatch_price_max();

		int industrial_HV_num = Power_market_inform.agent_profiles.industrial.HV.size();
		for(int agent_iter = 0; agent_iter < industrial_HV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.industrial.HV[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
			double marginal_price = Power_market_inform.TSO_Market.confirmed.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];
			double original_price = Power_market_inform.International_Market.confirmed.price(tick, bz_ID);
			int original_price_ID = Power_market_inform.price_map.price_ID[original_price];

			// Calculate scheduled results
			agent_scheduled_results_calculation(bz_ID, node_ID, marginal_price_ID, original_price_ID, Power_market_inform, Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids, Power_market_inform.agent_profiles.industrial.HV[agent_iter].results);

			double imbalance = Power_market_inform.agent_profiles.industrial.HV[agent_iter].results.actual_demand;
			imbalance *= Power_network_inform.points.imbalance_field(point_ID, tick);
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids.balancing_demand(price_interval + 1) += (imbalance >= 0.) * imbalance;
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids.balancing_supply(0) -= (imbalance < 0.) * imbalance;

			// Settlement in EOM
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement.volume_supply.EOM += Power_market_inform.agent_profiles.industrial.HV[agent_iter].results.confirmed_supply;
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement.volume_demand.EOM += Power_market_inform.agent_profiles.industrial.HV[agent_iter].results.confirmed_demand;
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement.price.EOM += Power_market_inform.agent_profiles.industrial.HV[agent_iter].results.confirmed_demand * original_price;
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement.utility.EOM += Power_market_inform.agent_profiles.industrial.HV[agent_iter].results.confirmed_supply * original_price;
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement.utility.EOM += Power_market_inform.agent_profiles.industrial.HV[agent_iter].results.confirmed_demand * Power_market_inform.price_map.bidded_price(price_interval + 1);

			// Settlement of redispatch
			agent_redispatch_settlement_calculation(tick, node_ID, original_price, Power_market_inform, Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids, Power_market_inform.agent_profiles.industrial.HV[agent_iter].results, Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement);
		}
	}

	void power_supplier_balancing_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int price_interval = power_market::parameters::price_interval();
		double redispatch_price_max = power_market::parameters::redispatch_price_max();

		int hydro_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant.size();
		for(int agent_iter = 0; agent_iter < hydro_HV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
			double marginal_price = Power_market_inform.TSO_Market.confirmed.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];
			double original_price = Power_market_inform.International_Market.confirmed.price(tick, bz_ID);
			int original_price_ID = Power_market_inform.price_map.price_ID[original_price];

			// Calculate scheduled results
			agent_scheduled_results_calculation(bz_ID, node_ID, marginal_price_ID, original_price_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].results);

			// Settlement in EOM
			Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].settlement.volume_supply.EOM += Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].results.confirmed_supply;
			Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].settlement.volume_demand.EOM += Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].results.confirmed_demand;
			Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].settlement.price.EOM += Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].results.confirmed_demand * original_price;
			Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].settlement.utility.EOM += Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].results.confirmed_supply * original_price;
			Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].settlement.utility.EOM += Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].results.confirmed_demand * Power_market_inform.price_map.bidded_price(price_interval + 1);

			// Settlement of redispatch
			agent_redispatch_settlement_calculation(tick, node_ID, original_price, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].settlement);
		}

		int hydro_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant.size();
		for(int agent_iter = 0; agent_iter < hydro_LV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
			double marginal_price = Power_market_inform.TSO_Market.confirmed.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];
			double original_price = Power_market_inform.International_Market.confirmed.price(tick, bz_ID);
			int original_price_ID = Power_market_inform.price_map.price_ID[original_price];

			// Calculate scheduled results
			agent_scheduled_results_calculation(bz_ID, node_ID, marginal_price_ID, original_price_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].results);

			// Settlement in EOM
			Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].settlement.volume_supply.EOM += Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].results.confirmed_supply;
			Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].settlement.volume_demand.EOM += Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].results.confirmed_demand;
			Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].settlement.price.EOM += Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].results.confirmed_demand * original_price;
			Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].settlement.utility.EOM += Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].results.confirmed_supply * original_price;
			Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].settlement.utility.EOM += Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].results.confirmed_demand * Power_market_inform.price_map.bidded_price(price_interval + 1);

			// Settlement of redispatch
			agent_redispatch_settlement_calculation(tick, node_ID, original_price, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].settlement);
		}

		int wind_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant.size();
		for(int agent_iter = 0; agent_iter < wind_HV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
			double marginal_price = Power_market_inform.TSO_Market.confirmed.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];
			double original_price = Power_market_inform.International_Market.confirmed.price(tick, bz_ID);
			int original_price_ID = Power_market_inform.price_map.price_ID[original_price];

			// Calculate scheduled results
			agent_scheduled_results_calculation(bz_ID, node_ID, marginal_price_ID, original_price_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].results);

			// Settlement in EOM
			Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].settlement.volume_supply.EOM += Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].results.confirmed_supply;
			Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].settlement.volume_demand.EOM += Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].results.confirmed_demand;
			Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].settlement.price.EOM += Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].results.confirmed_demand * original_price;
			Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].settlement.utility.EOM += Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].results.confirmed_supply * original_price;
			Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].settlement.utility.EOM += Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].results.confirmed_demand * Power_market_inform.price_map.bidded_price(price_interval + 1);

			// Settlement of redispatch
			agent_redispatch_settlement_calculation(tick, node_ID, original_price, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].settlement);
		}

		int wind_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant.size();
		for(int agent_iter = 0; agent_iter < wind_LV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
			double marginal_price = Power_market_inform.TSO_Market.confirmed.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];
			double original_price = Power_market_inform.International_Market.confirmed.price(tick, bz_ID);
			int original_price_ID = Power_market_inform.price_map.price_ID[original_price];

			// Calculate scheduled results
			agent_scheduled_results_calculation(bz_ID, node_ID, marginal_price_ID, original_price_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results);

			// Settlement in EOM
			Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].settlement.volume_supply.EOM += Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results.confirmed_supply;
			Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].settlement.volume_demand.EOM += Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results.confirmed_demand;
			Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].settlement.price.EOM += Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results.confirmed_demand * original_price;
			Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].settlement.utility.EOM += Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results.confirmed_supply * original_price;
			Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].settlement.utility.EOM += Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results.confirmed_demand * Power_market_inform.price_map.bidded_price(price_interval + 1);

			// Settlement of redispatch
			agent_redispatch_settlement_calculation(tick, node_ID, original_price, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].settlement);
		}
	}

	void agents_redispatch_settlement(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int node_num = Power_market_inform.TSO_Market.network.num_vertice;
		Eigen::VectorXd redispatch_price_demand(node_num);

		// Redispatch price per energy at each transmission node
		for(int node_iter = 0; node_iter < node_num; ++ node_iter){
			double redispatched_qaun_demand = Power_market_inform.TSO_Market.confirmed.demand(tick, node_iter);
			redispatched_qaun_demand -= Power_market_inform.TSO_Market.redispatch.demand_up(tick, node_iter);
			redispatch_price_demand(node_iter) = Power_market_inform.TSO_Market.redispatch.cost_demand(tick, node_iter);
			redispatch_price_demand(node_iter) /= redispatched_qaun_demand;
			std::cout << Power_market_inform.TSO_Market.redispatch.cost_demand(tick, node_iter) << "\t";
			std::cout << redispatched_qaun_demand + Power_market_inform.TSO_Market.redispatch.demand_down(tick, node_iter) << "\t";
			std::cout << Power_market_inform.TSO_Market.confirmed.demand(tick, node_iter) << "\t" << redispatch_price_demand(node_iter) << "\n";
		}
		//std::cout << redispatch_price_demand;
		std::cout << "\n";
	}

	void end_user_status_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, bool control_reserve_flag){
		int point_num = Power_network_inform.points.bidding_zone.size();
		int sample_num = agent::end_user::parameters::sample_num();
		int load_shift_time = agent::end_user::parameters::load_shift_time();
		int price_interval = power_market::parameters::price_interval();

		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int node_ID = Power_network_inform.points.node(point_iter);
			double marginal_price = Power_market_inform.TSO_Market.actual.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				if(control_reserve_flag){
					if(marginal_price_ID > 0){
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_supply += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply.head(marginal_price_ID).sum();
					}
					if(marginal_price_ID < price_interval + 1){
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_demand += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand.tail(price_interval + 1 - marginal_price_ID).sum();
					}
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_supply += Power_market_inform.TSO_Market.actual.ratio_supply(node_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply(marginal_price_ID);
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_demand += Power_market_inform.TSO_Market.actual.ratio_demand(node_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand(marginal_price_ID);
				}
				else{
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_supply = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.confirmed_supply;
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_demand = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.confirmed_demand;
				}

				double gap = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_demand;
				gap -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_supply;
				gap -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_inflex.sum();
				gap -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_flex.sum();
				gap += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_inflex.sum();
				gap += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_flex.sum();

				// Actual demand greater than initially planned
				while(gap > 0.){
					// Reduce BESS charge
					double BESS_flex_lb = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.soc - Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.self_consumption;
					BESS_flex_lb = -std::min(BESS_flex_lb, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.capacity_scale);
					double BESS_flex = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity;
					if(BESS_flex > 0.){
						BESS_flex = std::max(0., BESS_flex - gap * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency);
						BESS_flex = std::max(BESS_flex, BESS_flex_lb);
						gap -= (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity - BESS_flex) / Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency;
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity = BESS_flex;
						if(gap == 0.){
							break;
						}
					}
					if(BESS_flex > BESS_flex_lb){
						BESS_flex = std::max(BESS_flex_lb, BESS_flex - gap / Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency);
						gap -= (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity - BESS_flex) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency;
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity = BESS_flex;
						if(gap == 0.){
							break;
						}
					}

					// Reduce EV charge
					double EV_flex_lb = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc - Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.self_consumption;
					EV_flex_lb -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.default_demand_profile(0);
					EV_flex_lb = -std::min(EV_flex_lb, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.capacity_scale);
					double EV_flex = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity;
					if(EV_flex > 0.){
						EV_flex = std::max(0., EV_flex - gap * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency);
						EV_flex = std::max(EV_flex, EV_flex_lb);
						gap -= (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity - EV_flex) / Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency;
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity = EV_flex;
						if(gap == 0.){
							break;
						}
					}
					if(EV_flex > EV_flex_lb){
						EV_flex = std::max(EV_flex_lb, EV_flex - gap / Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency);
						gap -= (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity - EV_flex) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency;
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity = EV_flex;
						if(gap == 0.){
							break;
						}
					}

					// Reduce smart appliance demand
					for(int tock = 0; tock < 2 * load_shift_time + 1; ++ tock){
						int tock_ID = 2 * load_shift_time - tock;
						double sa_flex = std::max(0., Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.scheduled_demand(tock_ID) - gap);
						gap -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.scheduled_demand(tock_ID) - sa_flex;
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.scheduled_demand(tock_ID) = sa_flex;
						if(gap == 0.){
							break;
						}
					}

					break;
				}

				// Actual demand smaller than initially planned
				while(gap < 0.){
					// Increase BESS charge
					double BESS_flex_ub = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.energy_scale;
					BESS_flex_ub -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.soc;
					BESS_flex_ub += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.self_consumption;
					BESS_flex_ub = std::min(BESS_flex_ub, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.capacity_scale);
					double BESS_flex = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity;
					if(BESS_flex < 0.){
						BESS_flex = std::min(0., BESS_flex - gap / Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency);
						gap -= (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity - BESS_flex) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency;
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity = BESS_flex;
						if(gap == 0.){
							break;
						}
					}
					BESS_flex = std::min(BESS_flex_ub, BESS_flex - gap * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency);
					gap -= (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity - BESS_flex) / Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency;
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity = BESS_flex;
					if(gap == 0.){
						break;
					}

					// Increase EV charge
					double EV_flex_ub = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.energy_scale;
					EV_flex_ub -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc;
					EV_flex_ub += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.self_consumption;
					EV_flex_ub += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.default_demand_profile(0);
					EV_flex_ub = std::min(EV_flex_ub, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.capacity_scale);
					double EV_flex = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity;
					if(EV_flex < 0.){
						EV_flex = std::min(0., EV_flex - gap / Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency);
						gap -= (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity - EV_flex) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency;
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity = EV_flex;
						if(gap == 0.){
							break;
						}
					}
					EV_flex = std::min(EV_flex_ub, EV_flex - gap * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency);
					gap -= (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity - EV_flex) / Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency;
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity = EV_flex;

					break;
				}

				// Update state variables of end-users
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.soc += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.soc -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.self_consumption;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.self_consumption;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.default_demand_profile(0);
				for(int tock = 0; tock < 2 * load_shift_time + 1; ++ tock){
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand(tock) -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.scheduled_demand(tock);
				}

				// Balancing settlement
			}
		}
	}

	void aggregator_price_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int foresight_time = agent::aggregator::parameters::foresight_time();
		int aggregator_num = Power_market_inform.agent_profiles.aggregators.size();

		for(int agent_iter = 0; agent_iter < aggregator_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.aggregators[agent_iter].point_ID;
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
			Eigen::VectorXd bid_vec = Power_market_inform.International_Market.merit_order_curve.col(bz_ID);
			bid_vec *= Power_network_inform.plants.hydro.cap(agent_iter);
			bid_vec /= (Power_market_inform.International_Market.merit_order_curve.col(bz_ID).sum());

			Power_market_inform.agent_profiles.aggregators[agent_iter].price_expected_profile = Power_market_inform.International_Market.confirmed.price.col(bz_ID).segment(tick, foresight_time);
			Power_market_inform.agent_profiles.aggregators[agent_iter].price_demand_profile = Power_market_inform.International_Market.confirmed.price.col(bz_ID).segment(tick, foresight_time);
			Power_market_inform.agent_profiles.aggregators[agent_iter].price_supply_profile = Power_market_inform.International_Market.confirmed.price.col(bz_ID).segment(tick, foresight_time);
		}
	}

	void end_user_submit_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int foresight_time = agent::end_user::parameters::foresight_time();
		int point_num = Power_network_inform.points.bidding_zone.size();
		int sample_num = agent::end_user::parameters::sample_num();
		int price_interval = power_market::parameters::price_interval();
		double residential_ratio = agent::parameters::residential_ratio();

		// Update of forecast demand profile and operation strategies
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				int load_shift_time = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.shift_time;

				// Update of input profiles
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.house_default_period = Eigen::VectorXi::Ones(foresight_time);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.default_demand_profile = Eigen::VectorXd::Zero(foresight_time);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.default_demand_profile *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].investment.decision.EV_self_charging;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.default_demand_profile = Power_network_inform.points.nominal_mean_demand_field.row(point_iter).segment(tick, foresight_time);
				//Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.default_demand_profile -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.default_demand_profile;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand.head(foresight_time + load_shift_time - 1) = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand.tail(foresight_time + load_shift_time - 1);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand(foresight_time + load_shift_time - 1) = Power_network_inform.points.nominal_mean_demand_field(point_iter, tick + foresight_time - 1);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand(foresight_time + load_shift_time - 1) *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].investment.decision.smart_appliance * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.scale;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.default_demand_profile *= 1. - Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].investment.decision.smart_appliance * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.scale;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.default_PV_profile = Power_network_inform.points.solar_cf.row(point_iter).segment(tick, foresight_time);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.default_PV_profile *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.PV_scale;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.price_demand_profile = Power_market_inform.agent_profiles.aggregators[point_iter].price_demand_profile;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.price_supply_profile = Power_market_inform.agent_profiles.aggregators[point_iter].price_supply_profile;

				// Set bids and results information
				agent_bids_initialization(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids);
				agent_results_set(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results);

				// Optimization and update process variables
				agent::end_user::end_user_LP_optimize(tick, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter]);

				// Scale the bids correctly
				double scale = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.weight;
				scale *= agent::parameters::residential_ratio();
				scale *= Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
				agent_submitted_bids_scale(scale, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids);
			}
		}
	}

	void industrial_submit_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int industrial_HV_num = Power_market_inform.agent_profiles.industrial.HV.size();
		int price_interval = power_market::parameters::price_interval();

		for(int agent_iter = 0; agent_iter < industrial_HV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.industrial.HV[agent_iter].point_ID;

			double bid_inflex_industrial = Power_network_inform.points.nominal_mean_demand_field(point_ID, tick);
			bid_inflex_industrial *= Power_network_inform.points.population_density(point_ID) * Power_network_inform.points.point_area / 1000.;
			bid_inflex_industrial *= 1. - agent::parameters::residential_ratio();
			double bid_flex_industrial = bid_inflex_industrial;
			bid_inflex_industrial *= 1. - agent::industrial::flexible_ratio();
			bid_flex_industrial *= agent::industrial::flexible_ratio();
			bid_flex_industrial /= price_interval;

			// Set bids information
			agent_bids_initialization(Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids);
			agent_results_set(Power_market_inform.agent_profiles.industrial.HV[agent_iter].results);
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids.submitted_demand_flex(price_interval + 1) = bid_inflex_industrial;
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids.submitted_demand_flex.segment(1, price_interval) = Eigen::VectorXd::Constant(price_interval, bid_flex_industrial);
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids.redispatch_demand = Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids.submitted_demand_flex;
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids.filter_demand = Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids.submitted_demand_flex;
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids.balancing_demand = Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids.submitted_demand_flex;
		}
	}

	void power_supplier_submit_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int price_interval = power_market::parameters::price_interval();

		int hydro_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant.size();
		for(int agent_iter = 0; agent_iter < hydro_HV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].point_ID;
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
			Eigen::VectorXd bid_vec = Power_market_inform.International_Market.merit_order_curve.col(bz_ID);
			bid_vec *= Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].cap;
			bid_vec /= (Power_market_inform.International_Market.merit_order_curve.col(bz_ID).sum());

			// Set bids information
			agent_bids_initialization(Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].bids);
			agent_results_set(Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].results);
			Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].bids.submitted_supply_flex = bid_vec;
			Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].bids.redispatch_supply = bid_vec;
			Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].bids.filter_supply = bid_vec;
			Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].bids.balancing_supply = bid_vec;
		}

		int hydro_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant.size();
		for(int agent_iter = 0; agent_iter < hydro_LV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].point_ID;
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
			Eigen::VectorXd bid_vec = Power_market_inform.International_Market.merit_order_curve.col(bz_ID);
			bid_vec *= Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].cap;
			bid_vec /= (Power_market_inform.International_Market.merit_order_curve.col(bz_ID).sum());

			// Set bids information
			agent_bids_initialization(Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids);
			agent_results_set(Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].results);
			Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids.submitted_supply_flex = bid_vec;
			Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids.redispatch_supply = bid_vec;
			Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids.filter_supply = bid_vec;
			Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids.balancing_supply = bid_vec;
		}

		int wind_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant.size();
		for(int agent_iter = 0; agent_iter < wind_HV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].point_ID;
			int price_supply_flex_ID = 0;
			double bid_quan = Power_network_inform.points.wind_on_cf(point_ID, 0) * Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].cap;

			// Set bids information
			agent_bids_initialization(Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids);
			agent_results_set(Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].results);
			Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids.submitted_supply_flex(price_supply_flex_ID) = bid_quan;
			Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids.redispatch_supply(price_supply_flex_ID) = bid_quan;
			Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids.filter_supply(price_supply_flex_ID) = bid_quan;
			Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids.balancing_supply(price_supply_flex_ID) = bid_quan;
		}

		int wind_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant.size();
		for(int agent_iter = 0; agent_iter < wind_LV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].point_ID;
			int price_supply_flex_ID = 0;
			double bid_quan = Power_network_inform.points.wind_on_cf(point_ID, 0) * Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].cap;

			// Set bids information
			agent_bids_initialization(Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids);
			agent_results_set(Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results);
			Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.submitted_supply_flex(price_supply_flex_ID) = bid_quan;
			Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.redispatch_supply(price_supply_flex_ID) = bid_quan;
			Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.filter_supply(price_supply_flex_ID) = bid_quan;
			Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.balancing_supply(price_supply_flex_ID) = bid_quan;
		}

		int pump_HV_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV.size();
		for(int agent_iter = 0; agent_iter < pump_HV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].point_ID;
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
			Eigen::VectorXd bid_vec = Power_market_inform.International_Market.merit_order_curve.col(bz_ID);
			bid_vec *= Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].cap;
			bid_vec /= (Power_market_inform.International_Market.merit_order_curve.col(bz_ID).sum());

			// Set bids information
			agent_bids_initialization(Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].bids);
			agent_results_set(Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].results);
			Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].bids.submitted_supply_flex = bid_vec;
			Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].bids.redispatch_supply = bid_vec;
			Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].bids.filter_supply = bid_vec;
			Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].bids.balancing_supply = bid_vec;
		}

		int pump_LV_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
		for(int agent_iter = 0; agent_iter < pump_LV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
			Eigen::VectorXd bid_vec = Power_market_inform.International_Market.merit_order_curve.col(bz_ID);
			bid_vec *= Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].cap;
			bid_vec /= (Power_market_inform.International_Market.merit_order_curve.col(bz_ID).sum());

			// Set bids information
			agent_bids_initialization(Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids);
			agent_results_set(Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].results);
			Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.submitted_supply_flex = bid_vec;
			Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_supply = bid_vec;
			Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.filter_supply = bid_vec;
			Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.balancing_supply = bid_vec;
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

void agent::agents_filter_demand_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
	end_user_filter_demand_update(tick, Power_market_inform, Power_network_inform);
	power_supplier_filter_demand_update(tick, Power_market_inform, Power_network_inform);
}

void agent::agents_filter_supply_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
	end_user_filter_supply_update(tick, Power_market_inform, Power_network_inform);
	power_supplier_filter_supply_update(tick, Power_market_inform, Power_network_inform);
}

void agent::agents_balancing_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
	end_user_balancing_update(tick, Power_market_inform, Power_network_inform);
	industrial_balancing_update(tick, Power_market_inform, Power_network_inform);
	power_supplier_balancing_update(tick, Power_market_inform, Power_network_inform);
	agents_redispatch_settlement(tick, Power_market_inform, Power_network_inform);
}

void agent::agents_status_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, bool control_reserve_flag){
	end_user_status_update(tick, Power_market_inform, Power_network_inform, control_reserve_flag);
}

void agent::agents_submit_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
	aggregator_price_update(tick, Power_market_inform, Power_network_inform);
	end_user_submit_update(tick, Power_market_inform, Power_network_inform);
	industrial_submit_update(tick, Power_market_inform, Power_network_inform);
	power_supplier_submit_update(tick, Power_market_inform, Power_network_inform);
}
