// Source file for operational strategies of different agents
#include "agent_func.h"

namespace{
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
			int bz_ID = Power_network_inform.points.bidding_zone(point_iter);

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

				// Set bids information
				agent_bids_initialization(end_user_profiles[point_iter][sample_iter].operation.bids);

				// Optimization and update process variables
				agent::end_user::end_user_LP_optimize(0, end_user_profiles[point_iter][sample_iter]);
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
		for(int hydro_iter = 0; hydro_iter < hydro_num; ++ hydro_iter){
			int x_ID = int((Power_network_inform.plants.hydro.x(hydro_iter) - Power_network_inform.points.x.minCoeff()) / Power_network_inform.points.grid_length + .5);
			int y_ID = int((Power_network_inform.plants.hydro.y(hydro_iter) - Power_network_inform.points.y.minCoeff()) / Power_network_inform.points.grid_length + .5);
			int point_ID = Power_network_inform.points.coordinate_grid(x_ID, y_ID);
			if(point_ID == -1){
				continue;
			}
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);
			Eigen::VectorXd bid_vec = Power_market_inform.International_Market.merit_order_curve.col(bz_ID);
			bid_vec *= Power_network_inform.plants.hydro.cap(hydro_iter);
			bid_vec /= (Power_market_inform.International_Market.merit_order_curve.col(bz_ID).sum());

			if(Power_network_inform.plants.hydro.type(hydro_iter) < 4){
				agent::power_supplier::plant_profile plant_temp;
				plant_temp.point_ID = point_ID;
				plant_temp.cap = Power_network_inform.plants.hydro.cap(hydro_iter);

				// Set bids information
				agent_bids_initialization(plant_temp.bids);
				plant_temp.bids.submitted_supply_flex = bid_vec;
				plant_temp.bids.redispatch_supply = bid_vec;

				// High voltage power plants connect directly to transmission network
				if(plant_temp.cap >= cutoff_power){
					power_supplier_profiles.hydro.HV_plant.push_back(plant_temp);
				}
				// Low voltage power plants feed into distribution network
				else{
					power_supplier_profiles.hydro.LV_plant.push_back(plant_temp);
				}
			}
			else if(Power_network_inform.plants.hydro.type(hydro_iter) == 5){
				agent::power_supplier::storage_profile storage_temp;
				storage_temp.point_ID = point_ID;
				storage_temp.cap = Power_network_inform.plants.hydro.cap(hydro_iter);

				// Set bids information
				agent_bids_initialization(storage_temp.bids);
				storage_temp.bids.submitted_supply_flex = bid_vec;
				storage_temp.bids.redispatch_supply = bid_vec;

				// High voltage power plants connect directly to transmission network
				if(storage_temp.cap >= cutoff_power){
					power_supplier_profiles.pump_storage.HV.push_back(storage_temp);
				}
				// Low voltage power plants feed into distribution network
				else{
					power_supplier_profiles.pump_storage.LV.push_back(storage_temp);
				}
			}
		}

		power_supplier_profiles.wind.HV_hybrid.reserve(wind_num);
		power_supplier_profiles.wind.HV_plant.reserve(wind_num);
		power_supplier_profiles.wind.LV_hybrid.reserve(wind_num);
		power_supplier_profiles.wind.LV_plant.reserve(wind_num);
		for(int wind_iter = 0; wind_iter < wind_num; ++ wind_iter){
			int x_ID = int((Power_network_inform.plants.wind.x(wind_iter) - Power_network_inform.points.x.minCoeff()) / Power_network_inform.points.grid_length + .5);
			int y_ID = int((Power_network_inform.plants.wind.y(wind_iter) - Power_network_inform.points.y.minCoeff()) / Power_network_inform.points.grid_length + .5);
			int point_ID = Power_network_inform.points.coordinate_grid(x_ID, y_ID);
			if(point_ID == -1){
				continue;
			}
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);
			int price_supply_flex_ID = 0;
			double bid_quan = Power_network_inform.points.wind_on_cf(point_ID, 0) * Power_network_inform.plants.wind.cap(wind_iter);

			agent::power_supplier::plant_profile plant_temp;
			plant_temp.point_ID = point_ID;
			plant_temp.cap = Power_network_inform.plants.wind.cap(wind_iter);

			// Set bids information
			agent_bids_initialization(plant_temp.bids);
			plant_temp.bids.submitted_supply_flex(price_supply_flex_ID) = bid_quan;
			plant_temp.bids.redispatch_supply = plant_temp.bids.submitted_supply_flex;

			// High voltage power plants connect directly to transmission network
			if(plant_temp.cap >= cutoff_power){
				power_supplier_profiles.hydro.HV_plant.push_back(plant_temp);
			}
			// Low voltage power plants feed into distribution network
			else{
				power_supplier_profiles.hydro.LV_plant.push_back(plant_temp);
			}
		}

		return power_supplier_profiles;
	}
}

void agent::agents_set(power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
	Power_market_inform.agent_profiles.aggregators = aggregator_set(Power_market_inform.International_Market, Power_network_inform);
	Power_market_inform.agent_profiles.end_users = end_user_set(Power_market_inform, Power_network_inform);
	Power_market_inform.agent_profiles.industrial = industrial_set(Power_network_inform);
	power_supplier_set(Power_market_inform, Power_network_inform);
}
