// Source file for operational strategies of different agents
#include "agent_func.h"

namespace{
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
				end_user_profiles[point_iter][sample_iter].operation.bids.submitted_supply_inflex = Eigen::VectorXd::Zero(price_interval + 2);
				end_user_profiles[point_iter][sample_iter].operation.bids.submitted_demand_inflex = Eigen::VectorXd::Zero(price_interval + 2);
				end_user_profiles[point_iter][sample_iter].operation.bids.submitted_supply_flex = Eigen::VectorXd::Zero(price_interval + 2);
				end_user_profiles[point_iter][sample_iter].operation.bids.submitted_demand_flex = Eigen::VectorXd::Zero(price_interval + 2);

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
			int node_ID = Power_network_inform.points.node(point_iter);
			industrial_profiles.HV[point_iter].node_ID = node_ID;
			industrial_profiles.HV[point_iter].bids.submitted_demand_flex = Eigen::VectorXd::Zero(price_interval + 2);

			double bid_inflex_industrial = Power_network_inform.points.nominal_mean_demand_field(point_iter, 0);
			bid_inflex_industrial *= Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
			bid_inflex_industrial *= 1. - agent::parameters::residential_ratio();
			double bid_flex_industrial = bid_inflex_industrial;
			bid_inflex_industrial *= 1. - agent::industrial::flexible_ratio();
			bid_flex_industrial *= agent::industrial::flexible_ratio();
			bid_flex_industrial /= price_interval;

			industrial_profiles.HV[point_iter].bids.submitted_demand_flex(price_interval + 1) = bid_inflex_industrial;
			industrial_profiles.HV[point_iter].bids.submitted_demand_flex.segment(1, price_interval) = Eigen::VectorXd::Constant(price_interval, bid_flex_industrial);
			industrial_profiles.HV[point_iter].bids.redispatch_demand = industrial_profiles.HV[point_iter].bids.submitted_demand_flex;
		}

		return industrial_profiles;
	}
}

void agent::agents_set(power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
	Power_market_inform.agent_profiles.aggregators = aggregator_set(Power_market_inform.International_Market, Power_network_inform);
	Power_market_inform.agent_profiles.end_users = end_user_set(Power_market_inform, Power_network_inform);
	Power_market_inform.agent_profiles.industrial = industrial_set(Power_network_inform);
}
