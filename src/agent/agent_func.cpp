// Source file for operational strategies of different agents
#include "agent_func.h"

namespace{
	void agent_settlement_set(agent::settlement_struct &settlement){
		settlement.cost_demand.balancing = 0.;
		settlement.cost_demand.EOM = 0.;
		settlement.cost_demand.redispatch = 0.;
		settlement.cost_supply.balancing = 0.;
		settlement.cost_supply.EOM = 0.;
		settlement.cost_supply.redispatch = 0.;
		settlement.price.balancing = 0.;
		settlement.price.EOM = 0.;
		settlement.price.redispatch = 0.;
		settlement.reimburse.redispatch = 0.;
		settlement.reimburse.balancing = 0.;
		settlement.utility_demand.balancing = 0.;
		settlement.utility_demand.EOM = 0.;
		settlement.utility_demand.redispatch = 0.;
		settlement.utility_supply.balancing = 0.;
		settlement.utility_supply.EOM = 0.;
		settlement.utility_supply.redispatch = 0.;
//		settlement.volume_demand.balancing = 0.;
		settlement.volume_demand.EOM = 0.;
		settlement.volume_demand.BESS = 0.;
//		settlement.volume_supply.balancing = 0.;
		settlement.volume_supply.EOM = 0.;
		settlement.volume_supply.BESS = 0.;
		settlement.volume_demand_down.balancing = 0.;
		settlement.volume_demand_down.redispatch = 0.;
		settlement.volume_demand_down.imbalance = 0.;
		settlement.volume_demand_up.balancing = 0.;
		settlement.volume_demand_up.redispatch = 0.;
		settlement.volume_demand_up.imbalance = 0.;
		settlement.volume_supply_down.balancing = 0.;
		settlement.volume_supply_down.redispatch = 0.;
		settlement.volume_supply_down.imbalance = 0.;
		settlement.volume_supply_up.balancing = 0.;
		settlement.volume_supply_up.redispatch = 0.;
		settlement.volume_supply_up.imbalance = 0.;
	}

	void agent_results_set(agent::results_struct &results){
		results.cleared_supply = 0.;
		results.cleared_demand = 0.;
		results.confirmed_supply = 0.;
		results.confirmed_demand = 0.;
		results.imbalance_supply = 0.;
		results.imbalance_demand = 0.;
		results.actual_supply = 0.;
		results.actual_demand = 0.;
	}

	void agent_bids_initialization(agent::bids_struct &bids){
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
		bids.imbalance_supply = Eigen::VectorXd::Zero(price_interval + 2);
		bids.imbalance_demand = Eigen::VectorXd::Zero(price_interval + 2);
	}

	void agent_submitted_bids_scale(double scale, agent::bids_struct &bids){
		bids.submitted_supply_inflex *= scale;
		bids.submitted_demand_inflex *= scale;
		bids.submitted_supply_flex *= scale;
		bids.submitted_demand_flex *= scale;
	}

	void agent_scheduled_results_calculation(int bz_ID, int node_ID, int marginal_price_ID, int original_price_ID, power_market::market_whole_inform &Power_market_inform, agent::bids_struct &bids, agent::results_struct &results, bool cleared = 1){
		int price_interval = power_market::parameters::price_interval();

		if(marginal_price_ID > 0){
			results.confirmed_supply += bids.redispatch_supply.head(marginal_price_ID).sum();
		}
		if(marginal_price_ID < price_interval + 1){
			results.confirmed_demand += bids.redispatch_demand.tail(price_interval + 1 - marginal_price_ID).sum();
		}
		results.confirmed_supply += Power_market_inform.TSO_Market.confirmed.ratio_supply(node_ID) * bids.redispatch_supply(marginal_price_ID);
		results.confirmed_demand += Power_market_inform.TSO_Market.confirmed.ratio_demand(node_ID) * bids.redispatch_demand(marginal_price_ID);

		if(!cleared){
			return;
		}

		if(original_price_ID > 0){
			//results.cleared_supply += bids.filter_supply.head(original_price_ID).sum();
			results.cleared_supply += bids.submitted_supply_inflex.head(original_price_ID).sum();
			results.cleared_supply += bids.submitted_supply_flex.head(original_price_ID).sum();
		}
		if(original_price_ID < price_interval + 1){
			//results.cleared_demand += bids.filter_demand.tail(price_interval + 1 - original_price_ID).sum();
			results.cleared_demand += bids.submitted_demand_inflex.tail(price_interval + 1 - original_price_ID).sum();
			results.cleared_demand += bids.submitted_demand_flex.tail(price_interval + 1 - original_price_ID).sum();
		}
		results.cleared_supply += Power_market_inform.International_Market.confirmed.ratio_supply(bz_ID)  * bids.submitted_supply_inflex(original_price_ID);
		results.cleared_supply += Power_market_inform.International_Market.confirmed.ratio_supply(bz_ID)  * bids.submitted_supply_flex(original_price_ID);
		results.cleared_demand += Power_market_inform.International_Market.confirmed.ratio_demand(bz_ID) * bids.submitted_demand_inflex(original_price_ID);
		results.cleared_demand += Power_market_inform.International_Market.confirmed.ratio_demand(bz_ID) * bids.submitted_demand_flex(original_price_ID);
	}

	void agent_actual_results_calculation(int node_ID, int marginal_price_ID, power_market::market_whole_inform &Power_market_inform, agent::bids_struct &bids, agent::results_struct &results, bool control_reserve_flag){
		int price_interval = power_market::parameters::price_interval();

		if(control_reserve_flag){
			if(marginal_price_ID > 0){
				results.actual_supply += bids.balancing_supply.head(marginal_price_ID).sum();
			}
			if(marginal_price_ID < price_interval + 1){
				results.actual_demand += bids.balancing_demand.tail(price_interval + 1 - marginal_price_ID).sum();
			}
			results.actual_supply += Power_market_inform.TSO_Market.actual.ratio_supply(node_ID) * bids.balancing_supply(marginal_price_ID);
			results.actual_demand += Power_market_inform.TSO_Market.actual.ratio_demand(node_ID) * bids.balancing_demand(marginal_price_ID);
		}
		else{
			results.actual_supply = results.confirmed_supply;
			results.actual_demand= results.confirmed_demand;
		}
	}

	void agent_EOM_settlement_calculation(int tick, int node_ID, double original_price_supply, double original_price_demand, power_market::market_whole_inform &Power_market_inform, agent::bids_struct &bids, agent::results_struct &results, agent::settlement_struct &settlement, bool inflex_price = 0){
		int price_interval = power_market::parameters::price_interval();

		// Settlement of EOM
		// Supply side
		double cleared_supply_gap = results.cleared_supply;
		double margin_quan_supply;
		// Determine the markup cost if no supply quantity is served
		for(int price_iter = 0; price_iter < price_interval; ++ price_iter){
			double current_price = (1 - inflex_price) * Power_market_inform.price_map.bidded_price(price_iter);
			if(current_price >= 0.){
				break;
			}

			margin_quan_supply = bids.submitted_supply_inflex(price_iter);
			margin_quan_supply += bids.submitted_supply_flex(price_iter);
			settlement.cost_supply.EOM += -current_price * margin_quan_supply;
			Power_market_inform.TSO_Market.EOM.cost(tick, node_ID) += -current_price * margin_quan_supply;
		}
		for(int price_iter = 0; price_iter < price_interval; ++ price_iter){
			double current_price = (1 - inflex_price) * Power_market_inform.price_map.bidded_price(price_iter);

			margin_quan_supply = bids.submitted_supply_inflex(price_iter);
			margin_quan_supply += bids.submitted_supply_flex(price_iter);
			if(cleared_supply_gap > margin_quan_supply){
				cleared_supply_gap -= margin_quan_supply;
				settlement.volume_supply.EOM += margin_quan_supply;
				settlement.utility_supply.EOM += original_price_supply * margin_quan_supply;
				settlement.cost_supply.EOM += current_price * margin_quan_supply;
				Power_market_inform.TSO_Market.EOM.cost(tick, node_ID) += current_price * margin_quan_supply;
			}
			else{
				settlement.volume_supply.EOM += cleared_supply_gap;
				settlement.utility_supply.EOM += original_price_supply * cleared_supply_gap;
				settlement.cost_supply.EOM += current_price * cleared_supply_gap;
				Power_market_inform.TSO_Market.EOM.cost(tick, node_ID) += current_price * cleared_supply_gap;
				break;
			}
		}

		// Demand side
		double cleared_demand_gap = results.cleared_demand;
		double margin_quan_demand;
		int margin_ID_demand;
		for(int price_iter = price_interval + 1; price_iter >= 0; -- price_iter){
			double current_price;
			current_price = (1 - inflex_price) * Power_market_inform.price_map.bidded_price(price_iter);
			current_price += inflex_price * Power_market_inform.price_map.bidded_price(price_interval + 1);

			margin_quan_demand = bids.submitted_demand_inflex(price_iter);
			margin_quan_demand += bids.submitted_demand_flex(price_iter);

			if(cleared_demand_gap > margin_quan_demand){
				cleared_demand_gap -= margin_quan_demand;
				settlement.volume_demand.EOM += margin_quan_demand;
				settlement.utility_demand.EOM += current_price * margin_quan_demand;
				settlement.price.EOM += original_price_demand * margin_quan_demand;
				Power_market_inform.TSO_Market.EOM.utility(tick, node_ID) += current_price * margin_quan_demand;
			}
			else{
				settlement.volume_demand.EOM += cleared_demand_gap;
				settlement.utility_demand.EOM += current_price * cleared_demand_gap;
				settlement.price.EOM += original_price_demand * cleared_demand_gap;
				Power_market_inform.TSO_Market.EOM.utility(tick, node_ID) += current_price * cleared_demand_gap;
				break;
			}
		}
	}

	void end_user_redispatch_settlement_calculation(int tick, int node_ID, double original_price, power_market::market_whole_inform &Power_market_inform, agent::end_user::operation_struct &end_user, agent::aggregator::profile &aggregator, bool inflex_price = 0){
		int price_interval = power_market::parameters::price_interval();
		double redispatch_price_max = power_market::parameters::redispatch_price_max();

		// Settlement of redispatch
		// Supply side
		double cleared_supply_gap = end_user.results.cleared_supply;
		double confirmed_supply_gap = end_user.results.confirmed_supply;
		double min_supply_gap = std::min(cleared_supply_gap, confirmed_supply_gap);
		double max_supply_gap = std::max(cleared_supply_gap, confirmed_supply_gap);
		bool reduced_flag_supply = (confirmed_supply_gap == min_supply_gap);
		double margin_quan_supply;
		int margin_ID_supply;
		if(end_user.bids.submitted_supply_inflex.sum() < min_supply_gap){
			min_supply_gap -= end_user.bids.submitted_supply_inflex.sum();
			max_supply_gap -= end_user.bids.submitted_supply_inflex.sum();

			for(int price_iter = 0; price_iter < price_interval; ++ price_iter){
				margin_quan_supply = end_user.bids.submitted_supply_flex(price_iter);

				if(min_supply_gap > margin_quan_supply){
					min_supply_gap -= margin_quan_supply;
					max_supply_gap -= margin_quan_supply;
				}
				else{
					max_supply_gap -= min_supply_gap;
					margin_quan_supply -= min_supply_gap;
					margin_ID_supply = price_iter;
					break;
				}
			}
			for(int price_iter = margin_ID_supply; price_iter < price_interval; ++ price_iter){
				double current_price = Power_market_inform.price_map.bidded_price(price_iter);
				double real_price = (1 - inflex_price) * current_price;
				double redispatch_price = abs(original_price - current_price);
				if(reduced_flag_supply){
					redispatch_price = std::min(redispatch_price, redispatch_price_max);
				}

				if(price_iter > margin_ID_supply){
					margin_quan_supply = end_user.bids.submitted_supply_flex(price_iter);
				}

				if(max_supply_gap > margin_quan_supply){
					max_supply_gap -= margin_quan_supply;
					Power_market_inform.TSO_Market.redispatch.supply_down(tick, node_ID) += reduced_flag_supply * margin_quan_supply;
					Power_market_inform.TSO_Market.redispatch.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * margin_quan_supply;
					Power_market_inform.TSO_Market.redispatch.price_supply(tick, node_ID) += redispatch_price * margin_quan_supply;
					Power_market_inform.TSO_Market.redispatch.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * real_price * margin_quan_supply;
					end_user.settlement.volume_supply_down.redispatch +=  reduced_flag_supply * margin_quan_supply;
					end_user.settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * margin_quan_supply;
					end_user.settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * aggregator.price_supply_profile(0) * margin_quan_supply;
					end_user.settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * real_price * margin_quan_supply;
					aggregator.settlement.volume_supply_down.redispatch +=  reduced_flag_supply * margin_quan_supply;
					aggregator.settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * margin_quan_supply;
					aggregator.settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * original_price * margin_quan_supply;
					aggregator.settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * aggregator.price_supply_profile(0) * margin_quan_supply;
					aggregator.settlement.reimburse.redispatch += redispatch_price * margin_quan_supply;
				}
				else{
					Power_market_inform.TSO_Market.redispatch.supply_down(tick, node_ID) += reduced_flag_supply * max_supply_gap;
					Power_market_inform.TSO_Market.redispatch.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * max_supply_gap;
					Power_market_inform.TSO_Market.redispatch.price_supply(tick, node_ID) += redispatch_price * max_supply_gap;
					Power_market_inform.TSO_Market.redispatch.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * real_price * max_supply_gap;
					end_user.settlement.volume_supply_down.redispatch +=  reduced_flag_supply * max_supply_gap;
					end_user.settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * max_supply_gap;
					end_user.settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * aggregator.price_supply_profile(0) * max_supply_gap;
					end_user.settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * real_price * max_supply_gap;
					aggregator.settlement.volume_supply_down.redispatch +=  reduced_flag_supply * max_supply_gap;
					aggregator.settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * max_supply_gap;
					aggregator.settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * original_price * max_supply_gap;
					aggregator.settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * aggregator.price_supply_profile(0) * max_supply_gap;
					aggregator.settlement.reimburse.redispatch += redispatch_price * max_supply_gap;
					max_supply_gap = 0.;
					break;
				}
			}
		}
		else{
			for(int price_iter = 0; price_iter < price_interval; ++ price_iter){
				margin_quan_supply = end_user.bids.submitted_supply_inflex(price_iter);

				if(min_supply_gap > margin_quan_supply){
					min_supply_gap -= margin_quan_supply;
					max_supply_gap -= margin_quan_supply;
				}
				else{
					max_supply_gap -= min_supply_gap;
					margin_quan_supply -= min_supply_gap;
					margin_ID_supply = price_iter;
					break;
				}
			}
			for(int price_iter = margin_ID_supply; price_iter < price_interval; ++ price_iter){
				double current_price = Power_market_inform.price_map.bidded_price(price_iter);
				double real_price = (1 - inflex_price) * current_price;
				double redispatch_price = abs(original_price - current_price);
				if(reduced_flag_supply){
					redispatch_price = std::min(redispatch_price, redispatch_price_max);
				}

				if(price_iter > margin_ID_supply){
					margin_quan_supply = end_user.bids.submitted_supply_inflex(price_iter);
				}

				if(max_supply_gap > margin_quan_supply){
					max_supply_gap -= margin_quan_supply;
					Power_market_inform.TSO_Market.redispatch.supply_down(tick, node_ID) += reduced_flag_supply * margin_quan_supply;
					Power_market_inform.TSO_Market.redispatch.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * margin_quan_supply;
					Power_market_inform.TSO_Market.redispatch.price_supply(tick, node_ID) += redispatch_price * margin_quan_supply;
					Power_market_inform.TSO_Market.redispatch.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * real_price * margin_quan_supply;
					end_user.settlement.volume_supply_down.redispatch +=  reduced_flag_supply * margin_quan_supply;
					end_user.settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * margin_quan_supply;
					end_user.settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * aggregator.price_supply_profile(0) * margin_quan_supply;
					end_user.settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * real_price * margin_quan_supply;
					aggregator.settlement.volume_supply_down.redispatch +=  reduced_flag_supply * margin_quan_supply;
					aggregator.settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * margin_quan_supply;
					aggregator.settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * original_price * margin_quan_supply;
					aggregator.settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * aggregator.price_supply_profile(0) * margin_quan_supply;
					aggregator.settlement.reimburse.redispatch += redispatch_price * margin_quan_supply;
				}
				else{
					Power_market_inform.TSO_Market.redispatch.supply_down(tick, node_ID) += reduced_flag_supply * max_supply_gap;
					Power_market_inform.TSO_Market.redispatch.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * max_supply_gap;
					Power_market_inform.TSO_Market.redispatch.price_supply(tick, node_ID) += redispatch_price * max_supply_gap;
					Power_market_inform.TSO_Market.redispatch.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * real_price * max_supply_gap;
					end_user.settlement.volume_supply_down.redispatch +=  reduced_flag_supply * max_supply_gap;
					end_user.settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * max_supply_gap;
					end_user.settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * aggregator.price_supply_profile(0) * max_supply_gap;
					end_user.settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * real_price * max_supply_gap;
					aggregator.settlement.volume_supply_down.redispatch +=  reduced_flag_supply * max_supply_gap;
					aggregator.settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * max_supply_gap;
					aggregator.settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * original_price * max_supply_gap;
					aggregator.settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * aggregator.price_supply_profile(0) * max_supply_gap;
					aggregator.settlement.reimburse.redispatch += redispatch_price * max_supply_gap;
					max_supply_gap = 0.;
					break;
				}
			}

			if(max_supply_gap > 0.){
				for(int price_iter = 0; price_iter < price_interval; ++ price_iter){
					double current_price = Power_market_inform.price_map.bidded_price(price_iter);
					double real_price = (1 - inflex_price) * current_price;
					double redispatch_price = abs(original_price - current_price);
					if(reduced_flag_supply){
						redispatch_price = std::min(redispatch_price, redispatch_price_max);
					}

					margin_quan_supply = end_user.bids.submitted_supply_flex(price_iter);

					if(max_supply_gap > margin_quan_supply){
						max_supply_gap -= margin_quan_supply;
						Power_market_inform.TSO_Market.redispatch.supply_down(tick, node_ID) += reduced_flag_supply * margin_quan_supply;
						Power_market_inform.TSO_Market.redispatch.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * margin_quan_supply;
						Power_market_inform.TSO_Market.redispatch.price_supply(tick, node_ID) += redispatch_price * margin_quan_supply;
						Power_market_inform.TSO_Market.redispatch.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * real_price * margin_quan_supply;
						end_user.settlement.volume_supply_down.redispatch +=  reduced_flag_supply * margin_quan_supply;
						end_user.settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * margin_quan_supply;
						end_user.settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * aggregator.price_supply_profile(0) * margin_quan_supply;
						end_user.settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * real_price * margin_quan_supply;
						aggregator.settlement.volume_supply_down.redispatch +=  reduced_flag_supply * margin_quan_supply;
						aggregator.settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * margin_quan_supply;
						aggregator.settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * original_price * margin_quan_supply;
						aggregator.settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * aggregator.price_supply_profile(0) * margin_quan_supply;
						aggregator.settlement.reimburse.redispatch += redispatch_price * margin_quan_supply;
					}
					else{
						Power_market_inform.TSO_Market.redispatch.supply_down(tick, node_ID) += reduced_flag_supply * max_supply_gap;
						Power_market_inform.TSO_Market.redispatch.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * max_supply_gap;
						Power_market_inform.TSO_Market.redispatch.price_supply(tick, node_ID) += redispatch_price * max_supply_gap;
						Power_market_inform.TSO_Market.redispatch.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * real_price * max_supply_gap;
						end_user.settlement.volume_supply_down.redispatch +=  reduced_flag_supply * max_supply_gap;
						end_user.settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * max_supply_gap;
						end_user.settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * aggregator.price_supply_profile(0) * max_supply_gap;
						end_user.settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * real_price * max_supply_gap;
						aggregator.settlement.volume_supply_down.redispatch +=  reduced_flag_supply * max_supply_gap;
						aggregator.settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * max_supply_gap;
						aggregator.settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * original_price * max_supply_gap;
						aggregator.settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * aggregator.price_supply_profile(0) * max_supply_gap;
						aggregator.settlement.reimburse.redispatch += redispatch_price * max_supply_gap;
						max_supply_gap = 0.;
						break;
					}
				}
			}
		}

		// Demand side
		double cleared_demand_gap = end_user.results.cleared_demand;
		double confirmed_demand_gap = end_user.results.confirmed_demand;
		double min_demand_gap = std::min(cleared_demand_gap, confirmed_demand_gap);
		double max_demand_gap = std::max(cleared_demand_gap, confirmed_demand_gap);
		bool reduced_flag_demand = (confirmed_demand_gap == min_demand_gap);
		double margin_quan_demand;
		int margin_ID_demand;
		if(end_user.bids.submitted_demand_inflex.sum() < min_demand_gap){
			min_demand_gap -= end_user.bids.submitted_demand_inflex.sum();
			max_demand_gap -= end_user.bids.submitted_demand_inflex.sum();

			for(int price_iter = price_interval + 1; price_iter >= 0; -- price_iter){
				margin_quan_demand = end_user.bids.submitted_demand_flex(price_iter);

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
				double current_price = Power_market_inform.price_map.bidded_price(price_iter);
				double real_price = (1 - inflex_price) * current_price;
				real_price += inflex_price * Power_market_inform.price_map.bidded_price(price_interval + 1);
				double redispatch_price = abs(original_price - current_price);
				if(reduced_flag_demand){
					redispatch_price = std::min(redispatch_price, redispatch_price_max);
				}

				if(price_iter < margin_ID_demand){
					margin_quan_demand = end_user.bids.submitted_demand_flex(price_iter);
				}

				if(max_demand_gap > margin_quan_demand){
					max_demand_gap -= margin_quan_demand;
					Power_market_inform.TSO_Market.redispatch.demand_down(tick, node_ID) += reduced_flag_demand * margin_quan_demand;
					Power_market_inform.TSO_Market.redispatch.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * margin_quan_demand;
					Power_market_inform.TSO_Market.redispatch.price_demand(tick, node_ID) += redispatch_price * margin_quan_demand;
					Power_market_inform.TSO_Market.redispatch.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * real_price * margin_quan_demand;
					end_user.settlement.volume_demand_down.redispatch +=  reduced_flag_demand * margin_quan_demand;
					end_user.settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * margin_quan_demand;
					end_user.settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * real_price * margin_quan_demand;
					end_user.settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * aggregator.price_demand_profile(0) * margin_quan_demand;
					aggregator.settlement.volume_demand_down.redispatch +=  reduced_flag_demand * margin_quan_demand;
					aggregator.settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * margin_quan_demand;
					aggregator.settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * aggregator.price_demand_profile(0) * margin_quan_demand;
					aggregator.settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * original_price * margin_quan_demand;
					aggregator.settlement.reimburse.redispatch += redispatch_price * margin_quan_demand;
				}
				else{
					Power_market_inform.TSO_Market.redispatch.demand_down(tick, node_ID) += reduced_flag_demand * max_demand_gap;
					Power_market_inform.TSO_Market.redispatch.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * max_demand_gap;
					Power_market_inform.TSO_Market.redispatch.price_demand(tick, node_ID) += redispatch_price * max_demand_gap;
					Power_market_inform.TSO_Market.redispatch.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * real_price * max_demand_gap;
					end_user.settlement.volume_demand_down.redispatch +=  reduced_flag_demand * max_demand_gap;
					end_user.settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * max_demand_gap;
					end_user.settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * real_price * max_demand_gap;
					end_user.settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * aggregator.price_demand_profile(0) * max_demand_gap;
					aggregator.settlement.volume_demand_down.redispatch +=  reduced_flag_demand * max_demand_gap;
					aggregator.settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * max_demand_gap;
					aggregator.settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * aggregator.price_demand_profile(0) * max_demand_gap;
					aggregator.settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * original_price * max_demand_gap;
					aggregator.settlement.reimburse.redispatch += redispatch_price * max_demand_gap;
					max_demand_gap = 0.;
					break;
				}
			}
		}
		else{
			for(int price_iter = price_interval + 1; price_iter >= 0; -- price_iter){
				margin_quan_demand = end_user.bids.submitted_demand_inflex(price_iter);

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
				double current_price = Power_market_inform.price_map.bidded_price(price_iter);
				double real_price = (1 - inflex_price) * current_price;
				real_price += inflex_price * Power_market_inform.price_map.bidded_price(price_interval + 1);
				double redispatch_price = abs(original_price - current_price);
				if(reduced_flag_demand){
					redispatch_price = std::min(redispatch_price, redispatch_price_max);
				}

				if(price_iter < margin_ID_demand){
					margin_quan_demand = end_user.bids.submitted_demand_inflex(price_iter);
				}

				if(max_demand_gap > margin_quan_demand){
					max_demand_gap -= margin_quan_demand;
					Power_market_inform.TSO_Market.redispatch.demand_down(tick, node_ID) += reduced_flag_demand * margin_quan_demand;
					Power_market_inform.TSO_Market.redispatch.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * margin_quan_demand;
					Power_market_inform.TSO_Market.redispatch.price_demand(tick, node_ID) += redispatch_price * margin_quan_demand;
					Power_market_inform.TSO_Market.redispatch.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * real_price * margin_quan_demand;
					end_user.settlement.volume_demand_down.redispatch +=  reduced_flag_demand * margin_quan_demand;
					end_user.settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * margin_quan_demand;
					end_user.settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * real_price * margin_quan_demand;
					end_user.settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * aggregator.price_demand_profile(0) * margin_quan_demand;
					aggregator.settlement.volume_demand_down.redispatch +=  reduced_flag_demand * margin_quan_demand;
					aggregator.settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * margin_quan_demand;
					aggregator.settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * aggregator.price_demand_profile(0) * margin_quan_demand;
					aggregator.settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * original_price * margin_quan_demand;
					aggregator.settlement.reimburse.redispatch += redispatch_price * margin_quan_demand;
				}
				else{
					Power_market_inform.TSO_Market.redispatch.demand_down(tick, node_ID) += reduced_flag_demand * max_demand_gap;
					Power_market_inform.TSO_Market.redispatch.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * max_demand_gap;
					Power_market_inform.TSO_Market.redispatch.price_demand(tick, node_ID) += redispatch_price * max_demand_gap;
					Power_market_inform.TSO_Market.redispatch.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * real_price * max_demand_gap;
					end_user.settlement.volume_demand_down.redispatch +=  reduced_flag_demand * max_demand_gap;
					end_user.settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * max_demand_gap;
					end_user.settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * real_price * max_demand_gap;
					end_user.settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * aggregator.price_demand_profile(0) * max_demand_gap;
					aggregator.settlement.volume_demand_down.redispatch +=  reduced_flag_demand * max_demand_gap;
					aggregator.settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * max_demand_gap;
					aggregator.settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * aggregator.price_demand_profile(0) * max_demand_gap;
					aggregator.settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * original_price * max_demand_gap;
					aggregator.settlement.reimburse.redispatch += redispatch_price * max_demand_gap;
					max_demand_gap = 0.;
					break;
				}
			}

			if(max_demand_gap > 0.){
				for(int price_iter = price_interval + 1; price_iter >= 0; -- price_iter){
					double current_price = Power_market_inform.price_map.bidded_price(price_iter);
					double real_price = (1 - inflex_price) * current_price;
					real_price += inflex_price * Power_market_inform.price_map.bidded_price(price_interval + 1);
					double redispatch_price = abs(original_price - current_price);
					if(reduced_flag_demand){
						redispatch_price = std::min(redispatch_price, redispatch_price_max);
					}

					margin_quan_demand = end_user.bids.submitted_demand_flex(price_iter);

					if(max_demand_gap > margin_quan_demand){
						max_demand_gap -= margin_quan_demand;
						Power_market_inform.TSO_Market.redispatch.demand_down(tick, node_ID) += reduced_flag_demand * margin_quan_demand;
						Power_market_inform.TSO_Market.redispatch.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * margin_quan_demand;
						Power_market_inform.TSO_Market.redispatch.price_demand(tick, node_ID) += redispatch_price * margin_quan_demand;
						Power_market_inform.TSO_Market.redispatch.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * real_price * margin_quan_demand;
						end_user.settlement.volume_demand_down.redispatch +=  reduced_flag_demand * margin_quan_demand;
						end_user.settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * margin_quan_demand;
						end_user.settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * real_price * margin_quan_demand;
						end_user.settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * aggregator.price_demand_profile(0) * margin_quan_demand;
						aggregator.settlement.volume_demand_down.redispatch +=  reduced_flag_demand * margin_quan_demand;
						aggregator.settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * margin_quan_demand;
						aggregator.settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * aggregator.price_demand_profile(0) * margin_quan_demand;
						aggregator.settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * original_price * margin_quan_demand;
						aggregator.settlement.reimburse.redispatch += redispatch_price * margin_quan_demand;
					}
					else{
						Power_market_inform.TSO_Market.redispatch.demand_down(tick, node_ID) += reduced_flag_demand * max_demand_gap;
						Power_market_inform.TSO_Market.redispatch.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * max_demand_gap;
						Power_market_inform.TSO_Market.redispatch.price_demand(tick, node_ID) += redispatch_price * max_demand_gap;
						Power_market_inform.TSO_Market.redispatch.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * real_price * max_demand_gap;
						end_user.settlement.volume_demand_down.redispatch +=  reduced_flag_demand * max_demand_gap;
						end_user.settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * max_demand_gap;
						end_user.settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * real_price * max_demand_gap;
						end_user.settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * aggregator.price_demand_profile(0) * max_demand_gap;
						aggregator.settlement.volume_demand_down.redispatch +=  reduced_flag_demand * max_demand_gap;
						aggregator.settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * max_demand_gap;
						aggregator.settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * aggregator.price_demand_profile(0) * max_demand_gap;
						aggregator.settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * original_price * max_demand_gap;
						aggregator.settlement.reimburse.redispatch += redispatch_price * max_demand_gap;
						max_demand_gap = 0.;
						break;
					}
				}
			}
		}

		// Correct overestimated utility
//		double min_demand = std::min(end_user.results.cleared_demand, end_user.results.confirmed_demand);
//		double max_demand = std::max(end_user.results.cleared_demand, end_user.results.confirmed_demand);
//		min_demand = std::max(min_demand, end_user.direct_demand);
//		double over_est_utility_redispatch = (max_demand > min_demand) * (max_demand - min_demand);
//		over_est_utility_redispatch *= 1 - 2 * reduced_flag_demand;
//		over_est_utility_redispatch *= Power_market_inform.price_map.bidded_price(price_interval + 1);
//		end_user.settlement.utility_demand.redispatch -= over_est_utility_redispatch;
//		Power_market_inform.TSO_Market.redispatch.utility(tick, node_ID) -= over_est_utility_redispatch;
	}

	void agent_redispatch_settlement_calculation(int tick, int node_ID, double original_price, power_market::market_whole_inform &Power_market_inform, agent::bids_struct &bids, agent::results_struct &results, agent::settlement_struct&settlement, bool inflex_price = 0){
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
		if(bids.submitted_supply_inflex.sum() < min_supply_gap){
			min_supply_gap -= bids.submitted_supply_inflex.sum();
			max_supply_gap -= bids.submitted_supply_inflex.sum();

			for(int price_iter = 0; price_iter < price_interval; ++ price_iter){
				margin_quan_supply = bids.submitted_supply_flex(price_iter);

				if(min_supply_gap > margin_quan_supply){
					min_supply_gap -= margin_quan_supply;
					max_supply_gap -= margin_quan_supply;
				}
				else{
					max_supply_gap -= min_supply_gap;
					margin_quan_supply -= min_supply_gap;
					margin_ID_supply = price_iter;
					break;
				}
			}
			for(int price_iter = margin_ID_supply; price_iter < price_interval; ++ price_iter){
				double current_price = Power_market_inform.price_map.bidded_price(price_iter);
				double real_price = (1 - inflex_price) * current_price;
				double redispatch_price = abs(original_price - current_price);
				if(reduced_flag_supply){
					redispatch_price = std::min(redispatch_price, redispatch_price_max);
				}

				if(price_iter > margin_ID_supply){
					margin_quan_supply = bids.submitted_supply_flex(price_iter);
				}

				if(max_supply_gap > margin_quan_supply){
					max_supply_gap -= margin_quan_supply;
					Power_market_inform.TSO_Market.redispatch.supply_down(tick, node_ID) += reduced_flag_supply * margin_quan_supply;
					Power_market_inform.TSO_Market.redispatch.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * margin_quan_supply;
					Power_market_inform.TSO_Market.redispatch.price_supply(tick, node_ID) += redispatch_price * margin_quan_supply;
					Power_market_inform.TSO_Market.redispatch.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * real_price * margin_quan_supply;
					settlement.volume_supply_down.redispatch +=  reduced_flag_supply * margin_quan_supply;
					settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * margin_quan_supply;
					settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * original_price * margin_quan_supply;
					settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * real_price * margin_quan_supply;
					//settlement.reimburse.redispatch += (1 - 2 * reduced_flag_supply) * current_price * margin_quan_supply;
					settlement.reimburse.redispatch += redispatch_price * margin_quan_supply;
				}
				else{
					Power_market_inform.TSO_Market.redispatch.supply_down(tick, node_ID) += reduced_flag_supply * max_supply_gap;
					Power_market_inform.TSO_Market.redispatch.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * max_supply_gap;
					Power_market_inform.TSO_Market.redispatch.price_supply(tick, node_ID) += redispatch_price * max_supply_gap;
					Power_market_inform.TSO_Market.redispatch.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * real_price * max_supply_gap;
					settlement.volume_supply_down.redispatch +=  reduced_flag_supply * max_supply_gap;
					settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * max_supply_gap;
					settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * original_price * max_supply_gap;
					settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * real_price * max_supply_gap;
					//settlement.reimburse.redispatch += (1 - 2 * reduced_flag_supply) * current_price * max_supply_gap;
					settlement.reimburse.redispatch += redispatch_price * max_supply_gap;
					max_supply_gap = 0.;
					break;
				}
			}
		}
		else{
			for(int price_iter = 0; price_iter < price_interval; ++ price_iter){
				margin_quan_supply = bids.submitted_supply_inflex(price_iter);

				if(min_supply_gap > margin_quan_supply){
					min_supply_gap -= margin_quan_supply;
					max_supply_gap -= margin_quan_supply;
				}
				else{
					max_supply_gap -= min_supply_gap;
					margin_quan_supply -= min_supply_gap;
					margin_ID_supply = price_iter;
					break;
				}
			}
			for(int price_iter = margin_ID_supply; price_iter < price_interval; ++ price_iter){
				double current_price = Power_market_inform.price_map.bidded_price(price_iter);
				double real_price = (1 - inflex_price) * current_price;
				double redispatch_price = abs(original_price - current_price);
				if(reduced_flag_supply){
					redispatch_price = std::min(redispatch_price, redispatch_price_max);
				}

				if(price_iter > margin_ID_supply){
					margin_quan_supply = bids.submitted_supply_inflex(price_iter);
				}

				if(max_supply_gap > margin_quan_supply){
					max_supply_gap -= margin_quan_supply;
					Power_market_inform.TSO_Market.redispatch.supply_down(tick, node_ID) += reduced_flag_supply * margin_quan_supply;
					Power_market_inform.TSO_Market.redispatch.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * margin_quan_supply;
					Power_market_inform.TSO_Market.redispatch.price_supply(tick, node_ID) += redispatch_price * margin_quan_supply;
					Power_market_inform.TSO_Market.redispatch.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * real_price * margin_quan_supply;
					settlement.volume_supply_down.redispatch +=  reduced_flag_supply * margin_quan_supply;
					settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * margin_quan_supply;
					settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * original_price * margin_quan_supply;
					settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * real_price * margin_quan_supply;
					//settlement.reimburse.redispatch += (1 - 2 * reduced_flag_supply) * current_price * margin_quan_supply;
					settlement.reimburse.redispatch += redispatch_price * margin_quan_supply;
				}
				else{
					Power_market_inform.TSO_Market.redispatch.supply_down(tick, node_ID) += reduced_flag_supply * max_supply_gap;
					Power_market_inform.TSO_Market.redispatch.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * max_supply_gap;
					Power_market_inform.TSO_Market.redispatch.price_supply(tick, node_ID) += redispatch_price * max_supply_gap;
					Power_market_inform.TSO_Market.redispatch.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * real_price * max_supply_gap;
					settlement.volume_supply_down.redispatch +=  reduced_flag_supply * max_supply_gap;
					settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * max_supply_gap;
					settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * original_price * max_supply_gap;
					settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * real_price * max_supply_gap;
					settlement.price.redispatch += (1 - 2 * reduced_flag_supply) * current_price * max_supply_gap;
					//settlement.reimburse.redispatch += (1 - 2 * reduced_flag_supply) * current_price * max_supply_gap;
					settlement.reimburse.redispatch += redispatch_price * max_supply_gap;
					max_supply_gap = 0.;
					break;
				}
			}

			if(max_supply_gap > 0.){
				for(int price_iter = 0; price_iter < price_interval; ++ price_iter){
					double current_price = Power_market_inform.price_map.bidded_price(price_iter);
					double real_price = (1 - inflex_price) * current_price;
					double redispatch_price = abs(original_price - current_price);
					if(reduced_flag_supply){
						redispatch_price = std::min(redispatch_price, redispatch_price_max);
					}

					margin_quan_supply = bids.submitted_supply_flex(price_iter);

					if(max_supply_gap > margin_quan_supply){
						max_supply_gap -= margin_quan_supply;
						Power_market_inform.TSO_Market.redispatch.supply_down(tick, node_ID) += reduced_flag_supply * margin_quan_supply;
						Power_market_inform.TSO_Market.redispatch.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * margin_quan_supply;
						Power_market_inform.TSO_Market.redispatch.price_supply(tick, node_ID) += redispatch_price * margin_quan_supply;
						Power_market_inform.TSO_Market.redispatch.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * real_price * margin_quan_supply;
						settlement.volume_supply_down.redispatch +=  reduced_flag_supply * margin_quan_supply;
						settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * margin_quan_supply;
						settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * original_price * margin_quan_supply;
						settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * real_price * margin_quan_supply;
						//settlement.reimburse.redispatch += (1 - 2 * reduced_flag_supply) * current_price * margin_quan_supply;
						settlement.reimburse.redispatch += redispatch_price * margin_quan_supply;
					}
					else{
						Power_market_inform.TSO_Market.redispatch.supply_down(tick, node_ID) += reduced_flag_supply * max_supply_gap;
						Power_market_inform.TSO_Market.redispatch.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * max_supply_gap;
						Power_market_inform.TSO_Market.redispatch.price_supply(tick, node_ID) += redispatch_price * max_supply_gap;
						Power_market_inform.TSO_Market.redispatch.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * real_price * max_supply_gap;
						settlement.volume_supply_down.redispatch +=  reduced_flag_supply * max_supply_gap;
						settlement.volume_supply_up.redispatch +=  (1 - reduced_flag_supply) * max_supply_gap;
						settlement.utility_supply.redispatch += (1 - 2 * reduced_flag_supply) * original_price * max_supply_gap;
						settlement.cost_supply.redispatch += (1 - 2 * reduced_flag_supply) * real_price * max_supply_gap;
						//settlement.reimburse.redispatch += (1 - 2 * reduced_flag_supply) * current_price * margin_quan_supply;
						settlement.reimburse.redispatch += redispatch_price * max_supply_gap;
						max_supply_gap = 0.;
						break;
					}
				}
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
		if(bids.submitted_demand_inflex.sum() < min_demand_gap){
			min_demand_gap -= bids.submitted_demand_inflex.sum();
			max_demand_gap -= bids.submitted_demand_inflex.sum();

			for(int price_iter = price_interval + 1; price_iter >= 0; -- price_iter){
				margin_quan_demand = bids.submitted_demand_flex(price_iter);

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
				double current_price = Power_market_inform.price_map.bidded_price(price_iter);
				double real_price = (1 - inflex_price) * current_price;
				real_price += inflex_price * Power_market_inform.price_map.bidded_price(price_interval + 1);
				double redispatch_price = abs(original_price - current_price);
				if(reduced_flag_demand){
					redispatch_price = std::min(redispatch_price, redispatch_price_max);
				}

				if(price_iter < margin_ID_demand){
					margin_quan_demand = bids.submitted_demand_flex(price_iter);
				}

				if(max_demand_gap > margin_quan_demand){
					max_demand_gap -= margin_quan_demand;
					Power_market_inform.TSO_Market.redispatch.demand_down(tick, node_ID) += reduced_flag_demand * margin_quan_demand;
					Power_market_inform.TSO_Market.redispatch.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * margin_quan_demand;
					Power_market_inform.TSO_Market.redispatch.price_demand(tick, node_ID) += redispatch_price * margin_quan_demand;
					Power_market_inform.TSO_Market.redispatch.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * real_price * margin_quan_demand;
					settlement.volume_demand_down.redispatch +=  reduced_flag_demand * margin_quan_demand;
					settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * margin_quan_demand;
					settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * real_price * margin_quan_demand;
					settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * original_price * margin_quan_demand;
					settlement.reimburse.redispatch += redispatch_price * margin_quan_demand;
				}
				else{
					Power_market_inform.TSO_Market.redispatch.demand_down(tick, node_ID) += reduced_flag_demand * max_demand_gap;
					Power_market_inform.TSO_Market.redispatch.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * max_demand_gap;
					Power_market_inform.TSO_Market.redispatch.price_demand(tick, node_ID) += redispatch_price * max_demand_gap;
					Power_market_inform.TSO_Market.redispatch.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * real_price * max_demand_gap;
					settlement.volume_demand_down.redispatch +=  reduced_flag_demand * max_demand_gap;
					settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * max_demand_gap;
					settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * real_price * max_demand_gap;
					settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * original_price * max_demand_gap;
					settlement.reimburse.redispatch += redispatch_price * max_demand_gap;
					max_demand_gap = 0.;
					break;
				}
			}
		}
		else{
			for(int price_iter = price_interval + 1; price_iter >= 0; -- price_iter){
				margin_quan_demand = bids.submitted_demand_inflex(price_iter);

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
				double current_price = Power_market_inform.price_map.bidded_price(price_iter);
				double real_price = (1 - inflex_price) * current_price;
				real_price += inflex_price * Power_market_inform.price_map.bidded_price(price_interval + 1);
				double redispatch_price = abs(original_price - current_price);
				if(reduced_flag_demand){
					redispatch_price = std::min(redispatch_price, redispatch_price_max);
				}

				if(price_iter < margin_ID_demand){
					margin_quan_demand = bids.submitted_demand_inflex(price_iter);
				}

				if(max_demand_gap > margin_quan_demand){
					max_demand_gap -= margin_quan_demand;
					Power_market_inform.TSO_Market.redispatch.demand_down(tick, node_ID) += reduced_flag_demand * margin_quan_demand;
					Power_market_inform.TSO_Market.redispatch.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * margin_quan_demand;
					Power_market_inform.TSO_Market.redispatch.price_demand(tick, node_ID) += redispatch_price * margin_quan_demand;
					Power_market_inform.TSO_Market.redispatch.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * real_price * margin_quan_demand;
					settlement.volume_demand_down.redispatch +=  reduced_flag_demand * margin_quan_demand;
					settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * margin_quan_demand;
					settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * real_price * margin_quan_demand;
					settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * original_price * margin_quan_demand;
					settlement.reimburse.redispatch += redispatch_price * margin_quan_demand;
				}
				else{
					Power_market_inform.TSO_Market.redispatch.demand_down(tick, node_ID) += reduced_flag_demand * max_demand_gap;
					Power_market_inform.TSO_Market.redispatch.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * max_demand_gap;
					Power_market_inform.TSO_Market.redispatch.price_demand(tick, node_ID) += redispatch_price * max_demand_gap;
					Power_market_inform.TSO_Market.redispatch.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * real_price * max_demand_gap;
					settlement.volume_demand_down.redispatch +=  reduced_flag_demand * max_demand_gap;
					settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * max_demand_gap;
					settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * real_price * max_demand_gap;
					settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * original_price * max_demand_gap;
					settlement.reimburse.redispatch += redispatch_price * max_demand_gap;
					max_demand_gap = 0.;
					break;
				}
			}

			if(max_demand_gap > 0.){
				for(int price_iter = price_interval + 1; price_iter >= 0; -- price_iter){
					double current_price = Power_market_inform.price_map.bidded_price(price_iter);
					double real_price = (1 - inflex_price) * current_price;
					real_price += Power_market_inform.price_map.bidded_price(price_interval + 1);
					double redispatch_price = abs(original_price - current_price);
					if(reduced_flag_demand){
						redispatch_price = std::min(redispatch_price, redispatch_price_max);
					}

					margin_quan_demand = bids.submitted_demand_flex(price_iter);

					if(max_demand_gap > margin_quan_demand){
						max_demand_gap -= margin_quan_demand;
						Power_market_inform.TSO_Market.redispatch.demand_down(tick, node_ID) += reduced_flag_demand * margin_quan_demand;
						Power_market_inform.TSO_Market.redispatch.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * margin_quan_demand;
						Power_market_inform.TSO_Market.redispatch.price_demand(tick, node_ID) += redispatch_price * margin_quan_demand;
						Power_market_inform.TSO_Market.redispatch.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * real_price * margin_quan_demand;
						settlement.volume_demand_down.redispatch +=  reduced_flag_demand * margin_quan_demand;
						settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * margin_quan_demand;
						settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * real_price * margin_quan_demand;
						settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * original_price * margin_quan_demand;
						settlement.reimburse.redispatch += redispatch_price * margin_quan_demand;
					}
					else{
						Power_market_inform.TSO_Market.redispatch.demand_down(tick, node_ID) += reduced_flag_demand * max_demand_gap;
						Power_market_inform.TSO_Market.redispatch.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * max_demand_gap;
						Power_market_inform.TSO_Market.redispatch.price_demand(tick, node_ID) += redispatch_price * max_demand_gap;
						Power_market_inform.TSO_Market.redispatch.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * real_price * max_demand_gap;
						settlement.volume_demand_down.redispatch +=  reduced_flag_demand * max_demand_gap;
						settlement.volume_demand_up.redispatch +=  (1 - reduced_flag_demand) * max_demand_gap;
						settlement.utility_demand.redispatch += (1 - 2 * reduced_flag_demand) * real_price * max_demand_gap;
						settlement.cost_demand.redispatch += (1 - 2 * reduced_flag_demand) * original_price * max_demand_gap;
						settlement.reimburse.redispatch += redispatch_price * max_demand_gap;
						max_demand_gap = 0.;
						break;
					}
				}
			}
		}
	}

	void agent_balancing_settlement_calculation(int tick, int node_ID, power_market::market_whole_inform &Power_market_inform, agent::bids_struct &bids, agent::results_struct &results, agent::settlement_struct &settlement, bool inflex_price = 0){
		int price_interval = power_market::parameters::price_interval();
		double balancing_price_max = power_market::parameters::balancing_price_max();

		// Settlement of balancing
		// Supply side
		double imbalance_supply_ratio = results.imbalance_supply / (results.confirmed_supply + 1E-16);
		double confirmed_supply_gap = results.confirmed_supply * (1. + imbalance_supply_ratio);
		double actual_supply_gap = results.actual_supply;
		double min_supply_gap = std::min(confirmed_supply_gap, actual_supply_gap);
		double max_supply_gap = std::max(confirmed_supply_gap, actual_supply_gap);
		bool reduced_flag_supply = (actual_supply_gap == min_supply_gap);
		double margin_quan_supply;
		int margin_ID_supply;
		if(bids.submitted_supply_inflex.sum() * (1. + imbalance_supply_ratio) < min_supply_gap){
			min_supply_gap -= bids.submitted_supply_inflex.sum() * (1. + imbalance_supply_ratio);
			max_supply_gap -= bids.submitted_supply_inflex.sum() * (1. + imbalance_supply_ratio);

			for(int price_iter = 0; price_iter < price_interval; ++ price_iter){
				margin_quan_supply = bids.submitted_supply_flex(price_iter) * (1. + imbalance_supply_ratio);

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
				double current_price = (1 - inflex_price) * Power_market_inform.price_map.bidded_price(price_iter);
				double balancing_price = Power_market_inform.price_map.bidded_price(price_iter);
				double balancing_price_up = std::min(balancing_price, balancing_price_max);
				double balancing_price_down = std::max(balancing_price, -balancing_price_max);

				if(price_iter > margin_ID_supply){
					margin_quan_supply = bids.submitted_supply_flex(price_iter) * (1. + imbalance_supply_ratio);
				}

				if(max_supply_gap > margin_quan_supply){
					max_supply_gap -= margin_quan_supply;
					Power_market_inform.TSO_Market.balancing.supply_down(tick, node_ID) += reduced_flag_supply * margin_quan_supply;
					Power_market_inform.TSO_Market.balancing.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * margin_quan_supply;
					Power_market_inform.TSO_Market.balancing.price_up(tick, node_ID) += (1 - reduced_flag_supply) * balancing_price_up * margin_quan_supply;
					Power_market_inform.TSO_Market.balancing.price_down(tick, node_ID) += -reduced_flag_supply * balancing_price_down * margin_quan_supply;
					Power_market_inform.TSO_Market.balancing.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * current_price * margin_quan_supply;
					settlement.volume_supply_down.balancing +=  reduced_flag_supply * margin_quan_supply;
					settlement.volume_supply_up.balancing +=  (1 - reduced_flag_supply) * margin_quan_supply;
					settlement.cost_supply.balancing += (1 - 2 * reduced_flag_supply) * current_price * margin_quan_supply;
					settlement.reimburse.balancing += (1 - reduced_flag_supply) * balancing_price_up * margin_quan_supply;
					settlement.reimburse.balancing += -reduced_flag_supply * balancing_price_down * margin_quan_supply;
				}
				else{
					Power_market_inform.TSO_Market.balancing.supply_down(tick, node_ID) += reduced_flag_supply * max_supply_gap;
					Power_market_inform.TSO_Market.balancing.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * max_supply_gap;
					Power_market_inform.TSO_Market.balancing.price_up(tick, node_ID) += (1 - reduced_flag_supply) * balancing_price_up * max_supply_gap;
					Power_market_inform.TSO_Market.balancing.price_down(tick, node_ID) += -reduced_flag_supply * balancing_price_down * max_supply_gap;
					Power_market_inform.TSO_Market.balancing.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * current_price * max_supply_gap;
					settlement.volume_supply_down.balancing +=  reduced_flag_supply * max_supply_gap;
					settlement.volume_supply_up.balancing +=  (1 - reduced_flag_supply) * max_supply_gap;
					settlement.cost_supply.balancing += (1 - 2 * reduced_flag_supply) * current_price * max_supply_gap;
					settlement.reimburse.balancing += (1 - reduced_flag_supply) * balancing_price_up * max_supply_gap;
					settlement.reimburse.balancing += -reduced_flag_supply * balancing_price_down * max_supply_gap;
					max_supply_gap = 0.;
					break;
				}
			}
		}
		else{
			for(int price_iter = 0; price_iter < price_interval; ++ price_iter){
				margin_quan_supply = bids.submitted_supply_inflex(price_iter) * (1. + imbalance_supply_ratio);

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
				double current_price = (1 - inflex_price) * Power_market_inform.price_map.bidded_price(price_iter);
				double balancing_price = Power_market_inform.price_map.bidded_price(price_iter);
				double balancing_price_up = std::min(balancing_price, balancing_price_max);
				double balancing_price_down = std::max(balancing_price, -balancing_price_max);

				if(price_iter > margin_ID_supply){
					margin_quan_supply = bids.submitted_supply_inflex(price_iter) * (1. + imbalance_supply_ratio);
				}

				if(max_supply_gap > margin_quan_supply){
					max_supply_gap -= margin_quan_supply;
					Power_market_inform.TSO_Market.balancing.supply_down(tick, node_ID) += reduced_flag_supply * margin_quan_supply;
					Power_market_inform.TSO_Market.balancing.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * margin_quan_supply;
					Power_market_inform.TSO_Market.balancing.price_up(tick, node_ID) += (1 - reduced_flag_supply) * balancing_price_up * margin_quan_supply;
					Power_market_inform.TSO_Market.balancing.price_down(tick, node_ID) += -reduced_flag_supply * balancing_price_down * margin_quan_supply;
					Power_market_inform.TSO_Market.balancing.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * current_price * margin_quan_supply;
					settlement.volume_supply_down.balancing +=  reduced_flag_supply * margin_quan_supply;
					settlement.volume_supply_up.balancing +=  (1 - reduced_flag_supply) * margin_quan_supply;
					settlement.cost_supply.balancing += (1 - 2 * reduced_flag_supply) * current_price * margin_quan_supply;
					settlement.reimburse.balancing += (1 - reduced_flag_supply) * balancing_price_up * margin_quan_supply;
					settlement.reimburse.balancing += -reduced_flag_supply * balancing_price_down * margin_quan_supply;
				}
				else{
					Power_market_inform.TSO_Market.balancing.supply_down(tick, node_ID) += reduced_flag_supply * max_supply_gap;
					Power_market_inform.TSO_Market.balancing.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * max_supply_gap;
					Power_market_inform.TSO_Market.balancing.price_up(tick, node_ID) += (1 - reduced_flag_supply) * balancing_price_up * max_supply_gap;
					Power_market_inform.TSO_Market.balancing.price_down(tick, node_ID) += -reduced_flag_supply * balancing_price_down * max_supply_gap;
					Power_market_inform.TSO_Market.balancing.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * current_price * max_supply_gap;
					settlement.volume_supply_down.balancing +=  reduced_flag_supply * max_supply_gap;
					settlement.volume_supply_up.balancing +=  (1 - reduced_flag_supply) * max_supply_gap;
					settlement.cost_supply.balancing += (1 - 2 * reduced_flag_supply) * current_price * max_supply_gap;
					settlement.reimburse.balancing += (1 - reduced_flag_supply) * balancing_price_up * max_supply_gap;
					settlement.reimburse.balancing += -reduced_flag_supply * balancing_price_down * max_supply_gap;
					max_supply_gap = 0.;
					break;
				}
			}

			if(max_supply_gap > 0.){
				for(int price_iter = 0; price_iter < price_interval; ++ price_iter){
					double current_price = (1 - inflex_price) * Power_market_inform.price_map.bidded_price(price_iter);
					double balancing_price = Power_market_inform.price_map.bidded_price(price_iter);
					double balancing_price_up = std::min(balancing_price, balancing_price_max);
					double balancing_price_down = std::max(balancing_price, -balancing_price_max);

					margin_quan_supply = bids.submitted_supply_flex(price_iter) * (1. + imbalance_supply_ratio);

					if(max_supply_gap > margin_quan_supply){
						max_supply_gap -= margin_quan_supply;
						Power_market_inform.TSO_Market.balancing.supply_down(tick, node_ID) += reduced_flag_supply * margin_quan_supply;
						Power_market_inform.TSO_Market.balancing.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * margin_quan_supply;
						Power_market_inform.TSO_Market.balancing.price_up(tick, node_ID) += (1 - reduced_flag_supply) * balancing_price_up * margin_quan_supply;
						Power_market_inform.TSO_Market.balancing.price_down(tick, node_ID) += -reduced_flag_supply * balancing_price_down * margin_quan_supply;
						Power_market_inform.TSO_Market.balancing.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * current_price * margin_quan_supply;
						settlement.volume_supply_down.balancing +=  reduced_flag_supply * margin_quan_supply;
						settlement.volume_supply_up.balancing +=  (1 - reduced_flag_supply) * margin_quan_supply;
						settlement.cost_supply.balancing += (1 - 2 * reduced_flag_supply) * current_price * margin_quan_supply;
						settlement.reimburse.balancing += (1 - reduced_flag_supply) * balancing_price_up * margin_quan_supply;
						settlement.reimburse.balancing += -reduced_flag_supply * balancing_price_down * margin_quan_supply;
					}
					else{
						Power_market_inform.TSO_Market.balancing.supply_down(tick, node_ID) += reduced_flag_supply * max_supply_gap;
						Power_market_inform.TSO_Market.balancing.supply_up(tick, node_ID) += (1 - reduced_flag_supply) * max_supply_gap;
						Power_market_inform.TSO_Market.balancing.price_up(tick, node_ID) += (1 - reduced_flag_supply) * balancing_price_up * max_supply_gap;
						Power_market_inform.TSO_Market.balancing.price_down(tick, node_ID) += -reduced_flag_supply * balancing_price_down * max_supply_gap;
						Power_market_inform.TSO_Market.balancing.cost(tick, node_ID) += (1 - 2 * reduced_flag_supply) * current_price * max_supply_gap;
						settlement.volume_supply_down.balancing +=  reduced_flag_supply * max_supply_gap;
						settlement.volume_supply_up.balancing +=  (1 - reduced_flag_supply) * max_supply_gap;
						settlement.cost_supply.balancing += (1 - 2 * reduced_flag_supply) * current_price * max_supply_gap;
						settlement.reimburse.balancing += (1 - reduced_flag_supply) * balancing_price_up * max_supply_gap;
						settlement.reimburse.balancing += -reduced_flag_supply * balancing_price_down * max_supply_gap;
						max_supply_gap = 0.;
						break;
					}
				}
			}
		}

		// Demand side
		double imbalance_demand_ratio = results.imbalance_demand / (results.confirmed_demand + 1E-16);
		double confirmed_demand_gap = results.confirmed_demand * (1. + imbalance_demand_ratio);
		double actual_demand_gap = results.actual_demand;
		double min_demand_gap = std::min(confirmed_demand_gap, actual_demand_gap);
		double max_demand_gap = std::max(confirmed_demand_gap, actual_demand_gap);
		bool reduced_flag_demand = (actual_demand_gap == min_demand_gap);
		double margin_quan_demand;
		int margin_ID_demand;

		if(bids.submitted_demand_inflex.sum() * (1. + imbalance_demand_ratio) < min_demand_gap){
			min_demand_gap -= bids.submitted_demand_inflex.sum() * (1. + imbalance_demand_ratio);
			max_demand_gap -= bids.submitted_demand_inflex.sum() * (1. + imbalance_demand_ratio);

			for(int price_iter = price_interval + 1; price_iter >= 0; -- price_iter){
				margin_quan_demand = bids.submitted_demand_flex(price_iter) * (1. + imbalance_demand_ratio);

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
				double current_price = (1 - inflex_price) * Power_market_inform.price_map.bidded_price(price_iter);
				current_price += inflex_price * Power_market_inform.price_map.bidded_price(price_interval + 1);
				double balancing_price = Power_market_inform.price_map.bidded_price(price_iter);
				double balancing_price_up = std::min(balancing_price, balancing_price_max);
				double balancing_price_down = std::max(balancing_price, -balancing_price_max);

				if(price_iter < margin_ID_demand){
					margin_quan_demand = bids.submitted_demand_flex(price_iter) * (1. + imbalance_demand_ratio);
				}

				if(max_demand_gap > margin_quan_demand){
					max_demand_gap -= margin_quan_demand;
					Power_market_inform.TSO_Market.balancing.demand_down(tick, node_ID) += reduced_flag_demand * margin_quan_demand;
					Power_market_inform.TSO_Market.balancing.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * margin_quan_demand;
					Power_market_inform.TSO_Market.balancing.price_up(tick, node_ID) += reduced_flag_demand * balancing_price_up * margin_quan_demand;
					Power_market_inform.TSO_Market.balancing.price_down(tick, node_ID) += -(1 - reduced_flag_demand) * balancing_price_down * margin_quan_demand;
					Power_market_inform.TSO_Market.balancing.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * current_price * margin_quan_demand;
					settlement.volume_demand_down.balancing +=  reduced_flag_demand * margin_quan_demand;
					settlement.volume_demand_up.balancing +=  (1 - reduced_flag_demand) * margin_quan_demand;
					settlement.utility_demand.balancing += (1 - 2 * reduced_flag_demand) * current_price * margin_quan_demand;
					settlement.reimburse.balancing += reduced_flag_demand * balancing_price_up * margin_quan_demand;
					settlement.reimburse.balancing += -(1 - reduced_flag_demand) * balancing_price_down * margin_quan_demand;
				}
				else{
					Power_market_inform.TSO_Market.balancing.demand_down(tick, node_ID) += reduced_flag_demand * max_demand_gap;
					Power_market_inform.TSO_Market.balancing.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * max_demand_gap;
					Power_market_inform.TSO_Market.balancing.price_up(tick, node_ID) += reduced_flag_demand * balancing_price_up * max_demand_gap;
					Power_market_inform.TSO_Market.balancing.price_down(tick, node_ID) += -(1 - reduced_flag_demand) * balancing_price_down * max_demand_gap;
					Power_market_inform.TSO_Market.balancing.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * current_price * max_demand_gap;
					settlement.volume_demand_down.balancing +=  reduced_flag_demand * max_demand_gap;
					settlement.volume_demand_up.balancing +=  (1 - reduced_flag_demand) * max_demand_gap;
					settlement.utility_demand.balancing += (1 - 2 * reduced_flag_demand) * current_price * max_demand_gap;
					settlement.reimburse.balancing += reduced_flag_demand * balancing_price_up * max_demand_gap;
					settlement.reimburse.balancing += -(1 - reduced_flag_demand) * balancing_price_down * max_demand_gap;
					max_demand_gap = 0.;
					break;
				}
			}
		}
		else{
			for(int price_iter = price_interval + 1; price_iter >= 0; -- price_iter){
				margin_quan_demand = bids.submitted_demand_inflex(price_iter) * (1. + imbalance_demand_ratio);

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
				double current_price = (1 - inflex_price) * Power_market_inform.price_map.bidded_price(price_iter);
				current_price += inflex_price * Power_market_inform.price_map.bidded_price(price_interval + 1);
				double balancing_price = Power_market_inform.price_map.bidded_price(price_iter);
				double balancing_price_up = std::min(balancing_price, balancing_price_max);
				double balancing_price_down = std::max(balancing_price, -balancing_price_max);

				if(price_iter < margin_ID_demand){
					margin_quan_demand = bids.submitted_demand_inflex(price_iter) * (1. + imbalance_demand_ratio);
				}

				if(max_demand_gap > margin_quan_demand){
					max_demand_gap -= margin_quan_demand;
					Power_market_inform.TSO_Market.balancing.demand_down(tick, node_ID) += reduced_flag_demand * margin_quan_demand;
					Power_market_inform.TSO_Market.balancing.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * margin_quan_demand;
					Power_market_inform.TSO_Market.balancing.price_up(tick, node_ID) += reduced_flag_demand * balancing_price_up * margin_quan_demand;
					Power_market_inform.TSO_Market.balancing.price_down(tick, node_ID) += -(1 - reduced_flag_demand) * balancing_price_down * margin_quan_demand;
					Power_market_inform.TSO_Market.balancing.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * current_price * margin_quan_demand;
					settlement.volume_demand_down.balancing +=  reduced_flag_demand * margin_quan_demand;
					settlement.volume_demand_up.balancing +=  (1 - reduced_flag_demand) * margin_quan_demand;
					settlement.utility_demand.balancing += (1 - 2 * reduced_flag_demand) * current_price * margin_quan_demand;
					settlement.reimburse.balancing += reduced_flag_demand * balancing_price_up * margin_quan_demand;
					settlement.reimburse.balancing += -(1 - reduced_flag_demand) * balancing_price_down * margin_quan_demand;
				}
				else{
					Power_market_inform.TSO_Market.balancing.demand_down(tick, node_ID) += reduced_flag_demand * max_demand_gap;
					Power_market_inform.TSO_Market.balancing.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * max_demand_gap;
					Power_market_inform.TSO_Market.balancing.price_up(tick, node_ID) += reduced_flag_demand * balancing_price_up * max_demand_gap;
					Power_market_inform.TSO_Market.balancing.price_down(tick, node_ID) += -(1 - reduced_flag_demand) * balancing_price_down * max_demand_gap;
					Power_market_inform.TSO_Market.balancing.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * current_price * max_demand_gap;
					settlement.volume_demand_down.balancing +=  reduced_flag_demand * max_demand_gap;
					settlement.volume_demand_up.balancing +=  (1 - reduced_flag_demand) * max_demand_gap;
					settlement.utility_demand.balancing += (1 - 2 * reduced_flag_demand) * current_price * max_demand_gap;
					settlement.reimburse.balancing += reduced_flag_demand * balancing_price_up * max_demand_gap;
					settlement.reimburse.balancing += -(1 - reduced_flag_demand) * balancing_price_down * max_demand_gap;
					max_demand_gap = 0.;
					break;
				}
			}

			if(max_demand_gap > 0.){
				for(int price_iter = price_interval + 1; price_iter >= 0; -- price_iter){
					double current_price = (1 - inflex_price) * Power_market_inform.price_map.bidded_price(price_iter);
					current_price += inflex_price * Power_market_inform.price_map.bidded_price(price_interval + 1);
					double balancing_price = Power_market_inform.price_map.bidded_price(price_iter);
					double balancing_price_up = std::min(balancing_price, balancing_price_max);
					double balancing_price_down = std::max(balancing_price, -balancing_price_max);

					margin_quan_demand = bids.submitted_demand_flex(price_iter) * (1. + imbalance_demand_ratio);

					if(max_demand_gap > margin_quan_demand){
						max_demand_gap -= margin_quan_demand;
						Power_market_inform.TSO_Market.balancing.demand_down(tick, node_ID) += reduced_flag_demand * margin_quan_demand;
						Power_market_inform.TSO_Market.balancing.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * margin_quan_demand;
						Power_market_inform.TSO_Market.balancing.price_up(tick, node_ID) += reduced_flag_demand * balancing_price_up * margin_quan_demand;
						Power_market_inform.TSO_Market.balancing.price_down(tick, node_ID) += -(1 - reduced_flag_demand) * balancing_price_down * margin_quan_demand;
						Power_market_inform.TSO_Market.balancing.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * current_price * margin_quan_demand;
						settlement.volume_demand_down.balancing +=  reduced_flag_demand * margin_quan_demand;
						settlement.volume_demand_up.balancing +=  (1 - reduced_flag_demand) * margin_quan_demand;
						settlement.utility_demand.balancing += (1 - 2 * reduced_flag_demand) * current_price * margin_quan_demand;
						settlement.reimburse.balancing += reduced_flag_demand * balancing_price_up * margin_quan_demand;
						settlement.reimburse.balancing += -(1 - reduced_flag_demand) * balancing_price_down * margin_quan_demand;
					}
					else{
						Power_market_inform.TSO_Market.balancing.demand_down(tick, node_ID) += reduced_flag_demand * max_demand_gap;
						Power_market_inform.TSO_Market.balancing.demand_up(tick, node_ID) += (1 - reduced_flag_demand) * max_demand_gap;
						Power_market_inform.TSO_Market.balancing.price_up(tick, node_ID) += reduced_flag_demand * balancing_price_up * max_demand_gap;
						Power_market_inform.TSO_Market.balancing.price_down(tick, node_ID) += -(1 - reduced_flag_demand) * balancing_price_down * max_demand_gap;
						Power_market_inform.TSO_Market.balancing.utility(tick, node_ID) += (1 - 2 * reduced_flag_demand) * current_price * max_demand_gap;
						settlement.volume_demand_down.balancing +=  reduced_flag_demand * max_demand_gap;
						settlement.volume_demand_up.balancing +=  (1 - reduced_flag_demand) * max_demand_gap;
						settlement.utility_demand.balancing += (1 - 2 * reduced_flag_demand) * current_price * max_demand_gap;
						settlement.reimburse.balancing += reduced_flag_demand * balancing_price_up * max_demand_gap;
						settlement.reimburse.balancing += -(1 - reduced_flag_demand) * balancing_price_down * max_demand_gap;
						max_demand_gap = 0.;
						break;
					}
				}
			}
		}
	}

	void market_operation_update(int tick, int bz_ID, power_market::schedule &agent, agent::results_struct &results){
		agent.EOM(tick, bz_ID) += results.cleared_supply - results.cleared_demand;

		agent.redispatch(tick, bz_ID) += results.confirmed_supply - results.cleared_supply;
		agent.redispatch(tick, bz_ID) -= results.confirmed_demand - results.cleared_demand;

		agent.balancing(tick, bz_ID) += results.actual_supply - results.confirmed_supply;
		agent.balancing(tick, bz_ID) -= results.actual_demand - results.confirmed_demand;
		agent.balancing(tick, bz_ID) += results.imbalance_demand - results.imbalance_supply;
	}

	agent::aggregator::profiles aggregator_set(int start_time, power_market::market_inform &International_Market, power_network::network_inform &Power_network_inform){
		int foresight_time = agent::aggregator::parameters::foresight_time();
		int point_num = Power_network_inform.points.bidding_zone.size();
		double arbitrage_demand = agent::aggregator::parameters::arbitrage_demand();
		double arbitrage_supply = agent::aggregator::parameters::arbitrage_supply();

		agent::aggregator::profiles aggregator_profiles(point_num);
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int bz_ID = Power_network_inform.points.bidding_zone(point_iter);

			// Set bids and results information
			agent_bids_initialization(aggregator_profiles[point_iter].bids);
			agent_results_set(aggregator_profiles[point_iter].results);
			agent_settlement_set(aggregator_profiles[point_iter].settlement);

			aggregator_profiles[point_iter].point_ID = point_iter;
			aggregator_profiles[point_iter].arbitrage_demand = arbitrage_demand;
			aggregator_profiles[point_iter].arbitrage_supply = arbitrage_supply;
			aggregator_profiles[point_iter].price_expected_profile = International_Market.confirmed.price.col(bz_ID).segment(start_time, foresight_time);
			aggregator_profiles[point_iter].price_demand_profile = International_Market.confirmed.price.col(bz_ID).segment(start_time, foresight_time).array() + aggregator_profiles[point_iter].arbitrage_demand;
			aggregator_profiles[point_iter].price_supply_profile = International_Market.confirmed.price.col(bz_ID).segment(start_time, foresight_time).array() - aggregator_profiles[point_iter].arbitrage_supply;
		}

		return aggregator_profiles;
	}

	agent::cross_border::edge_profiles cross_border_set(power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int edge_num = Power_market_inform.International_Market.network.num_edges;
		agent::cross_border::edge_profiles cross_border_profiles(edge_num);
		for(int edge_iter = 0; edge_iter < edge_num; ++ edge_iter){
			cross_border_profiles[edge_iter].node_num = Power_network_inform.cbt.entry_node_num[edge_iter];
			cross_border_profiles[edge_iter].entry_bz_ID = Power_network_inform.cbt.entry_bz(edge_iter);

			int node_num = cross_border_profiles[edge_iter].node_num;
			cross_border_profiles[edge_iter].profiles = agent::cross_border::profiles_typedef(node_num);
			for(int node_iter = 0; node_iter < node_num; ++ node_iter){
				// Set entry node
				cross_border_profiles[edge_iter].profiles[node_iter].node_ID = Power_network_inform.cbt.entry_nodes(edge_iter, node_iter);

				// Set settlement and results information
				agent_results_set(cross_border_profiles[edge_iter].profiles[node_iter].results);
				agent_settlement_set(cross_border_profiles[edge_iter].profiles[node_iter].settlement);
			}
		}

		return cross_border_profiles;
	}

	agent::end_user::profiles end_user_set(int start_time, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, configuration::process_config &process_par){
		int foresight_time = agent::end_user::parameters::foresight_time();
		int point_num = Power_network_inform.points.bidding_zone.size();
		int load_shift_time = agent::end_user::parameters::load_shift_time();
		int price_interval = power_market::parameters::price_interval();
		double residential_ratio = agent::parameters::residential_ratio();

		agent::end_user::profiles end_user_profiles(point_num);
		//int sample_num = agent::end_user::parameters::sample_num();
		int sample_num = Power_market_inform.agent_profiles.end_user_type.sample_num;
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			end_user_profiles[point_iter] = std::vector <agent::end_user::profile> (sample_num);
		}

		// Initialization of forecast demand profile and operation strategies
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
            int node_ID = Power_network_inform.points.node(point_iter);

			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				// Initialization of investment parameters
				//end_user_profiles[point_iter][sample_iter].investment.decision.dynamic_tariff = Power_market_inform.agent_profiles.end_user_type.dynamic_tariff[sample_iter];
				end_user_profiles[point_iter][sample_iter].investment.decision.smart_appliance = Power_market_inform.agent_profiles.end_user_type.smart_appliance[sample_iter];
				end_user_profiles[point_iter][sample_iter].investment.decision.PV = (Power_market_inform.agent_profiles.end_user_type.PV_scale[sample_iter] != 0.);
				end_user_profiles[point_iter][sample_iter].investment.decision.BESS = (Power_market_inform.agent_profiles.end_user_type.BESS_energy[sample_iter] != 0.);
				end_user_profiles[point_iter][sample_iter].investment.decision.EV_self_charging = (Power_market_inform.agent_profiles.end_user_type.EV_energy[sample_iter] != 0.);
				end_user_profiles[point_iter][sample_iter].investment.decision.reverse_flow = 1;
				end_user_profiles[point_iter][sample_iter].investment.decision.redispatch = Power_market_inform.agent_profiles.end_user_type.redispatch[sample_iter];
				end_user_profiles[point_iter][sample_iter].investment.decision.control_reserve = Power_market_inform.agent_profiles.end_user_type.control_reserve[sample_iter];
                end_user_profiles[point_iter][sample_iter].investment.decision.contingency = Power_market_inform.agent_profiles.end_user_type.contingency[sample_iter];

				// Initialization of operational parameters
				end_user_profiles[point_iter][sample_iter].operation.foresight_time = foresight_time;
				end_user_profiles[point_iter][sample_iter].operation.weight = Power_market_inform.agent_profiles.end_user_type.weight[sample_iter];
				end_user_profiles[point_iter][sample_iter].operation.PV_scale = Power_market_inform.agent_profiles.end_user_type.PV_scale[sample_iter];
				int load_shift_time_temp = std::min(load_shift_time, foresight_time / 2);
				end_user_profiles[point_iter][sample_iter].operation.smart_appliance.shift_time = load_shift_time_temp;
				end_user_profiles[point_iter][sample_iter].operation.BESS.energy_scale = Power_market_inform.agent_profiles.end_user_type.BESS_energy[sample_iter];
				end_user_profiles[point_iter][sample_iter].operation.BESS.capacity_scale = Power_market_inform.agent_profiles.end_user_type.BESS_capacity[sample_iter];
				//default .5 * E_max
				end_user_profiles[point_iter][sample_iter].operation.BESS.soc = end_user_profiles[point_iter][sample_iter].operation.BESS.energy_scale / 2;
				//end_user_profiles[point_iter][sample_iter].operation.BESS.soc = 0.;
				end_user_profiles[point_iter][sample_iter].operation.BESS.soc *= end_user_profiles[point_iter][sample_iter].investment.decision.BESS;
				end_user_profiles[point_iter][sample_iter].operation.EV.BESS.energy_scale = Power_market_inform.agent_profiles.end_user_type.EV_energy[sample_iter];
				end_user_profiles[point_iter][sample_iter].operation.EV.BESS.capacity_scale = Power_market_inform.agent_profiles.end_user_type.EV_capacity[sample_iter];
				//default .5 * E_max
				end_user_profiles[point_iter][sample_iter].operation.EV.BESS.soc = end_user_profiles[point_iter][sample_iter].operation.EV.BESS.energy_scale / 2;
				//end_user_profiles[point_iter][sample_iter].operation.EV.BESS.soc = 0.;
				end_user_profiles[point_iter][sample_iter].operation.EV.BESS.soc *= end_user_profiles[point_iter][sample_iter].investment.decision.EV_self_charging;

				// Initialization of input profiles
				end_user_profiles[point_iter][sample_iter].operation.EV.house_default_period = end_user_profiles[point_iter][sample_iter].operation.EV.house_schedule(start_time);
				end_user_profiles[point_iter][sample_iter].operation.EV.default_demand_profile = end_user_profiles[point_iter][sample_iter].operation.EV.demand_profile(start_time);
				end_user_profiles[point_iter][sample_iter].operation.EV.default_demand_profile *= end_user_profiles[point_iter][sample_iter].investment.decision.EV_self_charging;
				end_user_profiles[point_iter][sample_iter].operation.default_demand_profile = Power_network_inform.points.nominal_mean_demand_field.row(point_iter).segment(start_time, foresight_time);
				end_user_profiles[point_iter][sample_iter].operation.default_demand_profile *= agent::parameters::residential_ratio();
				end_user_profiles[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand = Eigen::VectorXd::Zero(foresight_time + load_shift_time_temp);
				for(int tick = load_shift_time_temp; tick < foresight_time + load_shift_time_temp; ++ tick){
					if(start_time + tick - load_shift_time_temp < 0){
						end_user_profiles[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand(tick) = 0.;
					}
					else{
						end_user_profiles[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand(tick) = Power_network_inform.points.nominal_mean_demand_field(point_iter, start_time + tick - load_shift_time_temp) * agent::parameters::residential_ratio();
						end_user_profiles[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand(tick) *= end_user_profiles[point_iter][sample_iter].investment.decision.smart_appliance * end_user_profiles[point_iter][sample_iter].operation.smart_appliance.scale;
					}
				}
				if(end_user_profiles[point_iter][sample_iter].operation.EV.BESS.energy_scale != 0.){
					end_user_profiles[point_iter][sample_iter].operation.default_demand_profile *= 1.;
				}
				end_user_profiles[point_iter][sample_iter].operation.default_demand_profile *= 1. - end_user_profiles[point_iter][sample_iter].investment.decision.smart_appliance * end_user_profiles[point_iter][sample_iter].operation.smart_appliance.scale;
				end_user_profiles[point_iter][sample_iter].operation.default_PV_profile = Power_network_inform.points.solar_cf.row(point_iter).segment(start_time, foresight_time);
				end_user_profiles[point_iter][sample_iter].operation.default_PV_profile *= end_user_profiles[point_iter][sample_iter].operation.PV_scale;
				end_user_profiles[point_iter][sample_iter].operation.PV_output = end_user_profiles[point_iter][sample_iter].operation.default_PV_profile(0);
				end_user_profiles[point_iter][sample_iter].operation.price_demand_profile = Power_market_inform.agent_profiles.aggregators[point_iter].price_demand_profile;
				end_user_profiles[point_iter][sample_iter].operation.price_supply_profile = Power_market_inform.agent_profiles.aggregators[point_iter].price_supply_profile;

				// Store demand characteristic of each node
				double demand_inflex;
				double demand_shiftable;
                demand_inflex = end_user_profiles[point_iter][sample_iter].operation.default_demand_profile(0);
                demand_inflex *= Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
                demand_inflex *= end_user_profiles[point_iter][sample_iter].operation.weight;
                Power_market_inform.TSO_Market.flex_stat_end.demand_inflex(start_time, node_ID) += demand_inflex;
                demand_shiftable = demand_inflex;
                demand_shiftable *= end_user_profiles[point_iter][sample_iter].investment.decision.smart_appliance * end_user_profiles[point_iter][sample_iter].operation.smart_appliance.scale;
                demand_shiftable /= 1. - end_user_profiles[point_iter][sample_iter].investment.decision.smart_appliance * end_user_profiles[point_iter][sample_iter].operation.smart_appliance.scale;
                demand_shiftable *= end_user_profiles[point_iter][sample_iter].operation.weight;
                Power_market_inform.TSO_Market.flex_stat_end.demand_shiftable(start_time, node_ID) += demand_shiftable;
//                if(Power_market_inform.agent_profiles.end_user_type.contingency[sample_iter] != 0){
//                }
//                else{
//                    demand_inflex = 0.;
//                }

				// Set the LP problem
				agent::end_user::end_user_LP_set(end_user_profiles[point_iter][sample_iter]);

				// Set bids and results information
				agent_bids_initialization(end_user_profiles[point_iter][sample_iter].operation.bids);
				agent_results_set(end_user_profiles[point_iter][sample_iter].operation.results);
				agent_settlement_set(end_user_profiles[point_iter][sample_iter].operation.settlement);

				// Totally inflexible end-user, demand profile as default
				double scale = end_user_profiles[point_iter][sample_iter].operation.weight;
				scale *= Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
				Power_market_inform.TSO_Market.flex_stat_end.demand_inflex(start_time, node_ID) -= demand_inflex;
				if(!end_user_profiles[point_iter][sample_iter].investment.decision.smart_management){
//					end_user_profiles[point_iter][sample_iter].operation.bids.submitted_demand_inflex(price_interval + 1) = Power_network_inform.points.nominal_mean_demand_field(point_iter, start_time);
//					end_user_profiles[point_iter][sample_iter].operation.bids.submitted_demand_inflex(price_interval + 1) *= agent::parameters::residential_ratio();
					//end_user_profiles[point_iter][sample_iter].operation.direct_demand = Power_network_inform.points.nominal_mean_demand_field(point_iter, start_time);
					agent::end_user::end_user_no_LP(start_time, end_user_profiles[point_iter][sample_iter], demand_inflex, scale, process_par);
				}
				else{
//					if(point_iter != 0){
//						continue;
//					}
					// Optimization and update process variables
					agent::end_user::end_user_LP_optimize(start_time, end_user_profiles[point_iter][sample_iter], process_par);
				}
				// Update inflexible demand
//				if(Power_market_inform.agent_profiles.end_user_type.contingency[sample_iter] == 0){
//                    demand_inflex = 0.;
//				}
				Power_market_inform.TSO_Market.flex_stat_end.demand_inflex(start_time, node_ID) += demand_inflex;

				// Scale the bids correctly
				agent_submitted_bids_scale(scale, end_user_profiles[point_iter][sample_iter].operation.bids);
			}
		}

		return end_user_profiles;
	}

	agent::industrial::profiles industrial_set(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int point_num = Power_network_inform.points.bidding_zone.size();
		int price_interval = power_market::parameters::price_interval();

		agent::industrial::profiles industrial_profiles;
		industrial_profiles.HV.reserve(point_num);
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int point_ID = point_iter;
			int node_ID = Power_network_inform.points.node(point_ID);
			agent::industrial::profile profile_temp;
			profile_temp.point_ID = point_ID;

			double bid_inflex_industrial = Power_network_inform.points.nominal_mean_demand_field(point_iter, tick);
			bid_inflex_industrial *= Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
			bid_inflex_industrial *= 1. - agent::parameters::residential_ratio();
			double bid_flex_industrial = bid_inflex_industrial;
			bid_inflex_industrial *= 1. - agent::industrial::parameters::flexible_ratio();
			bid_flex_industrial *= agent::industrial::parameters::flexible_ratio();
//			if(Power_market_inform.agent_profiles.end_user_type.contingency[0] != 0){
//
//			}
            Power_market_inform.TSO_Market.flex_stat_end.demand_inflex(tick, node_ID) += bid_inflex_industrial;
            Power_market_inform.TSO_Market.flex_stat_end.demand_flex(tick, node_ID) += bid_flex_industrial;
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

	agent::power_supplier::profiles power_supplier_set(int start_time, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int point_num = Power_network_inform.points.bidding_zone.size();
		int price_interval = power_market::parameters::price_interval();
		int hydro_num = Power_network_inform.plants.hydro.cap.size();
		int wind_num = Power_network_inform.plants.wind.cap.size();
		double cutoff_power = agent::power_supplier::parameters::cutoff_power();

		// Hydro and pump storage
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
				profile_temp.original_ID = agent_iter;
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
				Power_market_inform.TSO_Market.flex_stat.supply_flex(start_time, node_ID) += bid_vec.sum();

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
				profile_temp.original_ID = agent_iter;
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
				Power_market_inform.TSO_Market.flex_stat.supply_flex(start_time, node_ID) += bid_vec.sum();

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

		// Wind power plants
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
			double bid_quan = Power_network_inform.points.wind_on_cf(point_ID, start_time) * Power_network_inform.plants.wind.cap(agent_iter);

			agent::power_supplier::plant_profile profile_temp;
			profile_temp.point_ID = point_ID;
			profile_temp.original_ID = agent_iter;
			profile_temp.cap = Power_network_inform.plants.wind.cap(agent_iter);

			// Set bids information
			agent_bids_initialization(profile_temp.bids);
			agent_results_set(profile_temp.results);
			agent_settlement_set(profile_temp.settlement);
			profile_temp.bids.submitted_supply_flex(price_supply_flex_ID) = bid_quan;
			profile_temp.bids.redispatch_supply = profile_temp.bids.submitted_supply_flex;
			profile_temp.bids.filter_supply = profile_temp.bids.submitted_supply_flex;
			profile_temp.bids.balancing_supply = profile_temp.bids.submitted_supply_flex;
			Power_market_inform.TSO_Market.flex_stat.supply_flex(start_time, node_ID) += bid_quan;

			// High voltage power plants connect directly to transmission network
			if(profile_temp.cap >= cutoff_power){
				power_supplier_profiles.wind.HV_plant.push_back(profile_temp);
			}
			// Low voltage power plants feed into distribution network
			else{
				power_supplier_profiles.wind.LV_plant.push_back(profile_temp);
			}
		}

		// Slack power plants
		power_supplier_profiles.slack.LV_plant.reserve(point_num);
		for(int agent_iter = 0; agent_iter < point_num; ++ agent_iter){
			int point_ID = agent_iter;
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);
			double bid_quan = Power_network_inform.points.nominal_mean_demand_field(point_ID, start_time);
			//double bid_quan = 0.;
			bid_quan *= Power_network_inform.points.population_density(point_ID);
			bid_quan *= Power_network_inform.points.point_area / 1000.;
			//bid_quan *= 1.;
//			double bid_quan = std::numeric_limits<double>::infinity();

			agent::power_supplier::plant_profile profile_temp;
			profile_temp.point_ID = point_ID;
			profile_temp.original_ID = agent_iter;
			profile_temp.fix_cost = 0.;
			profile_temp.var_cost = profile_temp.fix_cost;
			profile_temp.var_cost += Power_market_inform.International_Market.confirmed.price(start_time, bz_ID);
			int price_supply_flex_ID = Power_market_inform.price_map.price_ID[profile_temp.var_cost];
			int price_supply_max_ID = Power_market_inform.price_map.price_ID[99.5];
			int price_length = std::max(price_supply_max_ID - price_supply_flex_ID + 1, 1); // avoid negative length
			//std::cout << price_supply_flex_ID << "\t" << profile_temp.var_cost << "\n";

			// Set bids information
			agent_bids_initialization(profile_temp.bids);
			agent_results_set(profile_temp.results);
			profile_temp.bids.submitted_supply_flex.segment(price_supply_flex_ID, price_length) = Eigen::VectorXd::Ones(price_length) * bid_quan / price_length;
			profile_temp.bids.redispatch_supply = profile_temp.bids.submitted_supply_flex;
			profile_temp.bids.filter_supply = profile_temp.bids.submitted_supply_flex;
			profile_temp.bids.balancing_supply = profile_temp.bids.submitted_supply_flex;
			power_supplier_profiles.slack.LV_plant.push_back(profile_temp);
		}

		return power_supplier_profiles;
	}

	void cross_border_redispatch_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int edge_num = Power_market_inform.agent_profiles.cross_border.size();

		for(int edge_iter = 0; edge_iter < edge_num; ++ edge_iter){
			int node_num = Power_market_inform.agent_profiles.cross_border[edge_iter].node_num;
			int entry_bz_ID = Power_market_inform.agent_profiles.cross_border[edge_iter].entry_bz_ID;

			// Store cross-border transmission flow as inflexible supply / demand at entry nodes
			if(node_num == 0){
				continue;
			}
			bool sink_flag = Power_market_inform.International_Market.network.confirmed_power(tick, edge_iter) >= 0.;
			int price_inflex_ID = 0 + sink_flag * (Power_market_inform.International_Market.price_intervals - 0);
			double price_flex = Power_market_inform.International_Market.confirmed.price(tick, entry_bz_ID);
			int price_flex_ID = Power_market_inform.price_map.price_ID[price_flex];
			double source = -(1 - sink_flag) * Power_market_inform.International_Market.network.confirmed_power(tick, edge_iter);
			double sink = sink_flag * Power_market_inform.International_Market.network.confirmed_power(tick, edge_iter);

			for(int node_iter = 0; node_iter < node_num; ++ node_iter){
				int node_ID = Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].node_ID;
				int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);

				agent_bids_initialization(Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].bids);
				agent_results_set(Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results);

				Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results.cleared_supply = source / node_num;
				Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results.cleared_demand = sink / node_num;
				Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].bids.submitted_supply_flex(price_flex_ID) = source / node_num;
				Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].bids.submitted_demand_flex(price_flex_ID) = sink / node_num;
				Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].bids.redispatch_supply(price_inflex_ID) = source / node_num;
				Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].bids.redispatch_demand(price_inflex_ID) = sink / node_num;
				Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].bids.balancing_supply(price_inflex_ID) = source / node_num;
				Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].bids.balancing_demand(price_inflex_ID) = sink / node_num;

				// Settlement in EOM
				double original_price = Power_market_inform.International_Market.confirmed.price(tick, bz_ID);
				agent_EOM_settlement_calculation(tick, node_ID, original_price, original_price, Power_market_inform, Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].bids, Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results, Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].settlement);
			}
		}
	}

	void end_user_redispatch_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, configuration::process_config &process_par){
		int point_num = Power_network_inform.points.bidding_zone.size();
		//int sample_num = agent::end_user::parameters::sample_num();
		int sample_num = Power_market_inform.agent_profiles.end_user_type.sample_num;
		int price_interval = power_market::parameters::price_interval();
		power_market::parameters::price_ID_bimap bidded_price_map;
		power_market::parameters::bidded_price(bidded_price_map);

		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int bz_ID = Power_network_inform.points.bidding_zone(point_iter);
			double marginal_price = Power_market_inform.International_Market.confirmed.price(tick, bz_ID);
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
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(0) += (1. - Power_market_inform.International_Market.confirmed.ratio_demand(bz_ID)) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_inflex(marginal_price_ID);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(0) += Power_market_inform.International_Market.confirmed.ratio_supply(bz_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_inflex(marginal_price_ID);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(price_interval + 1) += Power_market_inform.International_Market.confirmed.ratio_demand(bz_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_inflex(marginal_price_ID);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(price_interval + 1) += (1. - Power_market_inform.International_Market.confirmed.ratio_supply(bz_ID)) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_inflex(marginal_price_ID);

				// Flexible bids have redispatcch priority in between
				if(process_par.encourage_redispatch){
					int price_supply_flex_BESS_ID = bidded_price_map.price_ID[Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.price_supply];
					int price_demand_flex_BESS_ID = bidded_price_map.price_ID[Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.price_demand];
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(price_demand_flex_BESS_ID) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_flex.sum();
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(price_supply_flex_BESS_ID) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_flex.sum();
					//Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_flex;
					//Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_flex;
				}
				else{
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_flex;
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_flex;
				}

				// Update filter bids
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.filter_demand = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.filter_supply = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply;
			}
		}
	}

	void end_user_filter_demand_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int point_num = Power_network_inform.points.bidding_zone.size();
		//int sample_num = agent::end_user::parameters::sample_num();
		int sample_num = Power_market_inform.agent_profiles.end_user_type.sample_num;
		int price_interval = power_market::parameters::price_interval();

		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int point_ID = Power_network_inform.points.in_cluster_ID(point_iter);
			int node_ID = Power_network_inform.points.node(point_iter);
			int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			double marginal_price = Power_market_inform.DSO_Markets[DSO_ID].confirmed.price(tick, point_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
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

			if(marginal_price_ID > 0){
				Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_demand.head(marginal_price_ID) *= 0.;
			}
			Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_demand(marginal_price_ID) *= Power_market_inform.DSO_Markets[DSO_ID].confirmed.ratio_demand(point_ID);
			Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.balancing_demand = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_demand;
		}
	}

	void end_user_filter_supply_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int point_num = Power_network_inform.points.bidding_zone.size();
		//int sample_num = agent::end_user::parameters::sample_num();
		int sample_num = Power_market_inform.agent_profiles.end_user_type.sample_num;
		int price_interval = power_market::parameters::price_interval();

		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int point_ID = Power_network_inform.points.in_cluster_ID(point_iter);
			int node_ID = Power_network_inform.points.node(point_iter);
			int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			double marginal_price = Power_market_inform.DSO_Markets[DSO_ID].confirmed.price(tick, point_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
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

			if(marginal_price_ID < price_interval + 1){
				Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_supply.tail(price_interval + 1 - marginal_price_ID) *= 0.;
			}
			Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_supply(marginal_price_ID) *= Power_market_inform.DSO_Markets[DSO_ID].confirmed.ratio_supply(point_ID);
			Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.balancing_supply = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_supply;
		}

		int slack_LV_num = Power_market_inform.agent_profiles.power_supplier.slack.LV_plant.size();
		for(int agent_iter = 0; agent_iter < slack_LV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int DSO_ID = Power_network_inform.nodes.cluster(node_ID);
			point_ID = Power_network_inform.points.in_cluster_ID(point_ID);
			double marginal_price = Power_market_inform.DSO_Markets[DSO_ID].confirmed.price(tick, point_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			if(marginal_price_ID < price_interval + 1){
				Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].bids.redispatch_supply.tail(price_interval + 1 - marginal_price_ID) *= 0.;
			}
			Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].bids.redispatch_supply(marginal_price_ID) *= Power_market_inform.DSO_Markets[DSO_ID].confirmed.ratio_supply(point_ID);
			Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].bids.balancing_supply = Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].bids.redispatch_supply;
		}
	}

	void cross_border_balancing_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int edge_num = Power_market_inform.agent_profiles.cross_border.size();

		for(int edge_iter = 0; edge_iter < edge_num; ++ edge_iter){
			int node_num = Power_market_inform.agent_profiles.cross_border[edge_iter].node_num;
			int entry_bz_ID = Power_market_inform.agent_profiles.cross_border[edge_iter].entry_bz_ID;

			if(node_num == 0){
				continue;
			}
			for(int node_iter = 0; node_iter < node_num; ++ node_iter){
				int node_ID = Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].node_ID;
				int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);
				double marginal_price = Power_market_inform.TSO_Market.confirmed.price(tick, node_ID);
				int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];
				double original_price = Power_market_inform.International_Market.confirmed.price(tick, bz_ID);
				int original_price_ID = Power_market_inform.price_map.price_ID[original_price];

				// Calculate scheduled results
				agent_scheduled_results_calculation(bz_ID, node_ID, marginal_price_ID, original_price_ID, Power_market_inform, Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].bids, Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results, 0);

				// Settlement of redispatch
				agent_redispatch_settlement_calculation(tick, node_ID, original_price, Power_market_inform, Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].bids, Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results, Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].settlement);
			}
		}
	}

	void end_user_balancing_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int point_num = Power_network_inform.points.bidding_zone.size();
		//int sample_num = agent::end_user::parameters::sample_num();
		int sample_num = Power_market_inform.agent_profiles.end_user_type.sample_num;
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
				// Calculate scheduled results
				agent_scheduled_results_calculation(bz_ID, node_ID, marginal_price_ID, original_price_ID, Power_market_inform, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results);

				if(!Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].investment.decision.control_reserve){
					if(marginal_price_ID > 0){
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand(0) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand.head(marginal_price_ID).sum();
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply(0) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply.head(marginal_price_ID).sum();
					}
					if(marginal_price_ID < price_interval + 1){
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand(price_interval + 1) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand.tail(price_interval + 1 - marginal_price_ID).sum();
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply(price_interval + 1) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply.tail(price_interval + 1 - marginal_price_ID).sum();
					}

					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand(0) += (1. - Power_market_inform.TSO_Market.confirmed.ratio_demand(node_ID)) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(marginal_price_ID);
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply(0) += Power_market_inform.TSO_Market.confirmed.ratio_supply(node_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(marginal_price_ID);
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand(price_interval + 1) += Power_market_inform.TSO_Market.confirmed.ratio_demand(node_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand(marginal_price_ID);
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply(price_interval + 1) += (1. - Power_market_inform.TSO_Market.confirmed.ratio_supply(node_ID)) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply(marginal_price_ID);
				}
				else{
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand;
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply;
				}

				// Imbalance accounting
				double imbalance = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.confirmed_demand;
				imbalance *= Power_network_inform.points.imbalance_field(point_iter, tick);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.imbalance_demand += imbalance;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement.volume_demand_up.imbalance += (imbalance > 0.) * imbalance;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement.volume_demand_down.imbalance += (imbalance < 0.) * imbalance;
				Power_market_inform.TSO_Market.imbalance.demand_up(tick, node_ID) += (imbalance > 0.) * imbalance;
				Power_market_inform.TSO_Market.imbalance.demand_down(tick, node_ID) += (imbalance < 0.) * imbalance;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.imbalance_demand = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.imbalance_demand *= Power_network_inform.points.imbalance_field(point_iter, tick);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.imbalance_demand;

				// Settlement in EOM
				// Aggregators
//				Power_market_inform.agent_profiles.aggregators[point_iter].settlement.utility_supply.EOM += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_demand * Power_market_inform.agent_profiles.aggregators[point_iter].price_demand_profile(0);
//				Power_market_inform.agent_profiles.aggregators[point_iter].settlement.utility_supply.EOM += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_supply * original_price;
//				Power_market_inform.agent_profiles.aggregators[point_iter].settlement.price.EOM += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_demand * original_price;
//				Power_market_inform.agent_profiles.aggregators[point_iter].settlement.price.EOM += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_supply * Power_market_inform.agent_profiles.aggregators[point_iter].price_supply_profile(0);
				// End-users
//				agent_EOM_settlement_calculation(tick, node_ID, Power_market_inform.agent_profiles.aggregators[point_iter].price_supply_profile(0), Power_market_inform.agent_profiles.aggregators[point_iter].price_demand_profile(0), Power_market_inform, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement, 1);
				agent_EOM_settlement_calculation(tick, node_ID, original_price, original_price, Power_market_inform, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement, 1);
//				double over_est_utility_EOM = (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_demand > Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.direct_demand);
//				over_est_utility_EOM *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_demand - Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.direct_demand;
//				over_est_utility_EOM *= Power_market_inform.price_map.bidded_price(price_interval + 1);
//				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement.utility_demand.EOM -= over_est_utility_EOM;
//				Power_market_inform.TSO_Market.EOM.utility(tick, node_ID) -= over_est_utility_EOM;
				// Utility from EV
//				double scale = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.weight;
//				scale *= agent::parameters::residential_ratio();
//				scale *= Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
//				double under_est_utility_EOM = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.default_demand_profile(0);
//				under_est_utility_EOM *= Power_market_inform.price_map.bidded_price(price_interval + 1) * scale;
//				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement.utility_demand.EOM += under_est_utility_EOM;
//				Power_market_inform.TSO_Market.EOM.utility(tick, node_ID) += under_est_utility_EOM;

				// Settlement of redispatch
				agent_redispatch_settlement_calculation(tick, node_ID, original_price, Power_market_inform, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement, 1);
//				if(point_iter == 3818 && sample_iter == 0){
//					std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_demand << "\t";
//					std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.confirmed_demand << "\n";
//				}
				//end_user_redispatch_settlement_calculation(tick, node_ID, original_price, Power_market_inform, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation, Power_market_inform.agent_profiles.aggregators[point_iter], 1);
			}
		}
		//std::cout << "\n";
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
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].results.imbalance_demand += imbalance;
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement.volume_demand_up.imbalance += (imbalance > 0.) * imbalance;
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement.volume_demand_down.imbalance += (imbalance < 0.) * imbalance;
			Power_market_inform.TSO_Market.imbalance.demand_up(tick, node_ID) += (imbalance > 0.) * imbalance;
			Power_market_inform.TSO_Market.imbalance.demand_down(tick, node_ID) += (imbalance < 0.) * imbalance;
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids.imbalance_demand = Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids.balancing_demand;
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids.imbalance_demand *= Power_network_inform.points.imbalance_field(point_ID, tick);
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids.balancing_demand += Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids.imbalance_demand;

			// Settlement in EOM
			agent_EOM_settlement_calculation(tick, node_ID, original_price, original_price, Power_market_inform, Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids, Power_market_inform.agent_profiles.industrial.HV[agent_iter].results, Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement);

			// Settlement of redispatch
			agent_redispatch_settlement_calculation(tick, node_ID, original_price, Power_market_inform, Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids, Power_market_inform.agent_profiles.industrial.HV[agent_iter].results, Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement);
		}
		//std::cout << "\n";
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
			agent_EOM_settlement_calculation(tick, node_ID, original_price, original_price, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].settlement, 0);

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
			agent_EOM_settlement_calculation(tick, node_ID, original_price, original_price, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].settlement);

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
			agent_EOM_settlement_calculation(tick, node_ID, original_price, original_price, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].settlement);

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
			agent_EOM_settlement_calculation(tick, node_ID, original_price, original_price, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].settlement);

			// Settlement of redispatch
			agent_redispatch_settlement_calculation(tick, node_ID, original_price, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].settlement);
		}

		int pump_HV_plant_num =  Power_market_inform.agent_profiles.power_supplier.pump_storage.HV.size();
		for(int agent_iter = 0; agent_iter < pump_HV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
			double marginal_price = Power_market_inform.TSO_Market.confirmed.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];
			double original_price = Power_market_inform.International_Market.confirmed.price(tick, bz_ID);
			int original_price_ID = Power_market_inform.price_map.price_ID[original_price];

			// Calculate scheduled results
			agent_scheduled_results_calculation(bz_ID, node_ID, marginal_price_ID, original_price_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].results);

			// Settlement in EOM
			agent_EOM_settlement_calculation(tick, node_ID, original_price, original_price, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].settlement);

			// Settlement of redispatch
			agent_redispatch_settlement_calculation(tick, node_ID, original_price, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].settlement);
		}

		int pump_LV_plant_num =  Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
		for(int agent_iter = 0; agent_iter < pump_LV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
			double marginal_price = Power_market_inform.TSO_Market.confirmed.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];
			double original_price = Power_market_inform.International_Market.confirmed.price(tick, bz_ID);
			int original_price_ID = Power_market_inform.price_map.price_ID[original_price];

			// Calculate scheduled results
			agent_scheduled_results_calculation(bz_ID, node_ID, marginal_price_ID, original_price_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].results);

			// Settlement in EOM
			agent_EOM_settlement_calculation(tick, node_ID, original_price, original_price, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].settlement);

			// Settlement of redispatch
			agent_redispatch_settlement_calculation(tick, node_ID, original_price, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].settlement);
		}

		int slack_LV_num = Power_market_inform.agent_profiles.power_supplier.slack.LV_plant.size();
		for(int agent_iter = 0; agent_iter < slack_LV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
			double marginal_price = Power_market_inform.TSO_Market.confirmed.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];
//			double original_price = Power_market_inform.International_Market.confirmed.price(tick, bz_ID);
//			int original_price_ID = Power_market_inform.price_map.price_ID[original_price];
			int original_price_ID = 0;
			double original_price = Power_market_inform.price_map.bidded_price[original_price_ID];

			// Calculate scheduled results
			agent_scheduled_results_calculation(bz_ID, node_ID, marginal_price_ID, original_price_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].results);
            Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].results.cleared_supply = 0.;
            Power_market_inform.TSO_Market.flex_stat.supply_flex(tick, node_ID) += Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].results.confirmed_supply;

			// Settlement in EOM
			//agent_EOM_settlement_calculation(tick, node_ID, original_price, original_price, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].settlement);

			// Settlement of redispatch
			agent_redispatch_settlement_calculation(tick, node_ID, original_price, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].settlement);
		}

		//std::cout << "\n";
	}

	void agents_redispatch_settlement(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int node_num = Power_market_inform.TSO_Market.network.num_vertice;
		int bz_num = Power_market_inform.International_Market.cross_border_zone_start;
		Eigen::VectorXd redispatch_confirmed_demand = Eigen::VectorXd::Zero(bz_num);
		Eigen::VectorXd redispatch_confirmed_supply = Eigen::VectorXd::Zero(bz_num);
		Eigen::VectorXd redispatch_cost_supply_self = Eigen::VectorXd::Zero(bz_num);
		Eigen::VectorXd redispatch_cost_supply_export = Eigen::VectorXd::Zero(bz_num);
		Eigen::VectorXd redispatch_quan_supply_import = Eigen::VectorXd::Zero(bz_num);

		// Redispatch price per energy at each transmission node
		for(int node_iter = 0; node_iter < node_num; ++ node_iter){
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_iter);

			Power_market_inform.International_Market.redispatch.price_demand(tick, bz_ID) += Power_market_inform.TSO_Market.redispatch.price_demand(tick, node_iter);
			redispatch_confirmed_demand(bz_ID) += Power_market_inform.TSO_Market.confirmed.demand(tick, node_iter);

			Power_market_inform.International_Market.redispatch.price_supply(tick, bz_ID) += Power_market_inform.TSO_Market.redispatch.price_supply(tick, node_iter);
			redispatch_confirmed_supply(bz_ID) += Power_market_inform.TSO_Market.confirmed.supply(tick, node_iter);
		}

		for(int zone_iter = 0; zone_iter < bz_num; ++ zone_iter){
			redispatch_cost_supply_self(zone_iter) = Power_market_inform.International_Market.redispatch.price_supply(tick, zone_iter) * std::min(redispatch_confirmed_demand(zone_iter), redispatch_confirmed_supply(zone_iter)) / redispatch_confirmed_supply(zone_iter);
			redispatch_quan_supply_import(zone_iter) = redispatch_confirmed_demand(zone_iter) - std::min(redispatch_confirmed_demand(zone_iter), redispatch_confirmed_supply(zone_iter));
			redispatch_cost_supply_export(zone_iter) = Power_market_inform.International_Market.redispatch.price_supply(tick, zone_iter) - redispatch_cost_supply_self(zone_iter);
		}

		for(int zone_iter = 0; zone_iter < bz_num; ++ zone_iter){
			Power_market_inform.International_Market.redispatch.price_demand(tick, zone_iter) += redispatch_cost_supply_self(zone_iter);
			Power_market_inform.International_Market.redispatch.price_demand(tick, zone_iter) += redispatch_quan_supply_import(zone_iter) / redispatch_quan_supply_import.sum() * redispatch_cost_supply_export.sum();
			Power_market_inform.International_Market.redispatch.price_demand(tick, zone_iter) /= std::min(Power_market_inform.International_Market.confirmed.demand(tick, zone_iter), redispatch_confirmed_demand(zone_iter));
		}
		//std::cout << Power_market_inform.International_Market.redispatch.price_demand.row(tick) << "\n\n";

		// End-users
		int point_num = Power_network_inform.points.bidding_zone.size();
		//int sample_num = agent::end_user::parameters::sample_num();
		int sample_num = Power_market_inform.agent_profiles.end_user_type.sample_num;
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int bz_ID = Power_network_inform.points.bidding_zone(point_iter);

			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				double redispatch_price =  Power_market_inform.International_Market.redispatch.price_demand(tick, bz_ID);
				redispatch_price *= std::min(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_demand, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.confirmed_demand);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement.price.redispatch += redispatch_price;
			}
		}

		// Industrial demand
		int industrial_HV_num = Power_market_inform.agent_profiles.industrial.HV.size();
		for(int agent_iter = 0; agent_iter < industrial_HV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.industrial.HV[agent_iter].point_ID;
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
			double redispatch_price =  Power_market_inform.International_Market.redispatch.price_demand(tick, bz_ID);
			redispatch_price *= std::min(Power_market_inform.agent_profiles.industrial.HV[agent_iter].results.cleared_demand, Power_market_inform.agent_profiles.industrial.HV[agent_iter].results.confirmed_demand);
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement.price.redispatch += redispatch_price;
		}

		// Pump hydro plants
		int pump_HV_plant_num =  Power_market_inform.agent_profiles.power_supplier.pump_storage.HV.size();
		for(int agent_iter = 0; agent_iter < pump_HV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].point_ID;
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
			double redispatch_price =  Power_market_inform.International_Market.redispatch.price_demand(tick, bz_ID);
			redispatch_price *= std::min(Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].results.cleared_demand, Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].results.confirmed_demand);
			Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].settlement.price.redispatch += redispatch_price;
		}

		int pump_LV_plant_num =  Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
		for(int agent_iter = 0; agent_iter < pump_LV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
			double redispatch_price =  Power_market_inform.International_Market.redispatch.price_demand(tick, bz_ID);
			redispatch_price *= std::min(Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].results.cleared_demand, Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].results.confirmed_demand);
			Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].settlement.price.redispatch += redispatch_price;
		}

		// Cross-border flow
		int edge_num = Power_market_inform.agent_profiles.cross_border.size();
		for(int edge_iter = 0; edge_iter < edge_num; ++ edge_iter){
			int node_num = Power_market_inform.agent_profiles.cross_border[edge_iter].node_num;
			int entry_bz_ID = Power_market_inform.agent_profiles.cross_border[edge_iter].entry_bz_ID;

			// Store cross-border transmission flow as inflexible supply / demand at entry nodes
			if(node_num == 0){
				continue;
			}

			for(int node_iter = 0; node_iter < node_num; ++ node_iter){
				int node_ID = Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].node_ID;
				int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);

				double redispatch_price =  Power_market_inform.International_Market.redispatch.price_demand(tick, bz_ID);
				redispatch_price *= std::min(Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results.cleared_demand, Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results.confirmed_demand);
				Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].settlement.price.redispatch += redispatch_price;
			}
		}
	}

	void cross_border_status_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, bool control_reserve_flag){
		int edge_num = Power_market_inform.agent_profiles.cross_border.size();

		for(int edge_iter = 0; edge_iter < edge_num; ++ edge_iter){
			int node_num = Power_market_inform.agent_profiles.cross_border[edge_iter].node_num;
			int entry_bz_ID = Power_market_inform.agent_profiles.cross_border[edge_iter].entry_bz_ID;

			if(node_num == 0){
				continue;
			}
			for(int node_iter = 0; node_iter < node_num; ++ node_iter){
				int node_ID = Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].node_ID;
				int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);
				double marginal_price = Power_market_inform.TSO_Market.actual.price(tick, node_ID);
				int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

				// Calculate actual results
				agent_actual_results_calculation(node_ID, marginal_price_ID, Power_market_inform, Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].bids, Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results, control_reserve_flag);

				// Balancing settlement
				agent_balancing_settlement_calculation(tick, node_ID, Power_market_inform, Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].bids, Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results, Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].settlement);

				// Market operation update
				market_operation_update(tick, bz_ID, Power_market_inform.International_Market.operation.cross_border, Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results);
			}
		}
	}

	void end_user_status_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, configuration::process_config &process_par){
		int point_num = Power_network_inform.points.bidding_zone.size();
		int sample_num = Power_market_inform.agent_profiles.end_user_type.sample_num;
		int load_shift_time = agent::end_user::parameters::load_shift_time();
		int price_interval = power_market::parameters::price_interval();
		int foresight_time = agent::end_user::parameters::foresight_time();
		power_market::parameters::price_ID_bimap bidded_price_map;
		power_market::parameters::bidded_price(bidded_price_map);

		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int node_ID = Power_network_inform.points.node(point_iter);
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);
			double marginal_price = Power_market_inform.TSO_Market.actual.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

//			if(point_iter != 0){
//				continue;
//			}

			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				if(process_par.control_reserve_flag){
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

				// Totally inflexible end-user don't need to update status
				if(!Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].investment.decision.smart_management){
					// Balancing settlement
					agent_balancing_settlement_calculation(tick, node_ID, Power_market_inform, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].investment.decision.control_reserve);

					// Market operation update
					market_operation_update(tick, bz_ID, Power_market_inform.International_Market.operation.end_user, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results);

					if(point_iter == 0 && sample_iter >= 2){
						std::cout << sample_iter << "\t";
						//std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand.transpose() << "\n";
						std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity << "\t";
						std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.soc << "\t";
						std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity << "\t";
						std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc << "\n\n";
					}

					continue;
				}

				// New
				double soc_min_EV = (tick % foresight_time >= 3) * (tick % foresight_time <= 6) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.energy_demand * (tick % foresight_time - 2) / 4;
				double scale = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.weight;
				scale *= Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
				double tol = 1E-12;

				double demand_remain = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_demand;
				demand_remain /= scale;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity = 0.;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity = 0.;
//				if(sample_iter == 2){
//					std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.default_demand_profile(0) << "\t" << demand_remain << "\n";
//				}

				// Fulfill default demand first
				demand_remain -= std::min(demand_remain,  Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.default_demand_profile(0));

				// Fulfill demand from highest marginal benefits
				for(int price_iter = price_interval + 1; price_iter >= 0; -- price_iter){
					// Fulfill smart appliance first, from earliest to latest
					for(int tock = 0; tock < 2 * load_shift_time + 1; ++ tock){
						if(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.price_demand(tock) == bidded_price_map.bidded_price(price_iter)){
							double sa_flex = std::min(demand_remain, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand(tock));
//							if(sample_iter == 2){
//								std::cout << bidded_price_map.bidded_price(price_iter) << "\t" << demand_remain << "\t" << sa_flex << "\n";
//							}
							Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand(tock) -= sa_flex;
							demand_remain -= sa_flex;
						}
					}

					// Fulfill EV demand second
                    double EV_flex = std::min(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.capacity_scale, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.energy_scale - Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc + Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.self_consumption);
                    EV_flex *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.house_default_period(0);
                    EV_flex /= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency;

                    double EV_inflex = std::max(soc_min_EV - Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc + Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.self_consumption, 0.);
                    EV_inflex *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.house_default_period(0);
                    EV_inflex /= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency;

                    // Inflexible EV demand has priority
					if(price_iter == price_interval + 1){
                        EV_inflex = std::min(EV_inflex, demand_remain);
                        demand_remain -= EV_inflex;
                        EV_inflex *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency;
                        Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity += EV_inflex;
                        Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc += EV_inflex;

//                        if(point_iter == 0 && sample_iter >= 2){
//                            std::cout << sample_iter << "\t";
//                            std::cout << EV_inflex << "\t";
//                        }
					}

					//if not encourage redispatch mode dispatch flexible EV demand
					if(!process_par.encourage_redispatch){
						if(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.price_demand == bidded_price_map.bidded_price(price_iter)){
							double EV_flex = std::min(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.capacity_scale, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.energy_scale - Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc + Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.self_consumption);
							EV_flex *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.house_default_period(0);
							EV_flex /= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency;
							//std::cout << price_iter << "\t" << EV_flex << "\t" << demand_remain << "\n";
							EV_flex = std::min(EV_flex, demand_remain);
							demand_remain -= EV_flex;
							EV_flex *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency;
							Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity += EV_flex;
							Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc += EV_flex;
	//						if(sample_iter == 2){
	//							std::cout << bidded_price_map.bidded_price(price_iter) << "\t" << demand_remain << "\t" << EV_flex << "\n";
	//						}
						}
					}

					// Fulfill BESS demand last if not encourage redispatch mode
					if(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.price_demand == bidded_price_map.bidded_price(price_iter)){
						double BESS_flex = std::min(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.capacity_scale, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.energy_scale - Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.soc + Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.self_consumption);
						BESS_flex /= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency;
						BESS_flex = std::min(BESS_flex, demand_remain);
						demand_remain -= BESS_flex;
						BESS_flex *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency;
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity += BESS_flex;
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.soc += BESS_flex;
//						if(sample_iter == 2){
//							std::cout << bidded_price_map.bidded_price(price_iter) << "\t" << demand_remain << "\t" << BESS_flex << "\n";
//						}
					}

					// Fulfill EV demand last if encourage redispatch mode
					if(process_par.encourage_redispatch){
						if(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.price_demand == bidded_price_map.bidded_price(price_iter)){
							double EV_flex = std::min(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.capacity_scale, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.energy_scale - Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc + Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.self_consumption);
							EV_flex *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.house_default_period(0);
							EV_flex /= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency;
							//std::cout << price_iter << "\t" << EV_flex << "\t" << demand_remain << "\n";
							EV_flex = std::min(EV_flex, 1.);
							EV_flex = std::min(EV_flex, demand_remain);
							demand_remain -= EV_flex;
							EV_flex *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency;
							Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity += EV_flex;
							Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc += EV_flex;
	//						if(sample_iter == 2){
	//							std::cout << bidded_price_map.bidded_price(price_iter) << "\t" << demand_remain << "\t" << EV_flex << "\n";
	//						}
						}
					}

					if(demand_remain < tol){
						break;
					}
				}

				double supply_remain = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_supply;
				supply_remain /= scale;

				// Fulfill supply from lowest marginal costs
				for(int price_iter = 0; price_iter < price_interval + 1; ++ price_iter){
					// Fulfill PV supply first

					// Fulfill BESS supply second
					if(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.price_demand == bidded_price_map.bidded_price(price_iter)){
						double BESS_flex = std::min(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.capacity_scale, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.soc - Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.self_consumption);
						BESS_flex *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency;
						BESS_flex = std::min(BESS_flex, supply_remain);
						supply_remain -= BESS_flex;
						BESS_flex /= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency;
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity -= BESS_flex;
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.soc -= BESS_flex;
//						if(sample_iter == 2){
//							std::cout << bidded_price_map.bidded_price(price_iter) << "\t" << supply_remain << "\t" << BESS_flex << "\n";
//						}
					}

					// Fulfill EV supply last
					if(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.price_demand == bidded_price_map.bidded_price(price_iter)){
						double EV_flex = std::min(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.capacity_scale, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc - Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.self_consumption - soc_min_EV);
						EV_flex = std::max(EV_flex, 0.);
						EV_flex *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.house_default_period(0);
						EV_flex *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency;
						EV_flex = std::min(EV_flex, supply_remain);
						supply_remain -= EV_flex;
						EV_flex /= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency;
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity -= EV_flex;
						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc -= EV_flex;
//						if(sample_iter == 2){
//							std::cout << bidded_price_map.bidded_price(price_iter) << "\t" << supply_remain << "\t" << EV_flex << "\n";
//						}
					}

					if(supply_remain < tol){
						break;
					}
				}
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.soc -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.self_consumption;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.self_consumption;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.default_demand_profile(0);

//
////				// Original
////				// Should comment this section after validate
//////				double gap = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_demand;
//////				gap -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_supply;
//////				gap -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_inflex.sum();
//////				gap -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_flex.sum();
//////				gap += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_inflex.sum();
//////				gap += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_flex.sum();
//////				gap /= Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
//////				gap /= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.weight;
//////				gap *= -1.;
//////				if(point_iter == 1111){
//////					if(sample_iter == 0 || sample_iter == 2){
//////						std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_demand << "\t";
//////						std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_supply << "\t";
//////						std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_inflex.sum() << "\t";
//////						std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_flex.sum() << "\t";
//////						std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_inflex.sum() << "\t";
//////						std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_flex.sum() << "\t";
//////						std::cout << gap << "\n";
//////					}
//////				}
//////				// Should comment this section after validate
//////
//////
////
////				// Should uncomment this section after validate
////				double gap = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_demand;
////				gap -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_supply;
////				gap -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_inflex.sum();
////				gap -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_demand_flex.sum();
////				gap += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_inflex.sum();
////				gap += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.submitted_supply_flex.sum();
////				gap /= Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
////				gap /= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.weight;
////				gap *= -1.;
////				// Should uncomment this section after validate
////
////				std::cout << point_iter << "\t" << sample_iter << ":\t";
////				std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_supply << "\t";
////				std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.cleared_demand << "\t";
////				std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity << "\t";
////				std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.default_demand_profile(0) << "\n";
////				if(point_iter == 1111){
////					std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.price_demand_profile.head(3).transpose() << "\n";
////					std::cout << Power_market_inform.agent_profiles.end_users[point_iter][0].operation.default_demand_profile.head(3).transpose() << "\n";
////					std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.scheduled_demand.transpose() << "\n";
////					std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity << "\n";
////				}
////
////				// Actual demand smaller than initially planned
////				while(gap > 0.){
////					// Reduce BESS charge from scheduled
////					double BESS_flex_lb = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.soc - Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.self_consumption;
////					BESS_flex_lb = -std::min(BESS_flex_lb, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.capacity_scale);
////					double BESS_flex = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity;
////					if(BESS_flex > 0.){
////						BESS_flex = std::max(0., BESS_flex - gap * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency);
////						BESS_flex = std::max(BESS_flex, BESS_flex_lb);
////						gap -= (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity - BESS_flex) / Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency;
////						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity = BESS_flex;
////						if(gap == 0.){
////							break;
////						}
////					}
////					if(BESS_flex > BESS_flex_lb){
////						BESS_flex = std::max(BESS_flex_lb, BESS_flex - gap / Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency);
////						gap -= (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity - BESS_flex) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency;
////						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity = BESS_flex;
////						if(gap == 0.){
////							break;
////						}
////					}
////
////					// Reduce EV charge from scheduled
////					double EV_flex_lb = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc - Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.self_consumption;
////					EV_flex_lb -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.default_demand_profile(0);
////					EV_flex_lb = -std::min(EV_flex_lb, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.capacity_scale);
////					double EV_flex = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity;
////					if(EV_flex > 0.){
////						EV_flex = std::max(0., EV_flex - gap * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency);
////						EV_flex = std::max(EV_flex, EV_flex_lb);
////						gap -= (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity - EV_flex) / Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency;
////						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity = EV_flex;
////						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.house_default_period(0);
////						if(gap == 0.){
////							break;
////						}
////					}
////					if(EV_flex > EV_flex_lb){
////						EV_flex = std::max(EV_flex_lb, EV_flex - gap / Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency);
////						gap -= (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity - EV_flex) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency;
////						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity = EV_flex;
////						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.house_default_period(0);
////						if(gap == 0.){
////							break;
////						}
////					}
////
////					// Reduce smart appliance demand from scheduled
////					for(int tock = 0; tock < 2 * load_shift_time + 1; ++ tock){
////						int tock_ID = 2 * load_shift_time - tock;
////						double sa_flex = std::max(0., Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.scheduled_demand(tock_ID) - gap);
////						gap -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.scheduled_demand(tock_ID) - sa_flex;
////						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.scheduled_demand(tock_ID) = sa_flex;
////						if(gap == 0.){
////							break;
////						}
////					}
////					break;
////				}
////
////				// Actual demand greater than initially planned
////				while(gap < 0.){
////					// Increase BESS charge from scheduled
////					double BESS_flex_ub = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.energy_scale;
////					BESS_flex_ub -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.soc;
////					BESS_flex_ub += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.self_consumption;
////					BESS_flex_ub = std::min(BESS_flex_ub, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.capacity_scale);
////					double BESS_flex = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity;
////					if(BESS_flex < 0.){
////						BESS_flex = std::min(0., BESS_flex - gap / Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency);
////						gap -= (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity - BESS_flex) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency;
////						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity = BESS_flex;
////						if(gap == 0.){
////							break;
////						}
////					}
////					BESS_flex = std::min(BESS_flex_ub, BESS_flex - gap * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency);
////					gap -= (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity - BESS_flex) / Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency;
////					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity = BESS_flex;
////					if(gap == 0.){
////						break;
////					}
////
////					// Increase EV charge from scheduled
////					double EV_flex_ub = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.energy_scale;
////					EV_flex_ub -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc;
////					EV_flex_ub += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.self_consumption;
////					EV_flex_ub += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.default_demand_profile(0);
////					EV_flex_ub = std::min(EV_flex_ub, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.capacity_scale);
////					double EV_flex = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity;
////					if(EV_flex < 0.){
////						EV_flex = std::min(0., EV_flex - gap / Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency);
////						gap -= (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity - EV_flex) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency;
////						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity = EV_flex;
////						Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.house_default_period(0);
////						if(gap == 0.){
////							break;
////						}
////					}
////					EV_flex = std::min(EV_flex_ub, EV_flex - gap * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency);
////					gap -= (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity - EV_flex) / Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.efficiency;
////					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity = EV_flex;
////					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.house_default_period(0);
////					if(gap == 0.){
////						break;
////					}
////
////					// Decrease PV output from scheduled
////					double PV_flex = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.PV_output;
////					PV_flex = std::min(PV_flex, -gap);
////					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.PV_output -= PV_flex;
////					break;
////				}
//
//				// Update state variables of end-users
//				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.soc += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity;
//				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.soc -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.self_consumption;
//				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity;
//				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.self_consumption;
//				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.default_demand_profile(0);
//				for(int tock = 0; tock < 2 * load_shift_time + 1; ++ tock){
//					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand(tock) -= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.scheduled_demand(tock);
//				}

				if(point_iter == 0 && sample_iter >= 2){
					std::cout << sample_iter << "\t";
					//std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand.transpose() << "\n";
					std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity << "\t";
					std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.soc << "\t";
					std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity << "\t";
					std::cout << Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.soc << "\n\n";
				}

				// Update storage settlement
				double vol_ch = (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity > 0.) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity;
				vol_ch += (Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity > 0.) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity;
				double vol_dc = -(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity < 0.) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.scheduled_capacity;
				vol_dc += -(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity < 0.) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.scheduled_capacity;
				vol_ch /= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency;
				vol_dc *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.BESS.efficiency;
				vol_ch *= scale;
				vol_dc *= scale;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement.volume_demand.BESS += vol_ch;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement.volume_supply.BESS += vol_dc;

				// Balancing settlement
				agent_balancing_settlement_calculation(tick, node_ID, Power_market_inform, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement, 1);

				// Market operation update
				market_operation_update(tick, bz_ID, Power_market_inform.International_Market.operation.end_user, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results);
			}
		}
		//std::cout << "\n";
	}

	void industrial_status_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, bool control_reserve_flag){
		int price_interval = power_market::parameters::price_interval();

		int industrial_HV_num = Power_market_inform.agent_profiles.industrial.HV.size();
		for(int agent_iter = 0; agent_iter < industrial_HV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.industrial.HV[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);

			double marginal_price = Power_market_inform.TSO_Market.actual.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			// Calculate actual results
			agent_actual_results_calculation(node_ID, marginal_price_ID, Power_market_inform, Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids, Power_market_inform.agent_profiles.industrial.HV[agent_iter].results, control_reserve_flag);

			// Balancing settlement
			agent_balancing_settlement_calculation(tick, node_ID, Power_market_inform, Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids, Power_market_inform.agent_profiles.industrial.HV[agent_iter].results, Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement);

			// Market operation update
			market_operation_update(tick, bz_ID, Power_market_inform.International_Market.operation.industrial, Power_market_inform.agent_profiles.industrial.HV[agent_iter].results);
		}
		//std::cout << "\n";
	}

	void power_supplier_status_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, bool control_reserve_flag){
		int price_interval = power_market::parameters::price_interval();

		int hydro_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant.size();
		for(int agent_iter = 0; agent_iter < hydro_HV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);

			double marginal_price = Power_market_inform.TSO_Market.actual.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			// Calculate actual results
			agent_actual_results_calculation(node_ID, marginal_price_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].results, control_reserve_flag);

			// Balancing settlement
			agent_balancing_settlement_calculation(tick, node_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].settlement);

			// Market operation update
			market_operation_update(tick, bz_ID, Power_market_inform.International_Market.operation.hydro, Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].results);
		}

		int hydro_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant.size();
		for(int agent_iter = 0; agent_iter < hydro_LV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);

			double marginal_price = Power_market_inform.TSO_Market.actual.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			// Calculate actual results
			agent_actual_results_calculation(node_ID, marginal_price_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].results, control_reserve_flag);

			// Balancing settlement
			agent_balancing_settlement_calculation(tick, node_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].settlement, 1);

			// Market operation update
			market_operation_update(tick, bz_ID, Power_market_inform.International_Market.operation.hydro, Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].results);
		}

		int wind_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant.size();
		for(int agent_iter = 0; agent_iter < wind_HV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);

			double marginal_price = Power_market_inform.TSO_Market.actual.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			// Calculate actual results
			agent_actual_results_calculation(node_ID, marginal_price_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].results, control_reserve_flag);

			// Balancing settlement
			agent_balancing_settlement_calculation(tick, node_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].settlement, 1);

			// Market operation update
			market_operation_update(tick, bz_ID, Power_market_inform.International_Market.operation.wind, Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].results);
		}

		int wind_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant.size();
		for(int agent_iter = 0; agent_iter < wind_LV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);

			double marginal_price = Power_market_inform.TSO_Market.actual.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			// Calculate actual results
			agent_actual_results_calculation(node_ID, marginal_price_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results, control_reserve_flag);

			// Balancing settlement
			agent_balancing_settlement_calculation(tick, node_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].settlement, 1);

			// Market operation update
			market_operation_update(tick, bz_ID, Power_market_inform.International_Market.operation.wind, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results);
		}

		int pump_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV.size();
		for(int agent_iter = 0; agent_iter < pump_HV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);

			double marginal_price = Power_market_inform.TSO_Market.actual.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			// Calculate actual results
			agent_actual_results_calculation(node_ID, marginal_price_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].results, control_reserve_flag);

			// Balancing settlement
			agent_balancing_settlement_calculation(tick, node_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].settlement, 1);

			// Market operation update
			market_operation_update(tick, bz_ID, Power_market_inform.International_Market.operation.pump_storage, Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].results);
		}

		int pump_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
		for(int agent_iter = 0; agent_iter < pump_LV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);

			double marginal_price = Power_market_inform.TSO_Market.actual.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			// Calculate actual results
			agent_actual_results_calculation(node_ID, marginal_price_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].results, control_reserve_flag);

			// Balancing settlement
			agent_balancing_settlement_calculation(tick, node_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].settlement, 1);

			// Market operation update
			market_operation_update(tick, bz_ID, Power_market_inform.International_Market.operation.pump_storage, Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].results);
		}

		int slack_LV_num = Power_market_inform.agent_profiles.power_supplier.slack.LV_plant.size();
		for(int agent_iter = 0; agent_iter < slack_LV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);

			double marginal_price = Power_market_inform.TSO_Market.actual.price(tick, node_ID);
			int marginal_price_ID = Power_market_inform.price_map.price_ID[marginal_price];

			// Calculate actual results
			agent_actual_results_calculation(node_ID, marginal_price_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].results, control_reserve_flag);

			// Balancing settlement
			agent_balancing_settlement_calculation(tick, node_ID, Power_market_inform, Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].bids, Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].results, Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].settlement, 1);

			// Market operation update
			market_operation_update(tick, bz_ID, Power_market_inform.International_Market.operation.slack, Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].results);
//			if(Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].results.actual_supply > 0.){
//				std::cout << Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].results.actual_supply << "\t";
//				std::cout << Power_market_inform.International_Market.operation.slack.redispatch(tick, bz_ID) << "\n";
//			}
		}
	}

	void agents_balancing_settlement(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int node_num = Power_market_inform.TSO_Market.network.num_vertice;
		int bz_num = Power_market_inform.International_Market.cross_border_zone_start;

		// Balancing price per energy at each transmission node
		Eigen::VectorXd balancing_price_total = Eigen::VectorXd::Zero(bz_num);
		Eigen::VectorXd imbalance_total = Eigen::VectorXd::Zero(bz_num);
		Eigen::VectorXd imbalance_supply = Eigen::VectorXd::Zero(bz_num);
		Eigen::VectorXd balancing_price_supply = Eigen::VectorXd::Zero(bz_num);
		for(int node_iter = 0; node_iter < node_num; ++ node_iter){
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_iter);

			balancing_price_total(bz_ID) += Power_market_inform.TSO_Market.balancing.price_down(tick, node_iter);
			balancing_price_total(bz_ID) += Power_market_inform.TSO_Market.balancing.price_up(tick, node_iter);

			imbalance_total(bz_ID) += Power_market_inform.TSO_Market.imbalance.supply_up(tick, node_iter);
			imbalance_supply(bz_ID) += Power_market_inform.TSO_Market.imbalance.supply_up(tick, node_iter);
			imbalance_total(bz_ID) += Power_market_inform.TSO_Market.imbalance.demand_down(tick, node_iter);
			imbalance_total(bz_ID) -= Power_market_inform.TSO_Market.imbalance.supply_down(tick, node_iter);
			imbalance_supply(bz_ID) -= Power_market_inform.TSO_Market.imbalance.supply_down(tick, node_iter);
			imbalance_total(bz_ID) -= Power_market_inform.TSO_Market.imbalance.demand_up(tick, node_iter);
		}
		balancing_price_total = balancing_price_total.array() / imbalance_total.array();
		balancing_price_supply = balancing_price_total.array() * imbalance_supply.array();

		// Supply side balancing cost
		Eigen::VectorXd balancing_actual_demand = Eigen::VectorXd::Zero(bz_num);
		Eigen::VectorXd balancing_actual_supply = Eigen::VectorXd::Zero(bz_num);
		Eigen::VectorXd balancing_cost_supply_self = Eigen::VectorXd::Zero(bz_num);
		Eigen::VectorXd balancing_cost_supply_export = Eigen::VectorXd::Zero(bz_num);
		Eigen::VectorXd balancing_quan_supply_import = Eigen::VectorXd::Zero(bz_num);
		Eigen::VectorXd balancing_additional_price = Eigen::VectorXd::Zero(bz_num);

		for(int node_iter = 0; node_iter < node_num; ++ node_iter){
			int bz_ID = Power_network_inform.nodes.bidding_zone(node_iter);

			balancing_actual_demand(bz_ID) += Power_market_inform.TSO_Market.actual.demand(tick, node_iter);
			balancing_actual_supply(bz_ID) += Power_market_inform.TSO_Market.actual.supply(tick, node_iter);
		}

		for(int zone_iter = 0; zone_iter < bz_num; ++ zone_iter){
			balancing_cost_supply_self(zone_iter) = balancing_price_supply(zone_iter) * std::min(balancing_actual_demand(zone_iter), balancing_actual_supply(zone_iter)) / balancing_actual_supply(zone_iter);
			balancing_quan_supply_import(zone_iter) = balancing_actual_demand(zone_iter) - std::min(balancing_actual_demand(zone_iter), balancing_actual_supply(zone_iter));
			balancing_cost_supply_export(zone_iter) = balancing_price_supply(zone_iter) - balancing_cost_supply_self(zone_iter);
		}

		for(int zone_iter = 0; zone_iter < bz_num; ++ zone_iter){
			balancing_additional_price(zone_iter) += balancing_cost_supply_self(zone_iter);
			balancing_additional_price(zone_iter) += balancing_quan_supply_import(zone_iter) / balancing_quan_supply_import.sum() * balancing_cost_supply_export.sum();
			balancing_additional_price(zone_iter) /= balancing_actual_demand(zone_iter);
		}

		// End-users
		int point_num = Power_network_inform.points.bidding_zone.size();
		//int sample_num = agent::end_user::parameters::sample_num();
		int sample_num = Power_market_inform.agent_profiles.end_user_type.sample_num;
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int bz_ID = Power_network_inform.points.bidding_zone(point_iter);

			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				double balancing_price = 0.;
				double balancing_quan = 0.;

				balancing_quan += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement.volume_supply_up.imbalance;
				balancing_quan += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement.volume_demand_down.imbalance;
				balancing_quan += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement.volume_supply_down.imbalance;
				balancing_quan += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement.volume_demand_up.imbalance;
				balancing_price += balancing_quan * balancing_price_total(bz_ID);
				balancing_price += balancing_additional_price(bz_ID) * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_demand;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement.price.balancing += balancing_price;
			}
		}

		// Industrial demand
		int industrial_HV_num = Power_market_inform.agent_profiles.industrial.HV.size();
		for(int agent_iter = 0; agent_iter < industrial_HV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.industrial.HV[agent_iter].point_ID;
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);

			double balancing_price = 0.;
			double balancing_quan = 0.;

			balancing_quan += Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement.volume_supply_up.imbalance;
			balancing_quan += Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement.volume_demand_down.imbalance;
			balancing_quan += Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement.volume_supply_down.imbalance;
			balancing_quan += Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement.volume_demand_up.imbalance;
			balancing_price += balancing_quan * balancing_price_total(bz_ID);
			balancing_price += balancing_additional_price(bz_ID) * Power_market_inform.agent_profiles.industrial.HV[agent_iter].results.actual_demand;
			Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement.price.balancing += balancing_price;
		}

		// Cross-border flow
		int edge_num = Power_market_inform.agent_profiles.cross_border.size();
		for(int edge_iter = 0; edge_iter < edge_num; ++ edge_iter){
			int node_num = Power_market_inform.agent_profiles.cross_border[edge_iter].node_num;
			int entry_bz_ID = Power_market_inform.agent_profiles.cross_border[edge_iter].entry_bz_ID;

			// Store cross-border transmission flow as inflexible supply / demand at entry nodes
			if(node_num == 0){
				continue;
			}

			for(int node_iter = 0; node_iter < node_num; ++ node_iter){
				int node_ID = Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].node_ID;
				int bz_ID = Power_network_inform.nodes.bidding_zone(node_ID);

				double balancing_price;
				double balancing_quan = 0.;

				balancing_quan += Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].settlement.volume_supply_up.imbalance;
				balancing_quan +=Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].settlement.volume_demand_down.imbalance;
				balancing_quan += Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].settlement.volume_supply_down.imbalance;
				balancing_quan += Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].settlement.volume_demand_up.imbalance;
				balancing_price += balancing_quan * balancing_price_total(bz_ID);
				balancing_price += balancing_additional_price(bz_ID) * Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results.actual_demand;
				Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].settlement.price.balancing += balancing_price;
			}
		}
	}

	void aggregator_price_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int foresight_time = agent::aggregator::parameters::foresight_time();
		int aggregator_num = Power_market_inform.agent_profiles.aggregators.size();

		for(int agent_iter = 0; agent_iter < aggregator_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.aggregators[agent_iter].point_ID;
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
//			Eigen::VectorXd bid_vec = Power_market_inform.International_Market.merit_order_curve.col(bz_ID);
//			std::cout << bid_vec << "\n\n";
//			std::cout << Power_network_inform.plants.hydro.cap(agent_iter) << "\n\n";
//			bid_vec *= Power_network_inform.plants.hydro.cap(agent_iter);
//			bid_vec /= (Power_market_inform.International_Market.merit_order_curve.col(bz_ID).sum());

			Power_market_inform.agent_profiles.aggregators[agent_iter].price_expected_profile = Power_market_inform.International_Market.confirmed.price.col(bz_ID).segment(tick, foresight_time);
			Power_market_inform.agent_profiles.aggregators[agent_iter].price_demand_profile = Power_market_inform.International_Market.confirmed.price.col(bz_ID).segment(tick, foresight_time).array() + Power_market_inform.agent_profiles.aggregators[agent_iter].arbitrage_demand;
			Power_market_inform.agent_profiles.aggregators[agent_iter].price_supply_profile = Power_market_inform.International_Market.confirmed.price.col(bz_ID).segment(tick, foresight_time).array() - Power_market_inform.agent_profiles.aggregators[agent_iter].arbitrage_supply;

			// Update bids and results information
			agent_bids_initialization(Power_market_inform.agent_profiles.aggregators[agent_iter].bids);
			agent_results_set(Power_market_inform.agent_profiles.aggregators[agent_iter].results);
		}
	}

	void end_user_submit_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, configuration::process_config &process_par){
		int foresight_time = agent::end_user::parameters::foresight_time();
		int point_num = Power_network_inform.points.bidding_zone.size();
		//int sample_num = agent::end_user::parameters::sample_num();
		int sample_num = Power_market_inform.agent_profiles.end_user_type.sample_num;
		int price_interval = power_market::parameters::price_interval();
		double residential_ratio = agent::parameters::residential_ratio();

		// Update of forecast demand profile and operation strategies
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
              int node_ID = Power_network_inform.points.node(point_iter);

			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				int load_shift_time = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.shift_time;

				// Update of input profiles
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.house_default_period = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.house_schedule(tick);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.default_demand_profile = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.demand_profile(tick);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.default_demand_profile *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].investment.decision.EV_self_charging;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.default_demand_profile = Power_network_inform.points.nominal_mean_demand_field.row(point_iter).segment(tick, foresight_time);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.default_demand_profile *= agent::parameters::residential_ratio();
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand.head(foresight_time + load_shift_time - 1) = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand.tail(foresight_time + load_shift_time - 1);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand(foresight_time + load_shift_time - 1) = Power_network_inform.points.nominal_mean_demand_field(point_iter, tick + foresight_time - 1) * agent::parameters::residential_ratio();
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.unfulfilled_demand(foresight_time + load_shift_time - 1) *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].investment.decision.smart_appliance * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.scale;
				if(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.EV.BESS.energy_scale != 0.){
					Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.default_demand_profile *= 1.;
				}
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.default_demand_profile *= 1. - Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].investment.decision.smart_appliance * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.scale;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.default_PV_profile = Power_network_inform.points.solar_cf.row(point_iter).segment(tick, foresight_time);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.default_PV_profile *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.PV_scale;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.PV_output = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.default_PV_profile(0);
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.price_demand_profile = Power_market_inform.agent_profiles.aggregators[point_iter].price_demand_profile;
				Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.price_supply_profile = Power_market_inform.agent_profiles.aggregators[point_iter].price_supply_profile;

				// Store demand characteristic of each node
				double demand_inflex;
				double demand_shiftable;
                demand_inflex = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.default_demand_profile(0);
                demand_inflex *= Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
                demand_inflex *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.weight;
                Power_market_inform.TSO_Market.flex_stat_end.demand_inflex(tick, node_ID) += demand_inflex;
                demand_shiftable = demand_inflex;
                demand_shiftable *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].investment.decision.smart_appliance * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.scale;
                demand_shiftable /= 1. -Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].investment.decision.smart_appliance * Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.smart_appliance.scale;
                demand_shiftable *= Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.weight;
                Power_market_inform.TSO_Market.flex_stat_end.demand_shiftable(tick, node_ID) += demand_shiftable;
//				if(Power_market_inform.agent_profiles.end_user_type.contingency[sample_iter] != 0){
//
//				}
//				else{
//                    demand_inflex = 0.;
//				}

				// Set bids and results information
				agent_bids_initialization(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids);
				agent_results_set(Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results);

				// Totally inflexible end-user, demand profile as default
				double scale = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.weight;
				scale *= Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
				Power_market_inform.TSO_Market.flex_stat_end.demand_inflex(tick, node_ID) -= demand_inflex;
				if(!Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].investment.decision.smart_management){
					agent::end_user::end_user_no_LP(tick, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter], demand_inflex, scale, process_par);
				}
				// Flexible end-user
				else{
                    // just for test (comment later)
					//agent::end_user::end_user_no_LP(tick, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter], demand_inflex, scale, process_par);
                    // just for test (comment later)

					// Optimization and update process variables
					agent::end_user::end_user_LP_optimize(tick, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter], process_par);
				}
				// Update inflexible demand
//				if(Power_market_inform.agent_profiles.end_user_type.contingency[sample_iter] == 0){
//                    demand_inflex = 0.;
//				}
				Power_market_inform.TSO_Market.flex_stat_end.demand_inflex(tick, node_ID) += demand_inflex;

				// Scale the bids correctly
				agent_submitted_bids_scale(scale, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids);
			}
		}
	}

	void industrial_submit_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
		int industrial_HV_num = Power_market_inform.agent_profiles.industrial.HV.size();
		int price_interval = power_market::parameters::price_interval();

		for(int agent_iter = 0; agent_iter < industrial_HV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.industrial.HV[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);

			double bid_inflex_industrial = Power_network_inform.points.nominal_mean_demand_field(point_ID, tick);
			bid_inflex_industrial *= Power_network_inform.points.population_density(point_ID) * Power_network_inform.points.point_area / 1000.;
			bid_inflex_industrial *= 1. - agent::parameters::residential_ratio();
			double bid_flex_industrial = bid_inflex_industrial;
			bid_inflex_industrial *= 1. - agent::industrial::parameters::flexible_ratio();
			bid_flex_industrial *= agent::industrial::parameters::flexible_ratio();
            Power_market_inform.TSO_Market.flex_stat_end.demand_inflex(tick, node_ID) += bid_inflex_industrial;
            Power_market_inform.TSO_Market.flex_stat_end.demand_flex(tick, node_ID) += bid_flex_industrial;
//			if(Power_market_inform.agent_profiles.end_user_type.contingency[0]){
//
//			}
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
			int node_ID = Power_network_inform.points.node(point_ID);
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
			Power_market_inform.TSO_Market.flex_stat.supply_flex(tick, node_ID) += bid_vec.sum();
		}

		int hydro_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant.size();
		for(int agent_iter = 0; agent_iter < hydro_LV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
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
			Power_market_inform.TSO_Market.flex_stat.supply_flex(tick, node_ID) += bid_vec.sum();
		}

		int wind_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant.size();
		for(int agent_iter = 0; agent_iter < wind_HV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int price_supply_flex_ID = 0;
			double bid_quan = Power_network_inform.points.wind_on_cf(point_ID, tick) * Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].cap;

			// Set bids information
			agent_bids_initialization(Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids);
			agent_results_set(Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].results);
			Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids.submitted_supply_flex(price_supply_flex_ID) = bid_quan;
			Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids.redispatch_supply(price_supply_flex_ID) = bid_quan;
			Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids.filter_supply(price_supply_flex_ID) = bid_quan;
			Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids.balancing_supply(price_supply_flex_ID) = bid_quan;
			Power_market_inform.TSO_Market.flex_stat.supply_flex(tick, node_ID) += bid_quan;
		}

		int wind_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant.size();
		for(int agent_iter = 0; agent_iter < wind_LV_plant_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
			int price_supply_flex_ID = 0;
			double bid_quan = Power_network_inform.points.wind_on_cf(point_ID, tick) * Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].cap;

			// Set bids information
			agent_bids_initialization(Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids);
			agent_results_set(Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results);
			Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.submitted_supply_flex(price_supply_flex_ID) = bid_quan;
			Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.redispatch_supply(price_supply_flex_ID) = bid_quan;
			Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.filter_supply(price_supply_flex_ID) = bid_quan;
			Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.balancing_supply(price_supply_flex_ID) = bid_quan;
			Power_market_inform.TSO_Market.flex_stat.supply_flex(tick, node_ID) += bid_quan;
		}

		int pump_HV_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV.size();
		for(int agent_iter = 0; agent_iter < pump_HV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
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
			Power_market_inform.TSO_Market.flex_stat.supply_flex(tick, node_ID) += bid_vec.sum();
		}

		int pump_LV_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
		for(int agent_iter = 0; agent_iter < pump_LV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);
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
			Power_market_inform.TSO_Market.flex_stat.supply_flex(tick, node_ID) += bid_vec.sum();
		}

		int slack_LV_num = Power_market_inform.agent_profiles.power_supplier.slack.LV_plant.size();
		for(int agent_iter = 0; agent_iter < slack_LV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].point_ID;
			int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
//			Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].fix_cost = 30. + Power_market_inform.International_Market.confirmed.price(tick, bz_ID);
//			int price_supply_flex_ID = Power_market_inform.price_map.price_ID[Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].fix_cost];
			Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].var_cost = Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].fix_cost;
			Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].var_cost += Power_market_inform.International_Market.confirmed.price(tick, bz_ID);
			int price_supply_flex_ID = Power_market_inform.price_map.price_ID[Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].var_cost];
			int price_supply_max_ID = Power_market_inform.price_map.price_ID[99.5];
			int price_length = std::max(price_supply_max_ID - price_supply_flex_ID + 1, 1);

			double bid_quan = Power_network_inform.points.nominal_mean_demand_field(point_ID, tick);
			bid_quan *= Power_network_inform.points.population_density(point_ID);
			bid_quan *= Power_network_inform.points.point_area / 1000.;
			//bid_quan *= 1.;
//			double bid_quan = std::numeric_limits<double>::infinity();
			Eigen::VectorXd bid_vec = Eigen::VectorXd::Zero(price_interval + 2);
			bid_vec.segment(price_supply_flex_ID, price_length) = Eigen::VectorXd::Ones(price_length) * bid_quan / price_length;

			// Set bids information
			agent_bids_initialization(Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].bids);
			agent_results_set(Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].results);
			Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].bids.submitted_supply_flex = bid_vec;
			Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].bids.redispatch_supply = bid_vec;
			Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].bids.filter_supply = bid_vec;
			Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].bids.balancing_supply = bid_vec;
		}
	}
}

void agent::agents_set(int start_time, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, std::string fin_name, configuration::process_config &process_par){
	auto fin_dim = basic::get_file_dim(fin_name, 1);
	auto end_user_type = basic::read_config_file(fin_name);
	Power_market_inform.agent_profiles.end_user_type.initialize(fin_dim[1]);
	for(int sample_iter = 0; sample_iter < fin_dim[1]; ++ sample_iter){
        Power_market_inform.agent_profiles.end_user_type.weight[sample_iter] = stod(end_user_type["ratio"][sample_iter]);
        //Power_market_inform.agent_profiles.end_user_type.dynamic_tariff[sample_iter] = (bool) stod(end_user_type["dynamic_tariff"][sample_iter]);
        Power_market_inform.agent_profiles.end_user_type.smart_management[sample_iter] = (bool) stod(end_user_type["smart_management"][sample_iter]);
        Power_market_inform.agent_profiles.end_user_type.smart_appliance[sample_iter] = (bool) stod(end_user_type["smart_appliance"][sample_iter]);
        Power_market_inform.agent_profiles.end_user_type.PV_scale[sample_iter] = stod(end_user_type["PV_scale"][sample_iter]);
        Power_market_inform.agent_profiles.end_user_type.BESS_energy[sample_iter] = stod(end_user_type["BESS_energy"][sample_iter]);
        Power_market_inform.agent_profiles.end_user_type.BESS_capacity[sample_iter] = stod(end_user_type["BESS_capacity"][sample_iter]);
        Power_market_inform.agent_profiles.end_user_type.EV_energy[sample_iter] = stod(end_user_type["EV_energy"][sample_iter]);
        Power_market_inform.agent_profiles.end_user_type.EV_capacity[sample_iter] = stod(end_user_type["EV_capacity"][sample_iter]);
        Power_market_inform.agent_profiles.end_user_type.redispatch[sample_iter] = (bool) stod(end_user_type["redispatch"][sample_iter]);
        Power_market_inform.agent_profiles.end_user_type.control_reserve[sample_iter] = (bool) stod(end_user_type["control_reserve"][sample_iter]);
        Power_market_inform.agent_profiles.end_user_type.contingency[sample_iter] = (bool) stod(end_user_type["contingency"][sample_iter]);
	}

	Power_market_inform.agent_profiles.aggregators = aggregator_set(start_time, Power_market_inform.International_Market, Power_network_inform);
	Power_market_inform.agent_profiles.cross_border = cross_border_set(Power_market_inform, Power_network_inform);
	Power_market_inform.agent_profiles.end_users = end_user_set(start_time, Power_market_inform, Power_network_inform, process_par);
	Power_market_inform.agent_profiles.industrial = industrial_set(start_time, Power_market_inform, Power_network_inform);
	Power_market_inform.agent_profiles.power_supplier = power_supplier_set(start_time, Power_market_inform, Power_network_inform);
}

void agent::agents_redispatch_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, configuration::process_config &process_par){
	cross_border_redispatch_update(tick, Power_market_inform, Power_network_inform);
	end_user_redispatch_update(tick, Power_market_inform, Power_network_inform, process_par);
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
	cross_border_balancing_update(tick, Power_market_inform, Power_network_inform);
	end_user_balancing_update(tick, Power_market_inform, Power_network_inform);
	industrial_balancing_update(tick, Power_market_inform, Power_network_inform);
	power_supplier_balancing_update(tick, Power_market_inform, Power_network_inform);
	agents_redispatch_settlement(tick, Power_market_inform, Power_network_inform);
}

void agent::agents_status_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, configuration::process_config &process_par){
	cross_border_status_update(tick, Power_market_inform, Power_network_inform, process_par.control_reserve_flag);
	end_user_status_update(tick, Power_market_inform, Power_network_inform, process_par);
	industrial_status_update(tick, Power_market_inform, Power_network_inform, process_par.control_reserve_flag);
	power_supplier_status_update(tick, Power_market_inform, Power_network_inform, process_par.control_reserve_flag);
	agents_balancing_settlement(tick, Power_market_inform, Power_network_inform);
}

void agent::agents_submit_update(int tick, power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, configuration::process_config &process_par){
	aggregator_price_update(tick, Power_market_inform, Power_network_inform);
	end_user_submit_update(tick, Power_market_inform, Power_network_inform, process_par);
	industrial_submit_update(tick, Power_market_inform, Power_network_inform);
	power_supplier_submit_update(tick, Power_market_inform, Power_network_inform);
}
