// Source file for the complete procedure of the power market
#include "power_market.h"

void Market_Initialization(market_inform &Market){
	// Initialization of process variables
	// Should re-initialize for every time slice
	Market.submitted_supply = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Market.submitted_demand = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
}

void Market_clearing_nodal(int tick, market_inform &Market, Eigen::VectorXi &default_price_ID, Eigen::MatrixXd &bidded_supply, Eigen::MatrixXd &bidded_demand){
	Eigen::VectorXi price_demand_ID = (Market.price_intervals + 1) * Eigen::VectorXi::Ones(Market.num_zone);
	Eigen::VectorXi price_supply_ID = Eigen::VectorXi::Zero(Market.num_zone);
	
	double trade_quantity;
	for(int zone_ID = 0; zone_ID < Market.num_zone; ++ zone_ID){
		
		while(price_demand_ID(zone_ID) > price_supply_ID(zone_ID)){
			// Check if there are demand bids at current price interval
			while(bidded_demand(price_demand_ID(zone_ID), zone_ID) == 0){
				if(price_demand_ID(zone_ID) > 0){
					price_demand_ID(zone_ID) -= 1;
				}
				else{
					// No available buyer left to buy electricity
					default_price_ID(zone_ID) = price_supply_ID(zone_ID);
					break;
				}
			}
			
			// Check if there are supply bids at current price interval
			while(bidded_supply(price_supply_ID(zone_ID), zone_ID) == 0){
				if(price_supply_ID(zone_ID) < Market.bidded_price.size() - 1){
					price_supply_ID(zone_ID) += 1;
				}
				else{
					// No available seller left to sell electricity
					default_price_ID(zone_ID) = price_demand_ID(zone_ID);
					break;				
				}			
			}
			
			if(price_demand_ID(zone_ID) > price_supply_ID(zone_ID)){
				trade_quantity = std::min(bidded_supply(price_supply_ID(zone_ID), zone_ID), bidded_demand(price_demand_ID(zone_ID), zone_ID));
				Market.confirmed_supply(tick, zone_ID) += trade_quantity;
				Market.confirmed_demand(tick, zone_ID) += trade_quantity;
				bidded_supply(price_supply_ID(zone_ID), zone_ID) -= trade_quantity;
				bidded_demand(price_demand_ID(zone_ID), zone_ID) -= trade_quantity;
			}
			else{
				if(bidded_supply(price_supply_ID(zone_ID), zone_ID) >= 0){
					default_price_ID(zone_ID) = price_supply_ID(zone_ID);
				}
				else{
					default_price_ID(zone_ID) = price_demand_ID(zone_ID);
				}
			}
		}
		
		// Record default market clearing prices
		Market.confirmed_price(tick, zone_ID) = Market.bidded_price(default_price_ID(zone_ID));		
	}	
}