// Source file for re-dispatch and tertiary control reserve market clearing of TSO in Norway
#include <iostream>
//#include <chrono>
//#include "../basic/LP_gpa.h"
#include "src/basic/rw_csv.h"
#include "power_market.h"

void power_market::TSO_Market_Set(market_inform &TSO_Market, power_network::network_inform &Power_network_inform, int Time){
	// Input parameters of TSO market
	TSO_Market.num_zone = Power_network_inform.nodes.bidding_zone.size();
	TSO_Market.time_intervals = Time;
	parameters::bidded_price(TSO_Market.bidded_price_map);

	// Set node admittance matrix and line capacity matrix
	TSO_Market.network.num_vertice = TSO_Market.num_zone;
	Eigen::MatrixXd admittance = Eigen::MatrixXd::Zero(TSO_Market.network.num_vertice, TSO_Market.network.num_vertice);
	Eigen::MatrixXd capacity = Eigen::MatrixXd::Zero(TSO_Market.network.num_vertice, TSO_Market.network.num_vertice);
	for(int edge_iter = 0; edge_iter < Power_network_inform.edges.distance.size(); ++ edge_iter){
		int from_ID = Power_network_inform.edges.from(edge_iter );
		int to_ID = Power_network_inform.edges.to(edge_iter);
		int voltage = Power_network_inform.edges.voltage_base(edge_iter);
		double y_ij = 1. / Power_network_inform.edges.distance(edge_iter);
		y_ij /= Power_network_inform.tech_parameters.z_trans_series.imag();
		y_ij *= Power_network_inform.tech_parameters.impedenace_base_levels[voltage];
		admittance(from_ID, to_ID) += y_ij ;
		admittance(to_ID, from_ID) += y_ij ;
		capacity(from_ID, to_ID) += Power_network_inform.tech_parameters.power_limit[voltage];
	}

	// Set compact incidence matrix and edge admittance matrix
	double tol = 1E-6;
	TSO_Market.network.incidence.reserve(TSO_Market.network.num_vertice * TSO_Market.network.num_vertice);
	TSO_Market.network.admittance.reserve(TSO_Market.network.num_vertice * TSO_Market.network.num_vertice);
	std::vector <double> power_limit;
	power_limit.reserve(TSO_Market.network.num_vertice * TSO_Market.network.num_vertice);
	for(int row_iter = 0; row_iter < TSO_Market.network.num_vertice - 1; ++ row_iter){
		for(int col_iter = row_iter + 1; col_iter < TSO_Market.network.num_vertice; ++ col_iter){
			if(abs(admittance(row_iter , col_iter)) > tol){
				TSO_Market.network.incidence.push_back(Eigen::Vector2i(row_iter, col_iter));
				TSO_Market.network.admittance.push_back(admittance(row_iter , col_iter));
				power_limit.push_back(capacity(row_iter , col_iter));
				//std::cout << TSO_Market.network.incidence.size() - 1 << ":\t" << admittance(row_iter , col_iter) << "\t" << capacity(row_iter , col_iter) << "\n";
			}
		}
	}
	TSO_Market.network.num_edges = TSO_Market.network.incidence.size();

	// Set voltage and power constraints at each edge
	TSO_Market.network.voltage_constraint = Eigen::MatrixXd::Ones(TSO_Market.network.num_vertice, 2);
	TSO_Market.network.voltage_constraint.col(0) *= -Power_network_inform.tech_parameters.theta_limit;
	TSO_Market.network.voltage_constraint.col(1) *= Power_network_inform.tech_parameters.theta_limit;
	TSO_Market.network.power_constraint = Eigen::MatrixXd (TSO_Market.network.num_edges, 2);
	TSO_Market.network.power_constraint.col(1) = Eigen::Map <Eigen::VectorXd> (power_limit.data(), power_limit.size());
	TSO_Market.network.power_constraint.col(1) /= Power_network_inform.tech_parameters.s_base;
	TSO_Market.network.power_constraint.col(0) = -TSO_Market.network.power_constraint.col(1);

	// Initialization of process variables
	power_market::Market_Initialization(TSO_Market);

	// Initialization of output variables
	TSO_Market.confirmed.supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed.demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.confirmed.price = Eigen::MatrixXd(Time, TSO_Market.num_zone);
	TSO_Market.EOM.cost = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.EOM.utility = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.redispatch.cost = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.redispatch.utility = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.redispatch.supply_up = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.redispatch.supply_down = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.redispatch.demand_up = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.redispatch.demand_down = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.redispatch.price_supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.redispatch.price_demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.balancing.supply_up = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.balancing.supply_down = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.balancing.demand_up = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.balancing.demand_down = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.balancing.price_supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.balancing.price_demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.actual.supply = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.actual.demand = Eigen::MatrixXd::Zero(Time, TSO_Market.num_zone);
	TSO_Market.actual.price = Eigen::MatrixXd(Time, TSO_Market.num_zone);
}

void power_market::Confirmed_bid_calculation(int tick, market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
	// Initialize submit bids of the TSO market
	Market_Initialization(Power_market_inform.TSO_Market);

	// Set reference price for redispatch settlement
	for(int node_iter = 0; node_iter < Power_market_inform.TSO_Market.network.num_vertice; ++ node_iter){
		int bz_ID = Power_network_inform.nodes.bidding_zone(node_iter);
		Power_market_inform.TSO_Market.reference_price(node_iter) = Power_market_inform.International_Market.confirmed.price(tick, bz_ID);
	}

	// Initialize boundary conditions with other bidding zones
	TSO_boundary_update(tick, Power_market_inform.TSO_Market, Power_market_inform.International_Market, Power_network_inform);

	// Residential demand
	int point_num = Power_network_inform.points.bidding_zone.size();
	int sample_num = Power_market_inform.agent_profiles.end_users[0].size();
	for(int point_iter = 0; point_iter < point_num; ++ point_iter){
		int node_ID = Power_network_inform.points.node(point_iter);

		for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
			Power_market_inform.TSO_Market.submitted_demand.col(node_ID) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_demand;
			Power_market_inform.TSO_Market.submitted_supply.col(node_ID) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.redispatch_supply;
		}
	}

	// Industrial demand
	int industrial_HV_num = Power_market_inform.agent_profiles.industrial.HV.size();
	for(int agent_iter = 0; agent_iter < industrial_HV_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.industrial.HV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		Power_market_inform.TSO_Market.submitted_demand.col(node_ID) += Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids.redispatch_demand;
	}

	// Hydroelectric power plants
	int hydro_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		Power_market_inform.TSO_Market.submitted_supply.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].bids.redispatch_supply;
		Power_market_inform.TSO_Market.submitted_demand.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].bids.redispatch_demand;
	}
	int hydro_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		Power_market_inform.TSO_Market.submitted_supply.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids.redispatch_supply;
		Power_market_inform.TSO_Market.submitted_demand.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids.redispatch_demand;
	}

	// Wind power plants
	int wind_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		Power_market_inform.TSO_Market.submitted_supply.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids.redispatch_supply;
		Power_market_inform.TSO_Market.submitted_demand.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids.redispatch_demand;
	}
	int wind_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		Power_market_inform.TSO_Market.submitted_supply.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.redispatch_supply;
		Power_market_inform.TSO_Market.submitted_demand.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.redispatch_demand;
	}

	// Pump storage
	int pump_HV_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV.size();
	for(int agent_iter = 0; agent_iter < pump_HV_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);
		Power_market_inform.TSO_Market.submitted_supply.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].bids.redispatch_supply;
		Power_market_inform.TSO_Market.submitted_demand.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].bids.redispatch_demand;
	}
	int pump_LV_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
	for(int agent_iter = 0; agent_iter < pump_LV_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		Power_market_inform.TSO_Market.submitted_supply.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_supply;
		Power_market_inform.TSO_Market.submitted_demand.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.redispatch_demand;
	}
}

void power_market::TSO_Market_Scheduled_Results_Get(int tick, market_inform &Market, alglib::minlpstate &Problem){
	alglib::real_1d_array sol;
	alglib::minlpreport rep;
	alglib::minlpresults(Problem, sol, rep);
	Eigen::VectorXd sol_vec = Eigen::Map <Eigen::VectorXd> (sol.getcontent(), sol.length());

	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		// Store power source / sink
		int row_start = 2 * Market.network.num_vertice + node_iter * (Market.price_intervals + 2);
		Market.confirmed.supply(tick, node_iter) = (sol_vec.segment(row_start, Market.price_intervals + 2).array().max(0)).sum();
		Market.confirmed.demand(tick, node_iter) = -(sol_vec.segment(row_start, Market.price_intervals + 2).array().min(0)).sum();

		// Store nodal prices
		Market.confirmed.price(tick, node_iter) = Market.bidded_price_map.bidded_price(0) + rep.lagbc[row_start];
		Market.confirmed.price(tick, node_iter) = int(Market.confirmed.price(tick, node_iter)) + .5;
		if(Market.confirmed.price(tick, node_iter) < Market.bidded_price_map.bidded_price(1)){
			Market.confirmed.price(tick, node_iter) = Market.bidded_price_map.bidded_price(0);
		}
		else if(Market.confirmed.price(tick, node_iter) > Market.bidded_price_map.bidded_price(Market.price_intervals)){
			Market.confirmed.price(tick, node_iter) = Market.bidded_price_map.bidded_price(Market.price_intervals + 1);
		}

		// Store ratio at nodes
		for(int price_iter = 0; price_iter < Market.price_intervals + 2; ++ price_iter){
			if(Market.bidded_price_map.bidded_price(price_iter) >= Market.confirmed.price(tick, node_iter) || price_iter == Market.price_intervals + 1){
				if(sol[row_start + price_iter] >= 0.){
					Market.confirmed.ratio_demand(node_iter) = std::min(Market.submitted_demand(price_iter, node_iter), Market.submitted_supply(price_iter, node_iter) - sol[row_start + price_iter]);
					Market.confirmed.ratio_supply(node_iter) = Market.confirmed.ratio_demand(node_iter) + sol[row_start + price_iter];
					Market.confirmed.ratio_demand(node_iter) /= Market.submitted_demand(price_iter, node_iter) + 1E-12;
					Market.confirmed.ratio_supply(node_iter) /= Market.submitted_supply(price_iter, node_iter) + 1E-12;
				}
				else{
					Market.confirmed.ratio_supply(node_iter) = std::min(Market.submitted_supply(price_iter, node_iter), Market.submitted_demand(price_iter, node_iter) + sol[row_start + price_iter]);
					Market.confirmed.ratio_demand(node_iter) = Market.confirmed.ratio_supply(node_iter) - sol[row_start + price_iter];
					Market.confirmed.ratio_demand(node_iter) /= Market.submitted_demand(price_iter, node_iter) + 1E-12;
					Market.confirmed.ratio_supply(node_iter) /= Market.submitted_supply(price_iter, node_iter) + 1E-12;
				}
				break;
			}
		}
	}

	std::cout << Market.confirmed.supply.row(tick).sum() << "\t" << Market.confirmed.demand.row(tick).sum() << "\n";
	std::cout << sol_vec.segment(Market.network.num_vertice, Market.network.num_vertice).minCoeff() << " " << sol_vec.segment(Market.network.num_vertice, Market.network.num_vertice).maxCoeff() << " " << .5 * sol_vec.segment(Market.network.num_vertice, Market.network.num_vertice).array().abs().sum() << "\n";
	std::cout << sol_vec.head(Market.network.num_vertice).minCoeff() << " " << sol_vec.head(Market.network.num_vertice).maxCoeff()  << "\n";
	std::cout << sol_vec.tail(Market.network.num_edges).minCoeff() << " " << sol_vec.tail(Market.network.num_edges).maxCoeff() << "\n\n";
}

void power_market::Balancing_bid_calculation(int tick, market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform){
	// Initialize submit bids of the TSO market
	Market_Initialization(Power_market_inform.TSO_Market);

	// Initialize boundary conditions with other bidding zones
	TSO_boundary_update(tick, Power_market_inform.TSO_Market, Power_market_inform.International_Market, Power_network_inform);

	// Residential demand
	int point_num = Power_network_inform.points.bidding_zone.size();
	int sample_num = Power_market_inform.agent_profiles.end_users[0].size();
	for(int point_iter = 0; point_iter < point_num; ++ point_iter){
		int node_ID = Power_network_inform.points.node(point_iter);

		for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
			Power_market_inform.TSO_Market.submitted_demand.col(node_ID) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_demand;
			Power_market_inform.TSO_Market.submitted_supply.col(node_ID) += Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.bids.balancing_supply;
		}
	}

	// Industrial demand
	int industrial_HV_num = Power_market_inform.agent_profiles.industrial.HV.size();
	for(int agent_iter = 0; agent_iter < industrial_HV_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.industrial.HV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		Power_market_inform.TSO_Market.submitted_demand.col(node_ID) += Power_market_inform.agent_profiles.industrial.HV[agent_iter].bids.balancing_demand;
	}

	// Hydroelectric power plants
	int hydro_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		Power_market_inform.TSO_Market.submitted_supply.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].bids.balancing_supply;
		Power_market_inform.TSO_Market.submitted_demand.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].bids.balancing_demand;
	}
	int hydro_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		Power_market_inform.TSO_Market.submitted_supply.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids.balancing_supply;
		Power_market_inform.TSO_Market.submitted_demand.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].bids.balancing_demand;
	}

	// Wind power plants
	int wind_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		Power_market_inform.TSO_Market.submitted_supply.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids.balancing_supply;
		Power_market_inform.TSO_Market.submitted_demand.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].bids.balancing_demand;
	}
	int wind_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		Power_market_inform.TSO_Market.submitted_supply.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.balancing_supply;
		Power_market_inform.TSO_Market.submitted_demand.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].bids.balancing_demand;
	}

	// Pump storage
	int pump_HV_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV.size();
	for(int agent_iter = 0; agent_iter < pump_HV_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);
		Power_market_inform.TSO_Market.submitted_supply.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].bids.balancing_supply;
		Power_market_inform.TSO_Market.submitted_demand.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].bids.balancing_demand;
	}
	int pump_LV_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
	for(int agent_iter = 0; agent_iter < pump_LV_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		Power_market_inform.TSO_Market.submitted_supply.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.balancing_supply;
		Power_market_inform.TSO_Market.submitted_demand.col(node_ID) += Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].bids.balancing_demand;
	}
}

void power_market::TSO_Market_Actual_Results_Get(int tick, market_inform &Market, alglib::minlpstate &Problem){
	alglib::real_1d_array sol;
	alglib::minlpreport rep;
	alglib::minlpresults(Problem, sol, rep);
	Eigen::VectorXd sol_vec = Eigen::Map <Eigen::VectorXd> (sol.getcontent(), sol.length());

	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		// Store power source / sink
		int row_start = 2 * Market.network.num_vertice + node_iter * (Market.price_intervals + 2);
		Market.actual.supply(tick, node_iter) = (sol_vec.segment(row_start, Market.price_intervals + 2).array().max(0)).sum();
		Market.actual.demand(tick, node_iter) = -(sol_vec.segment(row_start, Market.price_intervals + 2).array().min(0)).sum();

		// Store nodal prices
		Market.actual.price(tick, node_iter) = Market.bidded_price_map.bidded_price(0) + rep.lagbc[row_start];
		Market.actual.price(tick, node_iter) = int(Market.actual.price(tick, node_iter)) + .5;
		if(Market.actual.price(tick, node_iter) < Market.bidded_price_map.bidded_price(1)){
			Market.actual.price(tick, node_iter) = Market.bidded_price_map.bidded_price(0);
		}
		else if(Market.actual.price(tick, node_iter) > Market.bidded_price_map.bidded_price(Market.price_intervals)){
			Market.actual.price(tick, node_iter) = Market.bidded_price_map.bidded_price(Market.price_intervals + 1);
		}

		// Store ratio at nodes
		for(int price_iter = 0; price_iter < Market.price_intervals + 2; ++ price_iter){
			if(Market.bidded_price_map.bidded_price(price_iter) >= Market.actual.price(tick, node_iter) || price_iter == Market.price_intervals + 1){
				if(sol[row_start + price_iter] >= 0.){
					Market.actual.ratio_demand(node_iter) = std::min(Market.submitted_demand(price_iter, node_iter), Market.submitted_supply(price_iter, node_iter) - sol[row_start + price_iter]);
					Market.actual.ratio_supply(node_iter) = Market.actual.ratio_demand(node_iter) + sol[row_start + price_iter];
					Market.actual.ratio_demand(node_iter) /= Market.submitted_demand(price_iter, node_iter) + 1E-12;
					Market.actual.ratio_supply(node_iter) /= Market.submitted_supply(price_iter, node_iter) + 1E-12;
				}
				else{
					Market.actual.ratio_supply(node_iter) = std::min(Market.submitted_supply(price_iter, node_iter), Market.submitted_demand(price_iter, node_iter) + sol[row_start + price_iter]);
					Market.actual.ratio_demand(node_iter) = Market.actual.ratio_supply(node_iter) - sol[row_start + price_iter];
					Market.actual.ratio_demand(node_iter) /= Market.submitted_demand(price_iter, node_iter) + 1E-12;
					Market.actual.ratio_supply(node_iter) /= Market.submitted_supply(price_iter, node_iter) + 1E-12;
				}
				break;
			}
		}
		//std::cout << node_iter << ":\t" << Market.actual_ratio_demand(node_iter) << "\t" << Market.actual_ratio_supply(node_iter) << "\n";
	}
	//std::cout << "\n";

	std::cout << Market.actual.supply.row(tick).sum() << "\t" << Market.actual.demand.row(tick).sum() << "\n";
	std::cout << sol_vec.segment(Market.network.num_vertice, Market.network.num_vertice).minCoeff() << " " << sol_vec.segment(Market.network.num_vertice, Market.network.num_vertice).maxCoeff() << " " << .5 * sol_vec.segment(Market.network.num_vertice, Market.network.num_vertice).array().abs().sum() << "\n";
	std::cout << sol_vec.head(Market.network.num_vertice).minCoeff() << " " << sol_vec.head(Market.network.num_vertice).maxCoeff()  << "\n";
	std::cout << sol_vec.tail(Market.network.num_edges).minCoeff() << " " << sol_vec.tail(Market.network.num_edges).maxCoeff() << "\n\n";
}
