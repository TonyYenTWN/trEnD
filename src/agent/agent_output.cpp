// Source file for output settlement of different agents
#include <filesystem>
#include "agent_func.h"

namespace local{
	inline std::vector <std::string> var_names_set(){
		std::vector <std::string> var_names;
		var_names.push_back("ID");
		var_names.push_back("EOM_supply");
		var_names.push_back("EOM_demand");
		var_names.push_back("EOM_cost");
		var_names.push_back("EOM_utility");
		var_names.push_back("EOM_price");
		var_names.push_back("redispatch_supply_up");
		var_names.push_back("redispatch_supply_down");
		var_names.push_back("redispatch_demand_up");
		var_names.push_back("redispatch_demand_down");
		var_names.push_back("redispatch_cost");
		var_names.push_back("redispatch_utility");
		var_names.push_back("redispatch_price");
		var_names.push_back("redispatch_reimbursement");
		var_names.push_back("imbalance_supply_up");
		var_names.push_back("imbalance_supply_down");
		var_names.push_back("imbalance_demand_up");
		var_names.push_back("imbalance_demand_down");
		var_names.push_back("balancing_supply_up");
		var_names.push_back("balancing_supply_down");
		var_names.push_back("balancing_demand_up");
		var_names.push_back("balancing_demand_down");
		var_names.push_back("balancing_cost");
		var_names.push_back("balancing_utility");
		var_names.push_back("balancing_price");
		var_names.push_back("balancing_reimbursement");
		var_names.push_back("BESS_supply");
		var_names.push_back("BESS_demand");

		return var_names;
	}

	Eigen::RowVectorXd output_row_store(int col_num, agent::settlement_struct &settlement){
		Eigen::RowVectorXd output_row = Eigen::RowVectorXd::Zero(col_num);

		output_row(1) = settlement.volume_supply.EOM;
		output_row(2) = settlement.volume_demand.EOM;
		output_row(3) = settlement.cost_supply.EOM;
		output_row(4) = settlement.utility_demand.EOM;
		output_row(5) = settlement.price.EOM - settlement.utility_supply.EOM;
		output_row(6) = settlement.volume_supply_up.redispatch;
		output_row(7) = settlement.volume_supply_down.redispatch;
		output_row(8) = settlement.volume_demand_up.redispatch;
		output_row(9) = settlement.volume_demand_down.redispatch;
		output_row(10) = settlement.cost_supply.redispatch;
		output_row(11) = settlement.utility_demand.redispatch;
		output_row(12) = settlement.price.redispatch;
		output_row(13) = settlement.reimburse.redispatch + settlement.utility_supply.redispatch - settlement.cost_demand.redispatch;
		output_row(14) = settlement.volume_supply_up.imbalance;
		output_row(15) = settlement.volume_supply_down.imbalance;
		output_row(16) = settlement.volume_demand_up.imbalance;
		output_row(17) = settlement.volume_demand_down.imbalance;
		output_row(18) = settlement.volume_supply_up.balancing;
		output_row(19) = settlement.volume_supply_down.balancing;
		output_row(20) = settlement.volume_demand_up.balancing;
		output_row(21) = settlement.volume_demand_down.balancing;
		output_row(22) = settlement.cost_supply.balancing;
		output_row(23) = settlement.utility_demand.balancing;
		output_row(24) = settlement.price.balancing;
		output_row(25) = settlement.reimburse.balancing;
		output_row(26) = settlement.volume_supply.BESS;
		output_row(27) = settlement.volume_demand.BESS;

		return output_row;
	}

	void aggregators_settlement_print(power_market::market_whole_inform &Power_market_inform, std::string dir_name){
		int point_num = Power_market_inform.agent_profiles.aggregators.size();
		std::vector <std::string> var_names = var_names_set();
		int var_num = var_names.size();

		Eigen::MatrixXd Output_data = Eigen::MatrixXd::Zero(point_num, var_num);
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			Output_data(point_iter, 0) = point_iter;
			Output_data.row(point_iter) = output_row_store(var_num, Power_market_inform.agent_profiles.aggregators[point_iter].settlement);
		}

		std::string fout_name =  dir_name + "/aggregator.csv";
		basic::write_file(Output_data, fout_name, var_names);
	}

	void cross_border_settlement_print(power_market::market_whole_inform &Power_market_inform, std::string dir_name){
		int edge_num = Power_market_inform.agent_profiles.cross_border.size();
		std::vector <std::string> var_names = var_names_set();
		int var_num = var_names.size();
		std::vector <Eigen::RowVectorXd> stored_information;
		stored_information.reserve(2 * edge_num);

		for(int edge_iter = 0; edge_iter < edge_num; ++ edge_iter){
			int node_num = Power_market_inform.agent_profiles.cross_border[edge_iter].node_num;

			if(node_num == 0){
				continue;
			}
			for(int node_iter = 0; node_iter < node_num; ++ node_iter){
				stored_information.push_back(output_row_store(var_num, Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].settlement));
				stored_information[stored_information.size() - 1](0) = Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].node_ID;
			}
		}

		Eigen::MatrixXd Output_data = Eigen::MatrixXd::Zero(stored_information.size(), var_num);
		for(int row_iter = 0; row_iter < stored_information.size(); ++ row_iter){
			Output_data.row(row_iter) = stored_information[row_iter];
		}
		std::string fout_name =  dir_name + "/cross_border.csv";
		basic::write_file(Output_data, fout_name, var_names);
	}

	void end_users_settlement_print(power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, std::string dir_name){
		int point_num = Power_network_inform.points.bidding_zone.size();
		//int sample_num = agent::end_user::parameters::sample_num();
		int sample_num = Power_market_inform.agent_profiles.end_user_type.sample_num;
		std::vector <std::string> var_names = var_names_set();
		int var_num = var_names.size();

		// Different types of end_users
		for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
			Eigen::MatrixXd Output_data = Eigen::MatrixXd::Zero(point_num, var_num);
			for(int point_iter = 0; point_iter < point_num; ++ point_iter){
				Output_data.row(point_iter) = output_row_store(var_num, Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.settlement);
				Output_data(point_iter, 0) = point_iter;
			}

			std::string fout_name =  dir_name + "/end_user_type_" + std::to_string(sample_iter) + ".csv";
			basic::write_file(Output_data, fout_name, var_names);
		}
	}

	void industrial_settlement_print(power_market::market_whole_inform &Power_market_inform, std::string dir_name){
		std::vector <std::string> var_names = var_names_set();
		int var_num = var_names.size();

		int industrial_HV_num = Power_market_inform.agent_profiles.industrial.HV.size();
		Eigen::MatrixXd Output_data = Eigen::MatrixXd::Zero(industrial_HV_num, var_num);
		for(int agent_iter = 0; agent_iter < industrial_HV_num; ++ agent_iter){
			Output_data.row(agent_iter) = output_row_store(var_num, Power_market_inform.agent_profiles.industrial.HV[agent_iter].settlement);
			Output_data(agent_iter, 0) = Power_market_inform.agent_profiles.industrial.HV[agent_iter].point_ID;
		}

		std::string fout_name =  dir_name + "/industrial_HV.csv";
		basic::write_file(Output_data, fout_name, var_names);
	}

	void power_supplier_settlement_print(power_market::market_whole_inform &Power_market_inform, std::string dir_name){
		std::vector <std::string> var_names = var_names_set();
		int var_num = var_names.size();
		Eigen::MatrixXd Output_data;
		std::string fout_name;

		int hydro_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant.size();
		Output_data = Eigen::MatrixXd::Zero(hydro_HV_plant_num, var_num);
		for(int agent_iter = 0; agent_iter < hydro_HV_plant_num; ++ agent_iter){
			Output_data.row(agent_iter) = output_row_store(var_num, Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].settlement);
			Output_data(agent_iter, 0) = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].original_ID;
		}
		fout_name =  dir_name + "/hydro_HV_plant.csv";
		basic::write_file(Output_data, fout_name, var_names);

		int hydro_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant.size();
		Output_data = Eigen::MatrixXd::Zero(hydro_LV_plant_num, var_num);
		for(int agent_iter = 0; agent_iter < hydro_LV_plant_num; ++ agent_iter){
			Output_data.row(agent_iter) = output_row_store(var_num, Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].settlement);
			Output_data(agent_iter, 0) = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].original_ID;
		}
		fout_name =  dir_name + "/hydro_LV_plant.csv";
		basic::write_file(Output_data, fout_name, var_names);

		int wind_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant.size();
		Output_data = Eigen::MatrixXd::Zero(wind_HV_plant_num, var_num);
		for(int agent_iter = 0; agent_iter < wind_HV_plant_num; ++ agent_iter){
			Output_data.row(agent_iter) = output_row_store(var_num, Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].settlement);
			Output_data(agent_iter, 0) = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].original_ID;
		}
		fout_name =  dir_name + "/wind_HV_plant.csv";
		basic::write_file(Output_data, fout_name, var_names);

		int wind_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant.size();
		Output_data = Eigen::MatrixXd::Zero(wind_LV_plant_num, var_num);
		for(int agent_iter = 0; agent_iter < wind_LV_plant_num; ++ agent_iter){
			Output_data.row(agent_iter) = output_row_store(var_num, Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].settlement);
			Output_data(agent_iter, 0) = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].original_ID;
		}
		fout_name =  dir_name + "/wind_LV_plant.csv";
		basic::write_file(Output_data, fout_name, var_names);

		int pump_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV.size();
		Output_data = Eigen::MatrixXd::Zero(pump_HV_plant_num, var_num);
		for(int agent_iter = 0; agent_iter < pump_HV_plant_num; ++ agent_iter){
			Output_data.row(agent_iter) = output_row_store(var_num, Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].settlement);
			Output_data(agent_iter, 0) = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].original_ID;
		}
		fout_name =  dir_name + "/pump_storage_HV_plant.csv";
		basic::write_file(Output_data, fout_name, var_names);

		int pump_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
		Output_data = Eigen::MatrixXd::Zero(pump_LV_plant_num, var_num);
		for(int agent_iter = 0; agent_iter < pump_LV_plant_num; ++ agent_iter){
			Output_data.row(agent_iter) = output_row_store(var_num, Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].settlement);
			Output_data(agent_iter, 0) = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].original_ID;
		}
		fout_name =  dir_name +  "/pump_storage_LV_plant.csv";
		basic::write_file(Output_data, fout_name, var_names);

		int slack_LV_num = Power_market_inform.agent_profiles.power_supplier.slack.LV_plant.size();
		Output_data = Eigen::MatrixXd::Zero(slack_LV_num, var_num);
		for(int agent_iter = 0; agent_iter < slack_LV_num; ++ agent_iter){
			Output_data.row(agent_iter) = output_row_store(var_num, Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].settlement);
			Output_data(agent_iter, 0) = Power_market_inform.agent_profiles.power_supplier.slack.LV_plant[agent_iter].original_ID;
		}
		fout_name = dir_name + "/slack_LV_plant.csv";
		basic::write_file(Output_data, fout_name, var_names);
	}
}

void agent::agents_results_print(power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, configuration::process_config &process_par){
	// Create a folder to store the file
	std::string dir_name = "csv/case/" + process_par.folder_name + "/output/agent";
	std::filesystem::create_directories(dir_name);

	local::aggregators_settlement_print(Power_market_inform, dir_name);
	local::cross_border_settlement_print(Power_market_inform, dir_name);
	local::end_users_settlement_print(Power_market_inform, Power_network_inform, dir_name);
	local::industrial_settlement_print(Power_market_inform, dir_name);
	local::power_supplier_settlement_print(Power_market_inform, dir_name);
}
