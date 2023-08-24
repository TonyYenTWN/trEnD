// Header file for setting the configuration of simulation
#pragma once
#include <iostream>
#include "src/basic/basic_definitions.h"

namespace configuration{
	namespace parameters{
		static inline int Time(){
			//int value = 8760;
			int value = 48;
			return value;
		}
	}

	struct process_config{
		bool default_flag;
		bool estimation_flag;
		bool estimation_demand_flag;
		bool estimation_wind_flag;;
		bool estimation_solar_flag;
		bool simulation_flag;
		bool DSO_filter_flag;
		bool control_reserve_flag;
		bool encourage_redispatch;
		bool rule_based;
		std::vector <int> time_boundary;

		void process_default_get(){
			std::cout << "Default procedure?        Yes: 1 / No: 0 | ";
			std::cin >> this->default_flag;
			std::cout << "\n";
		}

		void process_bool_set(){
			this->estimation_flag = 0;
            this->estimation_demand_flag = 0;
            this->estimation_wind_flag = 0;
            this->estimation_solar_flag = 0;
			this->simulation_flag = 1;
			this->DSO_filter_flag = 0;
			this->control_reserve_flag = 0;
			this->encourage_redispatch = 0;
			this->rule_based = 0;
			this->time_boundary.push_back(0);
			this->time_boundary.push_back(168);
		}

		void process_bool_input(){
			std::cout << "Estimate spatial fields?  Yes: 1 / No: 0 | ";
			std::cin >> this->estimation_flag;
			std::cout << "\n";

			if(this->estimation_flag == 1){
                std::cout << "Estimate demand field?    Yes: 1 / No: 0 | ";
                std::cin >> this->estimation_demand_flag;
                std::cout << "\n";

                std::cout << "Estimate wind cf field?   Yes: 1 / No: 0 | ";
                std::cin >> this->estimation_wind_flag;
                std::cout << "\n";

                std::cout << "Estimate solar cf field?  Yes: 1 / No: 0 | ";
                std::cin >> this->estimation_solar_flag;
                std::cout << "\n";
			}

			std::cout << "Simulate operation?       Yes: 1 / No: 0 | ";
			std::cin >> this->simulation_flag;
			std::cout << "\n";

			if(this->simulation_flag == 1){
				std::cout << "DSOs filter bids?         Yes: 1 / No: 0 | ";
				std::cin >> this->DSO_filter_flag;
				std::cout << "\n";

				std::cout << "Control reserve?          Yes: 1 / No: 0 | ";
				std::cin >> this->control_reserve_flag;
				std::cout << "\n";

				std::cout << "Encourage Redispatch?     Yes: 1 / No: 0 | ";
				std::cin >> this->encourage_redispatch;
				std::cout << "\n";

				std::cout << "Rule based?               Yes: 1 / No: 0 | ";
				std::cin >> this->rule_based;
				std::cout << "\n";
			}

			if(this->estimation_flag || this->simulation_flag){
				std::cout << "Start tick?                              | ";
				int start_tick;
				std::cin >> start_tick;
				this->time_boundary.push_back(start_tick);
				std::cout << "\n";

				std::cout << "Run time length?                         | ";
				int time_length;
				std::cin >> time_length;
				this->time_boundary.push_back(time_length);
				std::cout << "\n";
			}
		}

		void process_bool_output(){
			std::cout << "Estimate spatial fields?  Yes: 1 / No: 0 | ";
			std::cout << this->estimation_flag;
			std::cout << "\n";

			std::cout << "Estimate demand field?  Yes: 1 / No: 0 | ";
			std::cout << this->estimation_demand_flag;
			std::cout << "\n";

			std::cout << "Estimate wind field?  Yes: 1 / No: 0 | ";
			std::cout << this->estimation_wind_flag;
			std::cout << "\n";

			std::cout << "Estimate solar field?  Yes: 1 / No: 0 | ";
			std::cout << this->estimation_solar_flag;
			std::cout << "\n";

			std::cout << "Simulate operation?       Yes: 1 / No: 0 | ";
			std::cout << this->simulation_flag;
			std::cout << "\n";

			std::cout << "DSOs filter bids?         Yes: 1 / No: 0 | ";
			std::cout << this->DSO_filter_flag;
			std::cout << "\n";

			std::cout << "Control reserve?          Yes: 1 / No: 0 | ";
			std::cout << this->control_reserve_flag;
			std::cout << "\n";

			std::cout << "Encourage Redispatch?     Yes: 1 / No: 0 | ";
			std::cout << this->encourage_redispatch;
			std::cout << "\n";

			std::cout << "Rule based?               Yes: 1 / No: 0 | ";
			std::cout << this->rule_based;
			std::cout << "\n";

			std::cout << "Start tick?                              | ";
			std::cout << this->time_boundary[0];
			std::cout << "\n";

			std::cout << "Run time length?                         | ";
			std::cout << this->time_boundary[1];
			std::cout << "\n\n";
		}
	};
}
