// Header file for setting the configuration of simulation
#pragma once
#include <iostream>

namespace configuration{
	namespace parameters{
		static inline int Time(){
			int value = 8760;
			return value;
		}
	}

	struct process_config{
		bool default_flag;
		bool estimation_flag;
		bool simulation_flag;
		bool DSO_filter_flag;
		bool control_reserve_flag;
		std::vector <int> time_boundary;

		void process_default_get(){
			std::cout << "Default procedure?        Yes: 1 / No: 0 | ";
			std::cin >> this->default_flag;
			std::cout << "\n";
		}

		void process_bool_set(){
			this->estimation_flag = 0;
			this->simulation_flag = 1;
			this->DSO_filter_flag = 0;
			this->control_reserve_flag = 1;
			this->time_boundary.push_back(0);
			this->time_boundary.push_back(168);
		}

		void process_bool_input(){
			std::cout << "Estimate spatial fields?  Yes: 1 / No: 0 | ";
			std::cin >> this->estimation_flag;
			std::cout << "\n";

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
			}

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
	};
}
