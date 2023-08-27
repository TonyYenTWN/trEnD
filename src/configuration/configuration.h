// Header file for setting the configuration of simulation
#pragma once
#include <iostream>
#include "src/basic/rw_csv.h"

namespace configuration{
    /** @brief Information of configuration for the process of the program*/
	struct process_config{
		/**
		* @name configuration for the initialization
		*/
		/*@{*/
		/** Indicate whether the program should use default configuration.*/
		bool default_flag;
		/** Indicate whether the user wish to manually set the configuration for the main process.*/
		bool set_flag;
		/*@{*/

		/**
		* @name configuration for the estimation step
		*/
		/*@{*/
		/** Indicate whether estimation step should be proceeded. If this flag is false, then no estimation will occur regardless of the values of the other flags.*/
		bool estimation_flag;
		/** Indicate whether estimation of demand should be proceeded.*/
		bool estimation_demand_flag;
		/** Indicate whether estimation of wind cf should be proceeded.*/
		bool estimation_wind_flag;
		/** Indicate whether estimation of solar cf should be proceeded.*/
		bool estimation_solar_flag;
		/*@{*/

		/**
		* @name configuration for the power market simulation step
		*/
		/*@{*/
		/** Indicate whether simulation should be proceeded. If this flag is false, then no simulation will occur regardless of the values of the other flags.*/
		bool simulation_flag;
		/** Indicate whether DSO filtering of bids on their network should occur.*/
		bool DSO_filter_flag;
		/** Indicate whether control reserve calculations should occur.*/
		bool control_reserve_flag;
		/** Indicate whether flexibility from end-users should be prioritized for redispatch.*/
		bool encourage_redispatch;
		/** Indicate whether active end-users use rule-based (and not optimization) when providing flexibility to the system.*/
		bool rule_based;
		/** Total time length of the input time series data (not the simulation time length!!)*/
		int total_time;
		/** Time boundary of simulation. 1st component = starting time; 2nd component = duration of the simulation.*/
		std::vector <int> time_boundary;
		/**  Name of the directory where the csv files are read, processed, and stored*/
		std::string folder_name;
		/*@{*/

		/**
		* @name configuration for the contingency analysis step
		*/
		/*@{*/
		/** Indicate whether a contingency analysis should be carried out. When flex_stat has been processed, this can be done without the simulation step.*/
        bool contingency_flag;
		/*@{*/

		void process_default_get(){
			std::cout << "Default procedure?            Yes: 1 / No: 0 | ";
			std::cin >> this->default_flag;
			std::cout << "\n";
		}

		void process_bool_set(){
		    this->set_flag = 0;
            this->estimation_flag = 0;
            this->estimation_demand_flag = 0;
            this->estimation_wind_flag = 0;
            this->estimation_solar_flag = 0;
            this->simulation_flag = 1;
            this->DSO_filter_flag = 0;
            this->control_reserve_flag = 0;
            this->encourage_redispatch = 0;
            this->rule_based = 0;
            this->total_time = 8760;
            this->time_boundary.push_back(0);
            this->time_boundary.push_back(336);
		}

		// Keep this function so users can manage the configuration file on console
		void process_bool_input(){
			std::cout << "Estimate spatial fields?      Yes: 1 / No: 0 | ";
			std::cin >> this->estimation_flag;
			std::cout << "\n";

			if(this->estimation_flag == 1){
                std::cout << "Estimate demand field?        Yes: 1 / No: 0 | ";
                std::cin >> this->estimation_demand_flag;
                std::cout << "\n";

                std::cout << "Estimate wind cf field?       Yes: 1 / No: 0 | ";
                std::cin >> this->estimation_wind_flag;
                std::cout << "\n";

                std::cout << "Estimate solar cf field?      Yes: 1 / No: 0 | ";
                std::cin >> this->estimation_solar_flag;
                std::cout << "\n";
			}

			std::cout << "Simulate operation?           Yes: 1 / No: 0 | ";
			std::cin >> this->simulation_flag;
			std::cout << "\n";

			if(this->simulation_flag == 1){
				std::cout << "DSOs filter bids?             Yes: 1 / No: 0 | ";
				std::cin >> this->DSO_filter_flag;
				std::cout << "\n";

				std::cout << "Control reserve?              Yes: 1 / No: 0 | ";
				std::cin >> this->control_reserve_flag;
				std::cout << "\n";

				std::cout << "Encourage Redispatch?         Yes: 1 / No: 0 | ";
				std::cin >> this->encourage_redispatch;
				std::cout << "\n";

				std::cout << "Rule based?                   Yes: 1 / No: 0 | ";
				std::cin >> this->rule_based;
				std::cout << "\n";
			}

			std::cout << "Contigency analysis?          Yes: 1 / No: 0 | ";
			std::cin >> this->contingency_flag;
			std::cout << "\n";

			if(this->estimation_flag + this->simulation_flag + this->contingency_flag > 0){
                std::cout << "Total length of input data time series?      | ";
                std::cin >> this->total_time;
                std::cout << "\n";

				std::cout << "Start tick?                                  | ";
				int start_tick;
				std::cin >> start_tick;
				this->time_boundary.push_back(start_tick);
				std::cout << "\n";

				std::cout << "Run time length?                             | ";
				int time_length;
				std::cin >> time_length;
				this->time_boundary.push_back(time_length);
				std::cout << "\n";
			}

			std::cout << "Folder name?                                 | ";
			std::cin >> this->folder_name;
			std::cout << "\n";
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
			std::cout << "\n";

			std::cout << "Contingency analysis?               Yes: 1 / No: 0 | ";
			std::cout << this->contingency_flag;
			std::cout << "\n";

			std::cout << "Folder name?                                 | ";
			std::cout << this->folder_name;
			std::cout << "\n\n";
		}
	};

    // Function for reading configuration data
    void process_config_input(process_config&, std::string);
}
