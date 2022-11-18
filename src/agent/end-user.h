// Header file for operation and investment of end-user
#pragma once
#include "agent_basic.h"

namespace agent{
	namespace end_user{
		namespace parameters{
			static inline int foresight_time(){
				int value = 24;
				return value;
			}

			static inline int load_shift_time(){
				int value = 1;
				return value;
			}

			static inline int sample_num(){
				int value = 1;
				return value;
			}

			static inline double power_factor(){
				double value = 1.;
				return value;
			}
		}

		/** @brief Decision variables representing possible investment combinations of an end-user.*/
		struct decision{
			bool dynamic_tariff;
			bool smart_management; // ideal optimization if true
			bool smart_appliance;
			bool PV;
			bool BESS;
			bool EV_self_charging;
			/** Whether the end-user can inject power flow back to grid*/
			bool reverse_flow;
			/** Whether the end-user can provide flexibility for redispatch; false when dynamic_tariff is false, or when end-user does not have PV + BESS, EV, or smart appliance*/
			bool redispatch;
			/** Whether the end-user can provide flexibility for control reserve; false when dynamic_tariff is false, or when end-user does not have PV + BESS, EV, or smart appliance*/
			bool control_reserve;
		};

		/** @brief Information of smart appliance of an end-user.*/
		struct smart_appliance_inform{
			/**
			* @name input parameters
			*/
			/*@{*/
			/** Indicates how much ratio of the total demand in the time interval can be shifted around flexibly;
			* assuming the default is constant profile before shifting.
			*/
			double scale = .1;
			/**  Maximum load shifting time length of a flexible demand;
			* indicates how flexible the smart appliances are.
			*/
			int shift_time;
			/*@{*/


			/**
			* @name process variables
			*/
			/*@{*/
			/** Flexible demand not yet fulfilled in the load shifting timeframe (kWh per hour per person).*/
			Eigen::VectorXd unfulfilled_demand;
			/** Scheduled flexible demand at current time step (kWh per hour per person).*/
			Eigen::VectorXd scheduled_demand;
			/*@{*/
		};

		/** @brief Information of battery storage system of an end-user.*/
		struct storage_inform{
			/**
			* @name input parameters
			*/
			/*@{*/
			/** Scale of energy level of storage (kWh per person).*/
			double energy_scale = .01;
			/** Scale of capacity level of storage (kW per person).*/
			double capacity_scale = .001;
			/** Conversion efficiency of charge / discharge.*/
			double efficiency = .95;
			/** Self-consumption (kWh per hour per person) of the BESS.*/
			double self_consumption = 0.;
			/*@{*/

			/**
			* @name process parameters
			*/
			/*@{*/
			/** State of charge level.*/
			double soc;
			/** Scheduled charge / discharge output at current time step (kWh per hour per person).*/
			double scheduled_capacity;
			/*@{*/
		};

		/** @brief Information of electric vehicle of an end-user.*/
		struct EV_inform{
			/**
			* @name input parameters
			*/
			/*@{*/
			/**kWh per person per hour of usage.*/
			double energy_demand = 7.671;
//			/** The time intervals when EV is actually used.*/
//			Eigen::VectorXi usage_default_period;
			/** The time intervals when EV is parked in the house.*/
			Eigen::VectorXi house_default_period;
			/** The time series of the default charging demand of the EV (kWh per person).*/
			Eigen::VectorXd default_demand_profile;
			/*@{*/

			/**
			* @name process parameters
			*/
			/*@{*/
			/** State of charge level.*/
			double soc;
			/** Scheduled charge / discharge output at current time step (kWh per hour per person).*/
			double scheduled_capacity;
			/*@{*/

			/**
			* @name mixed substructure
			*/
			/*@{*/
			storage_inform BESS;
			/*@{*/

			Eigen::VectorXi house_schedule(int tick){
				int foresight_time = agent::end_user::parameters::foresight_time();
				Eigen::VectorXi vec = Eigen::VectorXi::Zero(foresight_time);

				for(int tock = 0; tock < foresight_time; ++ tock){
					if((tick + tock) % foresight_time < 7 || (tick + tock) % foresight_time > 19){
						vec(tock) = 1;
					}
				}

				return vec;
			}

			Eigen::VectorXd demand_profile(int tick){
				int foresight_time = agent::end_user::parameters::foresight_time();
				Eigen::VectorXd vec = Eigen::VectorXd::Zero(foresight_time);

				for(int tock = 0; tock < foresight_time; ++ tock){
					if((tick + tock) % foresight_time == 7 || (tick + tock) % foresight_time == 19){
						vec(tock) = this->energy_demand;
					}
				}

				return vec;
			}
		};

		/** @brief Information of the investment strategies of an end-user.*/
		struct investment{
			// Input parameters
			decision decision;

			// Some social factors with complex systems modeling
		};

		/** @brief Information of the operation strategies of an end-user.*/
		struct operation{
			/**
			* @name input parameters
			*/
			/*@{*/
			/** Foresight time length for the end-user when schedule its demand*/
			int foresight_time;
			/** Weight of the sample representing the population distribution.
			* Summation at a sample point equals to 1.*/
			double weight;
			/** Scale of the PV system (kW / person)*/
			double PV_scale;
			/** Actual PV profile */
			double PV_output;
			/** Default demand profile; normalized to nominal value (kWh per hour per person)*/
			Eigen::VectorXd default_demand_profile;
			/** Default PV profile; normalized to nominal value (kWh per hour per person).*/
			Eigen::VectorXd default_PV_profile;
			/** Tariff for consumption.*/
			Eigen::VectorXd price_demand_profile;
			/** Tariff for feed-in supply.*/
			Eigen::VectorXd price_supply_profile;
			/*@{*/

			/**
			* @name process variables
			*/
			/*@{*/
			alglib::minlpstate Problem;
			bids bids;
			/*@{*/

			/**
			* @name output results
			*/
			/*@{*/
			results results;
			settlement settlement;
			/*@{*/

			/**
			* @name mixed substructure
			*/
			/*@{*/
			smart_appliance_inform smart_appliance;
			storage_inform BESS;
			EV_inform EV;
			/*@{*/
		};

		/** @brief Complete profile of an end-user.*/
		struct profile{
			investment investment;
			operation operation;
		};

		/** @brief A 2D vector of end-user profiles. The 1st dimension represents the spatial points,
		* the 2nd represents the samples at the point.*/
		typedef std::vector <std::vector <profile>> profiles;

		// Functions
		void end_user_LP_set(profile&);
		void end_user_LP_optimize(int tick, profile&);
	}
}
