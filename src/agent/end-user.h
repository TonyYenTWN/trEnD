// Header file for operation and investment of end-user
#pragma once
#include "agent.h"

namespace agent{
	namespace end_user{
		/** @brief Decision variables representing possible investment combinations of an end-user.*/
		struct decision{
			bool dynamic_tariff;
			bool smart_appliance;
			bool PV_BESS;
			bool EV_self_charging;
			/** Whether the end-user can inject power flow back to grid; false when active_flex is false*/
			bool reverse_flow;
			/** Whether the end-user can provide flexibility to the aggregator; false when dynamic_tariff is false, or when end-user does not have PV + BESS, EV, or smart appliance*/
			bool active_flex;
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
			double scale;
			/**  Maximum load shifting time length of a flexible demand;
			* indicates how flexible the smart appliances are.
			*/
			int shift_time;
//			/**  Indicates how flexible the smart appliances are;
//			* e.g. 1 / 2 = can concentrate the demand within half of the time interval.
//			*/
//			double flexibility_factor;
			/*@{*/


			/**
			* @name process variables
			*/
			/*@{*/
			/** Flexible demand not yet fulfilled in the load shifting timeframe (kWh per hour per person).*/
			Eigen::VectorXd unfulfilled_demand;
			/** Scheduled flexible demand at current time step (kWh per hour per person).*/
			double scheduled_demand;
			/*@{*/

//			/**
//			* @name output variables
//			*/
//			/*@{*/
//			/** Normalized schedule profile for smart appliance (kWh per hour per person).*/
//			Eigen::VectorXd normalized_scheduled_profile;
//			/*@{*/
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
//			/** Initial state of charge level.*/
//			double soc_ini;
//			/** Final state of charge level.*/
//			double soc_final;
			/*@{*/

			/**
			* @name process parameters
			*/
			/*@{*/
			/** State of charge level.*/
			double soc;
			/** Scheduled charge / discharge output at current time step (kWh per hour per person).*/
			double scheduled_capacity;
//			alglib::minlpstate Problem;
			/*@{*/

//			/**
//			* @name output variables
//			*/
//			/*@{*/
//			/** Normalized schedule profile for battery storage capacity (kWh per hour per person).*/
//			Eigen::VectorXd normalized_scheduled_capacity_profile;
//			/** Normalized schedule profile for battery storage state of charge level (kWh).*/
//			Eigen::VectorXd normalized_scheduled_soc_profile;
//			/*@{*/
		};

		/** @brief Information of electric vehicle of an end-user.*/
		struct EV_inform{
			/**
			* @name input parameters
			*/
			/*@{*/
			/**kWh per person per hour of usage.*/
			double energy_demand;
			/** The time intervals when EV is actually used.*/
			Eigen::VectorXi usage_default_period;
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
			/** Default demand profile; normalized to nominal value (kWh per hour per person)*/
			Eigen::VectorXd default_demand_profile;
			/** Default PV profile; normalized to nominal value (kWh per hour per person).*/
			Eigen::VectorXd default_PV_profile;
			/*@{*/

			/**
			* @name process variables
			*/
			/*@{*/
			alglib::minlpstate Problem;
			bids bids;
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

		void smart_appliance_schedule(sorted_vector, Eigen::VectorXd, smart_appliance_inform&);
		alglib::minlpstate storage_schedule_LP_mold(int);
		void storage_schedule_LP_optimize(int, sorted_vector, storage_inform&, bool fixed_end = 0);
		void EV_schedule(int, sorted_vector, EV_inform&);
	}
}
