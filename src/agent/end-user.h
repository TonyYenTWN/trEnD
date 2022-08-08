// Header file for operation and investment of end-user
#include "src/alglib/optimization.h"
#include "src/basic/basic_definitions.h"
#include "src/basic/eigen_sparse.h"

namespace agent{
	namespace parameters{
		static inline int foresight_time(){
			int value = 24;
			return value;
		}

		static inline int sample_num(){
			int value = 3;
			return value;
		}

		static inline double residential_ratio(){
			double value = 1.;
			return value;
		}
	}

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
			/**  Indicates how flexible the smart appliances are;
			* e.g. 1 / 2 = can concentrate the demand within half of the time interval.
			*/
			double flexibility_factor;
			/*@{*/

			/**
			* @name process variables
			*/
			/*@{*/
			/** Previous demand not yet fulfilled (kWh per hour per person).*/
			double unfulfilled_demand;
			/*@{*/

			/**
			* @name output variables
			*/
			/*@{*/
			/** Normalized schedule profile for smart appliance (kWh per hour per person).*/
			Eigen::VectorXd normalized_scheduled_profile;
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
			/** Initial state of charge level.*/
			double soc_ini;
			/** Final state of charge level.*/
			double soc_final;
			/*@{*/

			/**
			* @name process parameters
			*/
			/*@{*/
			alglib::minlpstate Problem;
			/*@{*/

			/**
			* @name output variables
			*/
			/*@{*/
			/** Normalized schedule profile for battery storage capacity (kWh per hour per person).*/
			Eigen::VectorXd normalized_scheduled_capacity_profile;
			/** Normalized schedule profile for battery storage state of charge level (kWh).*/
			Eigen::VectorXd normalized_scheduled_soc_profile;
			/*@{*/
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
			/** Weight of the sample representing the population distribution.
			* Summation at a sample point equals to 1.*/
			double weight;
			/** Scale of the PV system = capacity / normalized default base value of demand*/
			double PV_scale;
			/** Default demand profile; normalized to nominal value (kWh per hour per person)*/
			Eigen::VectorXd normalized_default_demand_profile;
			/** Default PV profile; normalized with the same base as demand.*/
			Eigen::VectorXd normalized_default_PV_profile;
			/*@{*/

			/**
			* @name output variables
			*/
			/*@{*/
			double supply_inflex_price;
			double supply_flex_price;
			double demand_inflex_price;
			double demand_flex_price;
			Eigen::VectorXd normalized_scheduled_residual_demand_inflex_profile;
			Eigen::VectorXd normalized_confirmed_residual_demand_inflex_profile;
			Eigen::VectorXd normalized_actual_residual_demand_inflex_profile;
			Eigen::VectorXd normalized_scheduled_residual_demand_flex_profile;
			Eigen::VectorXd normalized_confirmed_residual_demand_flex_profile;
			Eigen::VectorXd normalized_actual_residual_demand_flex_profile;
			//Eigen::VectorXd normalized_scheduled_pos_cr_profile;
			//Eigen::VectorXd normalized_scheduled_neg_cr_profile;
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
	}

	struct sorted_vector{
		Eigen::VectorXi id;
		Eigen::VectorXd value;
	};

	// Functions
	namespace end_user{
		void smart_appliance_schedule(sorted_vector, Eigen::VectorXd, smart_appliance_inform&);
		alglib::minlpstate storage_schedule_LP_mold(int);
		void storage_schedule_LP_optimize(int, sorted_vector, storage_inform&, bool fixed_end = 0);
		void EV_schedule(int, sorted_vector, EV_inform&);
	}
	sorted_vector sort(Eigen::VectorXd);
}
