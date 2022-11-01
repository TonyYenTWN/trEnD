// Header file for power plants
#pragma once
#include "agent_basic.h"

namespace agent{
	namespace power_supplier{
		namespace parameters{
			static inline double cutoff_power(){
				double value = 20.;
				return value;
			}
		}

		struct plant_profile{
			int original_ID;
			int point_ID;
			double cap;
			bids bids;
			results results;
			settlement settlement;
		};

		struct storage_profile{
			int original_ID;
			int point_ID;
			double energy;
			double cap;
			double soc;
			bids bids;
			results results;
			settlement settlement;
		};

		struct hybrid_profile{
			plant_profile plant;
			storage_profile storage;
		};

		struct plants{
			std::vector <plant_profile> HV_plant;
			std::vector <plant_profile> LV_plant;
			std::vector <hybrid_profile> HV_hybrid;
			std::vector <hybrid_profile> LV_hybrid;
		};

		struct storage{
			std::vector <storage_profile> HV;
			std::vector <storage_profile> LV;
		};

		struct profiles{
			plants hydro;
			plants wind;
			plants solar;
			storage pump_storage;
		};
	}
}
