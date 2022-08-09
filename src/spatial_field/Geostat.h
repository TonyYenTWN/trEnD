// Geostat functions header file
#pragma once

#ifndef GEOSTAT
#define GEOSTAT

#include "src/basic/basic_definitions.h"

namespace spatial_field{
	/** Constants of for the Earth spheroid*/
	class Earth_Constant{
		public:
			/** Short radius of Earth (in meters).*/
			double Earth_Radius_short = 6356752.3142;
			/** Long radius of Earth (in meters).*/
			double Earth_Radius_long = 6378137.;
			/** Flatness of Earth.*/
			double Earth_flat = 1 / 298.257223563;
			/** Eccentricity of Earth.*/
			double Earth_Ecc = pow(2 * Earth_flat - pow(Earth_flat, 2), .5);
	};

	// Functions for calculation of the geostatistic
	double geodist(Eigen::Vector2d, Eigen::Vector2d);
}

#endif
