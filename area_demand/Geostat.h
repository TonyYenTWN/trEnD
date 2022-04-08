// Geostat functions header file
#pragma once

#ifndef GEOSTAT
#define GEOSTAT

#include "Basic_Definitions.h"

// Constants of Earth Spheroid
class Earth_Constant{
	public:
		double Earth_Radius_short = 6356752.3142;							// long axis
	    double Earth_Radius_long = 6378137;									// short axis
	    double Earth_flat = 1 / 298.257223563; 								// flatness
	    double Earth_Ecc = pow(2 * Earth_flat - pow(Earth_flat, 2), .5);	// eccentricity
};

Vector3d xyz_transform(Vector2d);
Vector2d uv_transform(Vector3d);
double geodist(Vector2d, Vector2d);

#endif