// Geostat functions
#include "Geostat.h"

// Constants of Earth Spheroid
Earth_Constant Geostat_constants;

// Function definition
Eigen::Vector3d xyz_transform(Eigen::Vector2d uv){
	//	x = cos(u) * sin(v)
	//	y = sin(u) * sin(v)
	//	z = cos(v)
	
	Eigen::Vector3d xyz;
	xyz(0) = cos(uv(0)) * sin(uv(1));
	xyz(1) = sin(uv(0)) * sin(uv(1));
	xyz(2) = cos(uv(1));
	
	return xyz;
}

Eigen::Vector2d uv_transform(Eigen::Vector3d xyz){
	// u = atan2(y , x) 
	// v = acos(z)
	
	Eigen::Vector2d uv;
	uv(0) = atan2(xyz(1), xyz(0));
	uv(1) = acos(xyz(2));
	
	return uv;
}

double geodist(Eigen::Vector2d P_1, Eigen::Vector2d P_2){
	
	Eigen::Vector2d s_1 = P_1;
	Eigen::Vector2d s_2 = P_2;
	s_1(1) = pi / 2 - P_1(1);
	s_2(1) = pi / 2 - P_2(1);
	Eigen::Vector3d x_1 = xyz_transform(s_1);
	Eigen::Vector3d x_2 = xyz_transform(s_2);
	double theta = acos(x_1.dot(x_2));
	
	double lat_mean = .5 * (P_1(1) + P_2(1));
	double radius_mean = Geostat_constants.Earth_Radius_short / pow(1 - pow(Geostat_constants.Earth_Ecc * sin(lat_mean), 2), .5);

	return(radius_mean * theta);
}