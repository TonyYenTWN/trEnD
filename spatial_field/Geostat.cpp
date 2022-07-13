// Geostat functions
#include "Geostat.h"

// Function definition
inline Eigen::Vector3d xyz_transform(Eigen::Vector2d uv){
	//	x = cos(u) * sin(v)
	//	y = sin(u) * sin(v)
	//	z = cos(v)
	
	Eigen::Vector3d xyz;
	xyz(0) = cos(uv(0)) * sin(uv(1));
	xyz(1) = sin(uv(0)) * sin(uv(1));
	xyz(2) = cos(uv(1));
	
	return xyz;
}

inline Eigen::Vector2d uv_transform(Eigen::Vector3d xyz){
	// u = atan2(y , x) 
	// v = acos(z)
	
	Eigen::Vector2d uv;
	uv(0) = atan2(xyz(1), xyz(0));
	uv(1) = acos(xyz(2));
	
	return uv;
}

inline double geodist(Eigen::Vector2d P_1, Eigen::Vector2d P_2){
	// Constants of Earth Spheroid
	Earth_Constant Geostat_constants;	
	double pi = boost::math::constants::pi<double>();
	
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

// Find the distance and covariance between points
void point_distance_cov(points &point, double lambda){	
	double pi = boost::math::constants::pi<double>();
	
	#pragma omp parallel
	{
		#pragma omp for	
		for(int row_iter = 0; row_iter < point.lon.size() - 1; ++ row_iter){
			for(int col_iter = row_iter + 1; col_iter < point.lon.size(); ++ col_iter){
				point.distance(row_iter, col_iter) = geodist(Eigen::Vector2d(point.lon(row_iter) * pi / 180., point.lat(row_iter) * pi / 180.), Eigen::Vector2d(point.lon(col_iter) * pi / 180., point.lat(col_iter) * pi / 180.));
				point.distance(col_iter, row_iter) = point.distance(row_iter, col_iter);
				point.covariance(row_iter, col_iter) = exp(-pow(point.distance(row_iter, col_iter), 1.) / point.grid_length / lambda);
				point.covariance(col_iter, row_iter) = point.covariance(row_iter, col_iter);
			}
		}
	}
}