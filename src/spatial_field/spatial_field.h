// Header file for estimation of spatial fields

#ifndef SPATIAL_FIELD
#define SPATIAL_FIELD

#include <boost/math/distributions/normal.hpp>
#include "src/basic/rw_csv.h"
#include "src/power_market/power_market.h"
#include "src/power_network/power_network.h"
#include "geostat.h"

namespace spatial_field{
	// Spatial field objects
	struct estimation_inform{
		double alpha_iteration;
		boost::math::normal norm_dist = boost::math::normal(0.0, 1.0);
		Eigen::VectorXd mu_mean;
		Eigen::VectorXd mu_scale;
		Eigen::VectorXd x_scale;
		Eigen::VectorXd mu;
		Eigen::VectorXd x;
		Eigen::SparseMatrix <double> Constraint;
		Eigen::SparseMatrix <double> Conversion_Mat_1;
		Eigen::SparseMatrix <double> Conversion_Mat_2;
	};

	struct fin_field{
		std::string dir;
		std::string demand;
		std::string imbalance;
		std::string solar;
		std::string wind_on;
	};

	// Functions
	void demand_imbalance_estimation(power_network::network_inform&, power_market::market_inform&);
	void wind_on_cf_estimation(power_network::network_inform&);
	void solar_radiation_estimation(power_network::network_inform&);
	void spatial_field_store(power_network::network_inform&, fin_field, int);
	static inline double solar_cf_calculation(double solar_radiation){
		double value = solar_radiation * .0007;
		value = std::min(.7, value);
		value = std::max(0., value);
		return value;
	};
}

#endif
