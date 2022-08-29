// Header file for inference of spatial fields

#ifndef SPATIAL_FIELD
#define SPATIAL_FIELD

#include <boost/math/distributions/normal.hpp>
#include "src/basic/rw_csv.h"
#include "src/power_market/power_market.h"
#include "src/power_network/power_network.h"
#include "geostat.h"

namespace spatial_field{
	// Spatial field objects
	struct inference_inform{
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

	// Functions
	//void BME_copula(inference_inform&, power_network::network_inform&, Eigen::SparseMatrix <double>&, double);
	//void BME_linear(inference_inform&, Eigen::SparseMatrix <double>&);
	//void imbalance_inference(inference_inform &, power_network::network_inform&);
	void spatial_field_inference(power_network::network_inform&);
	void wind_on_cf_inference(power_network::network_inform&);
	void solar_radiation_inference(power_network::network_inform&);
	void spatial_field_store(power_network::network_inform&, std::string, std::string, std::string, int);
}

#endif
