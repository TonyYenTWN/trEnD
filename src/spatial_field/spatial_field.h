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
	};

	// Functions
	void BME_copula(inference_inform&, power_network::network_inform&, Eigen::SparseMatrix <double>&, double);
	void BME_linear(inference_inform&, power_network::network_inform&, Eigen::SparseMatrix <double>&, double);
	void nominal_demand_inference(power_network::network_inform&);
	void imbalance_inference(power_network::network_inform&);
	void spatial_field_store(power_network::network_inform&, std::string, int);
}

#endif
