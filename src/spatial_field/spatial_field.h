// Header file for inference of spatial fields

#ifndef SPATIAL_FIELD
#define SPATIAL_FIELD

#include <boost/math/distributions/normal.hpp>
#include "src/basic/rw_csv.h"
#include "src/power_market/power_market.h"
#include "src/power_network/power_network.h"
#include "geostat.h"

namespace spatial_field{
	void BME_area(power_network::network_inform&, Eigen::SparseMatrix <double>&, Eigen::MatrixXd&);

	Eigen::VectorXd nominal_demand_inference(power_network::network_inform&);

	void spatial_field_store(power_network::network_inform&, std::string, int);
}

#endif
