// Header file for inference of spatial fields

#ifndef SPATIAL_FIELD
#define SPATIAL_FIELD

#include <boost/math/distributions/normal.hpp>
#include "src/basic/rw_csv.h"
#include "src/power_network/power_network.h"
#include "geostat.h"

namespace spatial_field{
	Eigen::VectorXd BME(power_network::network_inform&);
	//Eigen::SparseMatrix <double>

	Eigen::VectorXd nominal_demand_inference();
}

#endif
