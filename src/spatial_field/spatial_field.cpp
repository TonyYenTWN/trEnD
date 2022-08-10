// Main source file for inference of spatial fields

#include "src/spatial_field/spatial_field.h"

void spatial_field::BME_area(power_network::network_inform &Power_network_inform, Eigen::SparseMatrix <double> &Constraint, Eigen::MatrixXd &x_ts){
	int bz_num = Power_network_inform.points.bidding_zone.maxCoeff() + 1;
	int point_num = Power_network_inform.points.bidding_zone.size();
	int Time = power_market::parameters::Time();
	double pi = boost::math::constants::pi<double>();
	double tol = 1E-12;
	boost::math::normal norm_dist(0.0, 1.0);

	// Initialization of Parameters
	Eigen::VectorXd x_0_mean(bz_num);
	x_0_mean = x_ts.colwise().sum().tail(bz_num);
	x_0_mean /= Time;
	double multiply_factor = 1.;
	x_0_mean *= multiply_factor;
	// Mean of weibull distribution = mu_0 * gamma(1 + 1 / k), here k = 2
	double mu_0_mean = x_0_mean.sum() / Constraint.sum() * 2. / pow(pi, .5);

	// ---------------------------------------------------------------------
	// Infer temporal average of the normalized mean demand field
	// ---------------------------------------------------------------------
	// MaxEnt, known covariance and mean = 0, reference equiprobable space = x
	// COV^(-1) * x* - Constraint * lambda = 0
	// x* = COV * Constraint * lambda
	// Constraint^T * x* = Constraint^T * COV * Constraint * lambda = demand
	// lambda = solve(Constraint^T * COV * Constraint, demand)
	//
	// Weibull Distribution using Gaussian Copula
	// ---------------------------------------------------------------------
	// Initialization of constants, vectors, and matrices
	int count;
	double alpha_iteration = 1.;
	Eigen::VectorXd mu_0 = Eigen::VectorXd::Constant(point_num, mu_0_mean);
	Eigen::VectorXd x_0(point_num);
	for(int item = 0; item < point_num; ++ item){
		x_0(item) = quantile(norm_dist, 1. - exp(-pow(mu_0(item) / mu_0_mean, 2.)));
	}
	Eigen::VectorXd mu_inv_0(point_num);
	Eigen::VectorXd lambda(bz_num);

//	// Iterations
//	count = 0;
//	while(count < 5000 && dmu_0.lpNorm<Eigen::Infinity>() > tol){
//		Eigen::VectorXd Conversion_vec(point_num);
//
//		for(int item = 0; item < point_num; ++ item){
//			Conversion_vec(item) = .5 * mu_0_mean * pow((-log(1. - cdf(norm_dist, x_0(item)))), -.5) * pow(1. - cdf(norm_dist, x_0(item)), -1.) * pdf(norm_dist, x_0(item));
//		}
//		Eigen::SparseMatrix <double> Conversion_Mat_1 = Conversion_vec.asDiagonal();
//		Eigen::SparseMatrix <double> Conversion_Mat_2 = Conversion_Mat_1 * Constraint;
//
//		lambda = (Conversion_Mat_2.transpose() * Covariance_Points * Conversion_Mat_2).colPivHouseholderQr().solve(Demand_0 - Constraint.transpose() * mu_0 + Conversion_Mat_2.transpose() * x_0);
//		Eigen::VectorXd dx_0 = Covariance_Points * Conversion_Mat_2 * lambda - x_0;
//		x_0 += alpha_iteration * dx_0;
//
//		Eigen::VectorXd dmu_0(point_num);
//		for(int item = 0; item < point_num; ++ item){
//			dmu_0(item) = mu_0_mean * pow(-log(1. - cdf(norm_dist, x_0(item))), .5) - mu_0(item);
//			mu_0(item) += dmu_0(item);
//		}
//		mu_inv_0 = pow(mu_0.array(), -1.);
//		dmu_0 *= mu_inv_0;
//
//		count += 1;
//	}
}

Eigen::VectorXd spatial_field::nominal_demand_inference(power_network::network_inform &Power_network_inform){
	int bz_num = Power_network_inform.points.bidding_zone.maxCoeff() + 1;
	int point_num = Power_network_inform.points.bidding_zone.size();

	// Read electricity demand data;
	std::string fin_demand = "csv/input/spatial_field/demand_actual_2021.csv";
	auto fin_demand_dim = basic::get_file_dim(fin_demand);
	auto Demand_ts = basic::read_file(fin_demand_dim[0], fin_demand_dim[1], fin_demand);

	// Set sparse matrix for equality constaints
	Eigen::SparseMatrix <double> Constraint (point_num, bz_num);
	std::vector<Eigen::TripletXd> Constraint_Trip;
	Constraint_Trip.reserve(point_num);

	for(int point_iter = 0; point_iter < point_num; ++ point_iter){
		int bz_ID = Power_network_inform.points.bidding_zone(point_iter);
		Constraint_Trip.push_back(Eigen::TripletXd(point_iter, bz_ID, Power_network_inform.points.population_density(point_iter)));
	}

	Constraint.setFromTriplets(Constraint_Trip.begin(), Constraint_Trip.end());

	// Spatial field inference
	BME_area(Power_network_inform, Constraint, Demand_ts);

//	Demand_0 = Demand_ts.colwise().sum().tail(bz_num);
//	Demand_0 /= Time;
//	double multiply_factor = 1.;
//	Demand_0 *= multiply_factor;
//	// Mean of weibull distribution = mu_0 * gamma(1 + 1 / k), here k = 2
//	double mu_0_mean = Demand_0.sum() / Constraint.sum() * 2. / pow(pi, .5);

	Eigen::VectorXd Demand_0(bz_num);
	return Demand_0;
}

// Function that stores processed mean demand field
void spatial_field::spatial_field_store(power_network::network_inform &Power_network_inform, std::string fin_demand, int Time){
	int row_num = Power_network_inform.points.bidding_zone.rows();
	Power_network_inform.points.nominal_mean_demand_field = Eigen::MatrixXd(row_num, Time);

	//for(int tick = 0; tick < Time; ++ tick){
	for(int tick = 0; tick < 100; ++ tick){
		// Find zeros before the number
		int count_zeros = 0;
		int tick_temp = tick;
		std::string digit_zeros;
		while(int (tick_temp / 10) != 0){
			count_zeros += 1;
			tick_temp /= 10;
		}
		for(int item = 0; item < 5 - count_zeros; ++item){
			digit_zeros += std::to_string(0);
		}

		// File name with enumeration
		std::string fin_demand_temp = fin_demand + digit_zeros + std::to_string(tick) + ".csv";

		Power_network_inform.points.nominal_mean_demand_field.col(tick) = basic::read_file(row_num, 1, fin_demand_temp);
	}
}
