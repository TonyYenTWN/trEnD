// Main source file for inference of spatial fields

#include "src/spatial_field/spatial_field.h"

void spatial_field::BME(power_network::network_inform &Power_network_inform, Eigen::SparseMatrix <double> &Constraint, Eigen::MatrixXd &mu_ts){
	int bz_num = Power_network_inform.points.bidding_zone.maxCoeff() + 1;
	int point_num = Power_network_inform.points.bidding_zone.size();
	int Time = power_market::parameters::Time();
	double pi = boost::math::constants::pi<double>();
	double tol = 1E-12;
	boost::math::normal norm_dist(0.0, 1.0);

	// Initialization of Parameters
	Eigen::VectorXd mu_0_mean(bz_num);
	mu_0_mean = mu_ts.colwise().sum().tail(bz_num);
	mu_0_mean /= Time;
	double multiply_factor = 1.;
	mu_0_mean *= multiply_factor;
	// Mean of weibull distribution = mu_0_scale * gamma(1 + 1 / k), here k = 2
	double mu_0_scale = mu_0_mean.sum() / Constraint.sum() * 2. / pow(pi, .5);

	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Infer temporal average of the normalized mean demand field
	// ------------------------------------------------------------------------------------------------------------------------------------------
	//
	// Weibull Distribution using Gaussian Copula
	// ------------------------------------------------------------------------------------------------------------------------------------------

	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Infer the annual average of the normalized mean demand field
	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Initialization of constants, vectors, and matrices
	double alpha_iteration = 1.;
	Eigen::VectorXd mu_0 = Eigen::VectorXd::Constant(point_num, mu_0_scale);
	Eigen::VectorXd x_0(point_num);
	for(int item = 0; item < point_num; ++ item){
		x_0(item) = quantile(norm_dist, 1. - exp(-pow(mu_0(item) / mu_0_scale, 2.)));
	}

	// Iterations
	int count = 0;
	Eigen::VectorXd dmu_0 = Eigen::VectorXd::Constant(point_num, 1.);
	while(count < 5000 && dmu_0.lpNorm<Eigen::Infinity>() > tol){
		// Matrices for the linear system
		Eigen::SparseMatrix <double> Conversion_Mat_1(point_num, point_num);
		std::vector <Eigen::TripletXd> Conversion_Mat_1_Trip;
		Conversion_Mat_1.reserve(point_num);

		for(int item = 0; item < point_num; ++ item){
			double diff_x = .5 * mu_0_scale * pow((-log(1. - cdf(norm_dist, x_0(item)))), -.5) * pow(1. - cdf(norm_dist, x_0(item)), -1.) * pdf(norm_dist, x_0(item));
			Conversion_Mat_1_Trip.push_back(Eigen::TripletXd(item, item, diff_x));
		}
		Conversion_Mat_1.setFromTriplets(Conversion_Mat_1_Trip.begin(), Conversion_Mat_1_Trip.end());

		Eigen::SparseMatrix <double> Conversion_Mat_2 = Conversion_Mat_1 * Constraint;

		// Solve the linear system
		Eigen::SparseLU <Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering<int>> solver;
		solver.compute(Conversion_Mat_2.transpose() * Power_network_inform.points.covariance * Conversion_Mat_2);
		Eigen::VectorXd lambda = solver.solve(mu_0_mean - Constraint.transpose() * mu_0 + Conversion_Mat_2.transpose() * x_0);
		Eigen::VectorXd dx_0 = Power_network_inform.points.covariance * Conversion_Mat_2 * lambda - x_0;
		x_0 += alpha_iteration * dx_0;

		// Update the solution
		for(int item = 0; item < point_num; ++ item){
			dmu_0(item) = mu_0_scale * pow(-log(1. - cdf(norm_dist, x_0(item))), .5) - mu_0(item);
			mu_0(item) += dmu_0(item);
		}
		Eigen::VectorXd mu_inv_0 = pow(mu_0.array(), -1.);
		dmu_0 *= mu_inv_0;

		count += 1;
	}

	// Output the annual average of normalized mean demand field
	std::string fout_name;
	fout_name = "csv/processed/spatial_field/nominal_mean_demand_field_10km_annual_mean.csv";
	std::vector <std::string> col_name;
	col_name.push_back("nominal_mean_demand");
	basic::write_file(mu_0, fout_name, col_name);

	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Infer the normalized mean demand field
	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Initialization of Constants, Vectors and Matrices
	alpha_iteration = .01;
	Eigen::VectorXd mu_scale = mu_0 * 2. / pow(pi, .5);
	Eigen::VectorXd x_scale(point_num);

	for(int item = 0; item < point_num; ++ item){
		x_scale(item) = quantile(norm_dist, 1. - exp(-pow(1., 2.)));
	}

	Eigen::VectorXd mu = mu_scale;
	Eigen::VectorXd x = x_scale;

	// Run the algorithm for each time slice
	for(int tick = 0; tick < 1; ++ tick){
		Eigen::VectorXd mu_mean = mu_ts.row(tick).tail(bz_num);
		Eigen::VectorXd dmu = Eigen::VectorXd::Constant(point_num, 1.);

		while(count < 5000 && dmu.lpNorm<Eigen::Infinity>() > pow(10., -3.)){
			// Matrices for the linear system
			Eigen::SparseMatrix <double> Conversion_Mat_1(point_num, point_num);
			std::vector <Eigen::TripletXd> Conversion_Mat_1_Trip;
			Conversion_Mat_1.reserve(point_num);

			for(int item = 0; item < point_num; ++ item){
				double diff_x = .5 * mu_scale(item) * pow((-log(1. - cdf(norm_dist, x(item)))), -.5) * pow(1. - cdf(norm_dist, x(item)), -1.) * pdf(norm_dist, x(item));
				Conversion_Mat_1_Trip.push_back(Eigen::TripletXd(item, item, diff_x));
			}
			Conversion_Mat_1.setFromTriplets(Conversion_Mat_1_Trip.begin(), Conversion_Mat_1_Trip.end());

			Eigen::SparseMatrix <double> Conversion_Mat_2 = Conversion_Mat_1 * Constraint;

			// Solve the linear system
			Eigen::SparseLU <Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering<int>> solver;
			solver.compute(Conversion_Mat_2.transpose() * Power_network_inform.points.covariance * Conversion_Mat_2);
			Eigen::VectorXd lambda = solver.solve(mu_mean - Constraint.transpose() * mu + Conversion_Mat_2.transpose() * x);
			Eigen::VectorXd dx = Power_network_inform.points.covariance * Conversion_Mat_2 * lambda - x;
			x += alpha_iteration * dx;

			// Update the solution
			for(int item = 0; item < point_num; ++ item){
				dmu(item) = mu_scale(item) * pow(-log(1. - cdf(norm_dist, x(item))), .5) - mu(item);
				mu(item) += dmu(item);
			}

			Eigen::VectorXd mu_inv = pow(mu.array(), -1.);
			dmu *= mu_inv;

			count += 1;
		}

		// Output normalized mean demand field
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
		fout_name = "csv/processed/spatial_field/nominal_mean_demand_field_10km_ts_" + digit_zeros + std::to_string(tick) + ".csv";
		basic::write_file(mu, fout_name, col_name);
		digit_zeros.clear();
	}
}

void spatial_field::nominal_demand_inference(power_network::network_inform &Power_network_inform){
	int bz_num = Power_network_inform.points.bidding_zone.maxCoeff() + 1;
	int point_num = Power_network_inform.points.bidding_zone.size();

	// Read electricity demand data
	std::string fin_demand = "csv/input/spatial_field/demand_actual_2021.csv";
	auto fin_demand_dim = basic::get_file_dim(fin_demand);
	auto Demand_ts = basic::read_file(fin_demand_dim[0], fin_demand_dim[1], fin_demand);

	// Set sparse matrix for equality constaints
	Eigen::SparseMatrix <double> Constraint (point_num, bz_num);
	std::vector<Eigen::TripletXd> Constraint_Trip;
	Constraint_Trip.reserve(point_num);

	for(int point_iter = 0; point_iter < point_num; ++ point_iter){
		int bz_ID = Power_network_inform.points.bidding_zone(point_iter);
		double value = Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area;
		Constraint_Trip.push_back(Eigen::TripletXd(point_iter, bz_ID, value));
	}

	Constraint.setFromTriplets(Constraint_Trip.begin(), Constraint_Trip.end());

	// Spatial field inference
	BME(Power_network_inform, Constraint, Demand_ts);
}

void spatial_field::imbalance_inference(power_network::network_inform &Power_network_inform){
	int bz_num = Power_network_inform.points.bidding_zone.maxCoeff() + 1;
	int point_num = Power_network_inform.points.bidding_zone.size();

	// Read actual imbalance data
	std::string fin_imbalance = "csv/input/power_market/control_reserve_activated_2021.csv";
	auto fin_imbalance_dim = basic::get_file_dim(fin_imbalance);
	auto Imbalnce_ts = basic::read_file(fin_imbalance_dim[0], fin_imbalance_dim[1], fin_imbalance);
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
