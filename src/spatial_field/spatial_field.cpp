// Main source file for inference of spatial fields
#include "src/spatial_field/spatial_field.h"

namespace{
	void BME_copula(spatial_field::inference_inform &inform, power_network::network_inform &Power_network_inform, Eigen::SparseMatrix <double> &Constraint, double tol){
		int bz_num = Constraint.cols();
		int point_num = Constraint.rows();

		// Iterations
		int count = 0;
		Eigen::VectorXd dmu = Eigen::VectorXd::Constant(point_num, 1.);
		while(count < 5000 && dmu.lpNorm<Eigen::Infinity>() > tol){
			// Matrices for the linear system
			Eigen::SparseMatrix <double> Conversion_Mat_1(point_num, point_num);
			std::vector <Eigen::TripletXd> Conversion_Mat_1_Trip;
			Conversion_Mat_1.reserve(point_num);

			for(int item = 0; item < point_num; ++ item){
				double diff_x = .5 * inform.mu_scale(item) * pow((-log(1. - cdf(inform.norm_dist, inform.x(item)))), -.5) * pow(1. - cdf(inform.norm_dist, inform.x(item)), -1.) * pdf(inform.norm_dist, inform.x(item));
				Conversion_Mat_1_Trip.push_back(Eigen::TripletXd(item, item, diff_x));
			}
			Conversion_Mat_1.setFromTriplets(Conversion_Mat_1_Trip.begin(), Conversion_Mat_1_Trip.end());

			Eigen::SparseMatrix <double> Conversion_Mat_2 = Conversion_Mat_1 * Constraint;

			// Solve the linear system
			Eigen::MatrixXd solver = Conversion_Mat_2.transpose() * Power_network_inform.points.covariance * Conversion_Mat_2;
			Eigen::VectorXd lambda = solver.colPivHouseholderQr().solve(inform.mu_mean - Constraint.transpose() * inform.mu + Conversion_Mat_2.transpose() * inform.x);
			Eigen::VectorXd dx = Power_network_inform.points.covariance * Conversion_Mat_2 * lambda - inform.x;
			inform.x += inform.alpha_iteration * dx;

			// Update the solution
			for(int item = 0; item < point_num; ++ item){
				dmu(item) = inform.mu_scale(item) * pow(-log(1. - cdf(inform.norm_dist, inform.x(item))), .5) - inform.mu(item);
				inform.mu(item) += dmu(item);
			}
			Eigen::VectorXd mu_inv_0 = pow(inform.mu.array(), -1.);
			dmu *= mu_inv_0;

			std::cout << inform.mu.minCoeff() << "\t" << inform.mu.maxCoeff() << "\t" << dmu.lpNorm<Eigen::Infinity>() << "\n";
			count += 1;
		}
	}

	void BME_linear(spatial_field::inference_inform &inform, Eigen::SparseMatrix <double> &Constraint){
		// ---------------------------------------------------------------------
		// MaxEnt, known covariance and mean = 0, reference equiprobable space = x
		// COV^(-1) * x* - Constraint * lambda = 0
		// x* = COV * Constraint * lambda
		// Constraint^T * (x* + x_0) = Constraint^T * COV * Constraint * lambda + Constraint^T * x_0 = demand
		// lambda = solve(Constraint^T * COV * Constraint, demand - Constraint^T * x_0)
		// ---------------------------------------------------------------------
		// Solve the linear system
		Eigen::VectorXd Constraint_sum = Constraint.transpose() * inform.mu_scale;
		Eigen::SparseLU <Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering<int>> solver;
		solver.compute(inform.Conversion_Mat_1);
		Eigen::VectorXd lambda = solver.solve(inform.mu_mean - Constraint_sum);
		inform.mu = inform.Conversion_Mat_2 * lambda;
		inform.mu += inform.mu_scale;
		std::cout << inform.mu.minCoeff() << "\t" << inform.mu.maxCoeff() << "\n\n";
	}
}

void spatial_field::spatial_field_inference(power_network::network_inform &Power_network_inform){
	int bz_num = Power_network_inform.points.bidding_zone.maxCoeff() + 1;
	int point_num = Power_network_inform.points.bidding_zone.size();
	int Time = power_market::parameters::Time();
	double pi = boost::math::constants::pi<double>();

	// Read electricity demand data
	std::string fin_demand = "csv/input/spatial_field/demand_actual_2021.csv";
	auto fin_demand_dim = basic::get_file_dim(fin_demand);
	auto Demand_ts = basic::read_file(fin_demand_dim[0], fin_demand_dim[1], fin_demand);
	Demand_ts *= 1000; // MWh -> kWh

	// Read actual imbalance data
	bool fin_imbalance_row_name = 1;
	std::string fin_imbalance = "csv/input/power_market/control_reserve_activated_2021.csv";
	auto fin_imbalance_dim = basic::get_file_dim(fin_imbalance);
	auto Imbalance_ts_raw = basic::read_file(fin_imbalance_dim[0], fin_imbalance_dim[1], fin_imbalance);
	Eigen::MatrixXd Imbalance_ts(Time, bz_num);
	for(int zone_iter = 0; zone_iter < bz_num; ++ zone_iter){
		Imbalance_ts.col(zone_iter) = Imbalance_ts_raw.col(2 * zone_iter + fin_imbalance_row_name) - Imbalance_ts_raw.col(2 * zone_iter + fin_imbalance_row_name + 1);
	}

	// Read onshore wind power data
	bool fin_wind_on_row_name = 1;
	std::string fin_wind_on = "csv/input/power_market/generation_wind_onshore_forecast_2021.csv";
	auto fin_wind_on_dim = basic::get_file_dim(fin_wind_on);
	auto Wind_on_ts = basic::read_file(fin_wind_on_dim[0], fin_wind_on_dim[1], fin_wind_on);

//	// Set sparse matrix for equality constraints for nominal demand
//	Eigen::SparseMatrix <double> Constraint_demand (point_num, bz_num);
//	std::vector<Eigen::TripletXd> Constraint_demand_Trip;
//	Constraint_demand_Trip.reserve(point_num);
//	for(int point_iter = 0; point_iter < point_num; ++ point_iter){
//		int bz_ID = Power_network_inform.points.bidding_zone(point_iter);
//		double value = Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area;
//		Constraint_demand_Trip.push_back(Eigen::TripletXd(point_iter, bz_ID, value));
//	}
//	Constraint_demand.setFromTriplets(Constraint_demand_Trip.begin(), Constraint_demand_Trip.end());
//	//std::cout << Constraint_demand.topRows(10) << "\n\n";
//
//	// ------------------------------------------------------------------------------------------------------------------------------------------
//	// Infer the annual average of the normalized mean demand field
//	// ------------------------------------------------------------------------------------------------------------------------------------------
//	// Initialization of parameters
//	inference_inform nominal_demand;
//	nominal_demand.alpha_iteration = 1.;
//	nominal_demand.mu_mean = Demand_ts.colwise().sum().tail(bz_num);
//	nominal_demand.mu_mean /= Time;
//	// Mean of weibull distribution = mu_0_scale * gamma(1 + 1 / k), here k = 2
//	double mu_0_scale = nominal_demand.mu_mean.sum() / Constraint_demand.sum() * 2. / pow(pi, .5);
//	nominal_demand.mu_scale = Eigen::VectorXd::Constant(point_num, mu_0_scale);
//	nominal_demand.mu = Eigen::VectorXd::Constant(point_num, mu_0_scale);
//	nominal_demand.x = Eigen::VectorXd(point_num);
//	for(int item = 0; item < point_num; ++ item){
//		nominal_demand.x(item) = quantile(nominal_demand.norm_dist, 1. - exp(-pow(nominal_demand.mu(item) / nominal_demand.mu_scale(item), 2.)));
//	}
//
//	// Inference step
//	BME_copula(nominal_demand, Power_network_inform, Constraint_demand, 1E-12);
//
//	// Output the annual average of normalized mean demand field
//	std::string fout_name;
//	fout_name = "csv/processed/spatial_field/nominal_mean_demand_field_10km_annual_mean.csv";
//	std::vector <std::string> col_name;
//	col_name.push_back("nominal_mean_demand");
//	basic::write_file(nominal_demand.mu, fout_name, col_name);
//
//	// Set sparse matrix for equality constraints for imbalance
//	// must come after nominal demand!!
//	Eigen::SparseMatrix <double> Constraint_imbalance (point_num, bz_num);
//	std::vector<Eigen::TripletXd> Constraint_imbalance_Trip;
//	Constraint_imbalance_Trip.reserve(point_num);
//	for(int point_iter = 0; point_iter < point_num; ++ point_iter){
//		int bz_ID = Power_network_inform.points.bidding_zone(point_iter);
//		double value = nominal_demand.mu(point_iter) * Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
//		Constraint_imbalance_Trip.push_back(Eigen::TripletXd(point_iter, bz_ID, value));
//	}
//	Constraint_imbalance.setFromTriplets(Constraint_imbalance_Trip.begin(), Constraint_imbalance_Trip.end());
//
//	// ------------------------------------------------------------------------------------------------------------------------------------------
//	// Infer the annual average of the imbalance field
//	// ------------------------------------------------------------------------------------------------------------------------------------------
//	// Initialization of parameters
//	inference_inform imbalance;
//	imbalance.mu_mean = Imbalance_ts.colwise().sum().segment(fin_imbalance_row_name, bz_num);;
//	imbalance.mu_mean /= Time;
//	mu_0_scale = imbalance.mu_mean.sum() / Constraint_imbalance.sum();
//	imbalance.mu_scale = Eigen::VectorXd::Constant(point_num, mu_0_scale);
//	imbalance.Conversion_Mat_1 = Constraint_imbalance.transpose() * Power_network_inform.points.covariance * Constraint_imbalance;
//	imbalance.Conversion_Mat_2 = Power_network_inform.points.covariance * Constraint_imbalance;
//
//	// Inference step
//	BME_linear(imbalance, Constraint_imbalance);
//	std::cout << imbalance.mu_mean.transpose() << "\n\n";
//	std::cout << imbalance.mu.transpose() * Constraint_imbalance << "\n\n";
//
//	// Output the annual average of normalized mean demand field
//	fout_name = "csv/processed/spatial_field/imbalance_field_10km_annual_mean.csv";
//	col_name.clear();
//	col_name.push_back("imbalance");
//	basic::write_file(imbalance.mu, fout_name, col_name);
//
//	// ------------------------------------------------------------------------------------------------------------------------------------------
//	// Infer the normalized mean demand / imbalance field
//	// ------------------------------------------------------------------------------------------------------------------------------------------
//	// Initialization of parameters
//	nominal_demand.alpha_iteration = .01;
//	nominal_demand.mu_scale = nominal_demand.mu * 2. / pow(pi, .5);
//	nominal_demand.x_scale = Eigen::VectorXd(point_num);
//	for(int item = 0; item < point_num; ++ item){
//		nominal_demand.x_scale(item) = quantile(nominal_demand.norm_dist, 1. - exp(-pow(1., 2.)));
//	}
//	nominal_demand.mu = nominal_demand.mu_scale;
//	nominal_demand.x = nominal_demand.x_scale;
//	imbalance.mu_scale = imbalance.mu;
//
//	// Run the algorithm for each time slice
//	for(int tick = 0; tick < 1; ++ tick){
//		nominal_demand.mu_mean = Demand_ts.row(tick).tail(bz_num);
//
//		// Inference step
//		BME_copula(nominal_demand, Power_network_inform, Constraint_demand, 1E-3);
//		//std::cout << nominal_demand.mu.transpose() * Constraint_demand << "\n\n";
//
//		// Output normalized mean demand field
//		int count_zeros = 0;
//		int tick_temp = tick;
//		std::string digit_zeros;
//		while(int (tick_temp / 10) != 0){
//			count_zeros += 1;
//			tick_temp /= 10;
//		}
//		for(int item = 0; item < 5 - count_zeros; ++item){
//			digit_zeros += std::to_string(0);
//		}
//		col_name.clear();
//		col_name.push_back("nominal_mean_demand");
//		fout_name = "csv/processed/spatial_field/nominal_mean_demand_field_10km_ts_" + digit_zeros + std::to_string(tick) + ".csv";
//		basic::write_file(nominal_demand.mu, fout_name, col_name);
//
//		// Set sparse matrix for equality constraints for imbalance
//		// must come after nominal demand!!
//		std::vector<Eigen::TripletXd> Constraint_imbalance_Trip;
//		Constraint_imbalance_Trip.reserve(point_num);
//		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
//			int bz_ID = Power_network_inform.points.bidding_zone(point_iter);
//			double value = nominal_demand.mu(point_iter) * Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
//			Constraint_imbalance_Trip.push_back(Eigen::TripletXd(point_iter, bz_ID, value));
//		}
//		Constraint_imbalance.setFromTriplets(Constraint_imbalance_Trip.begin(), Constraint_imbalance_Trip.end());
//
//		imbalance.mu_mean = Imbalance_ts.row(tick).tail(bz_num);
//
//		// Inference step
//		BME_linear(imbalance, Constraint_imbalance);
//		std::cout << imbalance.mu.transpose() * Constraint_imbalance << "\n\n";
//
//		// Output the annual average of normalized mean demand field
//		fout_name = "csv/processed/spatial_field/imbalance_field_10km_ts_" + digit_zeros + std::to_string(tick) + ".csv";
//		col_name.clear();
//		col_name.push_back("imbalance");
//		basic::write_file(imbalance.mu, fout_name, col_name);
//		digit_zeros.clear();
//	}

	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Capacity factor fields
	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Set sparse matrix for equality constraints for nominal demand
	Eigen::MatrixXd Constraint_wind_on_dense = Eigen::MatrixXd::Zero(point_num, bz_num);
	for(int wind_iter = 0; wind_iter < Power_network_inform.plants.wind.node.size(); ++ wind_iter){
		int x_ID = int((Power_network_inform.plants.wind.x(wind_iter) - Power_network_inform.points.x.minCoeff()) / Power_network_inform.points.grid_length + .5);
		int y_ID = int((Power_network_inform.plants.wind.y(wind_iter) - Power_network_inform.points.y.minCoeff()) / Power_network_inform.points.grid_length + .5);
		int point_ID = Power_network_inform.points.coordinate_grid(x_ID, y_ID);
		if(point_ID == -1){
			continue;
		}
		int bz_ID = Power_network_inform.points.bidding_zone(point_ID);
		Constraint_wind_on_dense(point_ID, bz_ID) += Power_network_inform.plants.wind.cap(wind_iter);
	}
	Eigen::SparseMatrix <double> Constraint_wind_on = Constraint_wind_on_dense.sparseView(1E-12);

	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Infer the annual average of onshore wind capacity factor field
	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Initialization of parameters
	inference_inform wind_on_cf;
	wind_on_cf.alpha_iteration = 1.;
	wind_on_cf.mu_mean = Wind_on_ts.colwise().sum().segment(fin_wind_on_row_name, bz_num);
	wind_on_cf.mu_mean /= Time;
	// Mean of weibull distribution = mu_0_scale * gamma(1 + 1 / k), here k = 2
	double mu_0_scale = wind_on_cf.mu_mean.sum() / Constraint_wind_on.sum() * 2. / pow(pi, .5);
	wind_on_cf.mu_scale = Eigen::VectorXd::Constant(point_num, mu_0_scale);
	wind_on_cf.mu = Eigen::VectorXd::Constant(point_num, mu_0_scale);
	wind_on_cf.x = Eigen::VectorXd(point_num);
	for(int item = 0; item < point_num; ++ item){
		wind_on_cf.x(item) = quantile(wind_on_cf.norm_dist, 1. - exp(-pow(wind_on_cf.mu(item) / wind_on_cf.mu_scale(item), 2.)));
	}
	wind_on_cf.Conversion_Mat_1 = Constraint_wind_on.transpose() * Power_network_inform.points.covariance * Constraint_wind_on;
	wind_on_cf.Conversion_Mat_2 = Power_network_inform.points.covariance * Constraint_wind_on;

	// Inference step
	//BME_copula(wind_on_cf, Power_network_inform, Constraint_wind_on, 1E-12);
	BME_linear(wind_on_cf, Constraint_wind_on);
	//std::cout << wind_on_cf.mu.transpose() << "\n\n";
}

// Function that stores processed mean demand field
void spatial_field::spatial_field_store(power_network::network_inform &Power_network_inform, std::string fin_demand, std::string fin_imbalance, int Time){
	int row_num = Power_network_inform.points.bidding_zone.rows();
	Power_network_inform.points.nominal_mean_demand_field = Eigen::MatrixXd(row_num, Time);
	Power_network_inform.points.imbalance_field = Eigen::MatrixXd(row_num, Time);

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
		std::string fin_imbalance_temp = fin_imbalance + digit_zeros + std::to_string(tick) + ".csv";

		Power_network_inform.points.nominal_mean_demand_field.col(tick) = basic::read_file(row_num, 1, fin_demand_temp);
		Power_network_inform.points.imbalance_field.col(tick) = basic::read_file(row_num, 1, fin_imbalance_temp);
	}
}
