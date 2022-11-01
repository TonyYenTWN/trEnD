// Main source file for estimation of spatial fields
#include "src/spatial_field/spatial_field.h"

namespace{
	void BME_copula(spatial_field::estimation_inform &inform, power_network::network_inform &Power_network_inform, Eigen::SparseMatrix <double> &Constraint, double tol){
		int bz_num = Constraint.cols();
		int point_num = Constraint.rows();

		// Iterations
		int count = 0;
		Eigen::VectorXd dmu = Eigen::VectorXd::Constant(point_num, 1.);
		//std::cout << inform.mu.minCoeff() << "\t" << inform.mu.maxCoeff() << "\t" << dmu.lpNorm<Eigen::Infinity>() << "\n";
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

			//std::cout << inform.mu.minCoeff() << "\t" << inform.mu.maxCoeff() << "\t" << dmu.lpNorm<Eigen::Infinity>() << "\n";
			count += 1;
		}
		//std::cout << inform.mu.minCoeff() << "\t" << inform.mu.maxCoeff() << "\t" << dmu.lpNorm<Eigen::Infinity>() << "\n\n";
	}

	void BME_linear(spatial_field::estimation_inform &inform, Eigen::SparseMatrix <double> &Constraint){
		// ---------------------------------------------------------------------
		// MaxEnt, known covariance and mean = 0, reference equiprobable space = x
		// COV^(-1) * x* - Constraint * lambda = 0
		// x* = COV * Constraint * lambda
		// Constraint^T * (x* + x_0) = Constraint^T * COV * Constraint * lambda + Constraint^T * x_0 = demand
		// lambda = solve(Constraint^T * COV * Constraint, demand - Constraint^T * x_0)
		// ---------------------------------------------------------------------
		// Solve the linear system
		Eigen::VectorXd Constraint_sum = Constraint.transpose() * inform.mu_scale;
		//std::cout << Constraint_sum.transpose() << "\n\n";
		Eigen::SparseLU <Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering<int>> solver;
		solver.compute(inform.Conversion_Mat_1);
		Eigen::VectorXd lambda = solver.solve(inform.mu_mean - Constraint_sum);
		inform.mu = inform.Conversion_Mat_2 * lambda;
		inform.mu += inform.mu_scale;
		//std::cout << inform.mu.minCoeff() << "\t" << inform.mu.maxCoeff() << "\n\n";
	}
}

void spatial_field::demand_imbalance_estimation(power_network::network_inform &Power_network_inform, power_market::market_inform &International_Market, configuration::process_config &process_par){
	int bz_num = Power_network_inform.points.bidding_zone.maxCoeff() + 1;
	int point_num = Power_network_inform.points.bidding_zone.size();
	int Time = configuration::parameters::Time();
	double pi = boost::math::constants::pi<double>();

	// Read electricity demand data
	Eigen::MatrixXd Demand_ts = International_Market.demand_default.leftCols(bz_num);
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

	// Set sparse matrix for equality constraints for nominal demand
	Eigen::SparseMatrix <double> Constraint_demand (point_num, bz_num);
	std::vector<Eigen::TripletXd> Constraint_demand_Trip;
	Constraint_demand_Trip.reserve(point_num);
	for(int point_iter = 0; point_iter < point_num; ++ point_iter){
		int bz_ID = Power_network_inform.points.bidding_zone(point_iter);
		double value = Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area;
		Constraint_demand_Trip.push_back(Eigen::TripletXd(point_iter, bz_ID, value));
	}
	Constraint_demand.setFromTriplets(Constraint_demand_Trip.begin(), Constraint_demand_Trip.end());
	//std::cout << Constraint_demand.topRows(10) << "\n\n";

	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Estimate the annual average of the normalized mean demand field
	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Initialization of parameters
	estimation_inform nominal_demand;
	nominal_demand.alpha_iteration = 1.;
	nominal_demand.mu_mean = Demand_ts.colwise().sum();
	nominal_demand.mu_mean /= Time;
	// Mean of weibull distribution = mu_0_scale * gamma(1 + 1 / k), here k = 2
	double mu_0_scale = nominal_demand.mu_mean.sum() / Constraint_demand.sum() * 2. / pow(pi, .5);
	nominal_demand.mu_scale = Eigen::VectorXd::Constant(point_num, mu_0_scale);
	nominal_demand.mu = Eigen::VectorXd::Constant(point_num, mu_0_scale);
	nominal_demand.x = Eigen::VectorXd(point_num);
	for(int item = 0; item < point_num; ++ item){
		nominal_demand.x(item) = quantile(nominal_demand.norm_dist, 1. - exp(-pow(nominal_demand.mu(item) / nominal_demand.mu_scale(item), 2.)));
	}

	// estimation step
	BME_copula(nominal_demand, Power_network_inform, Constraint_demand, 1E-12);
	std::cout << tick << ":\t" << nominal_demand.mu.transpose() * Constraint_demand << "\n\n";

	// Output the annual average of normalized mean demand field
	std::string fout_name;
	fout_name = "csv/processed/spatial_field/nominal_mean_demand_field_10km_annual_mean.csv";
	std::vector <std::string> col_name;
	col_name.push_back("nominal_mean_demand");
	basic::write_file(nominal_demand.mu, fout_name, col_name);

	// Set sparse matrix for equality constraints for imbalance
	// must come after nominal demand!!
	Eigen::SparseMatrix <double> Constraint_imbalance (point_num, bz_num);
	std::vector<Eigen::TripletXd> Constraint_imbalance_Trip;
	Constraint_imbalance_Trip.reserve(point_num);
	for(int point_iter = 0; point_iter < point_num; ++ point_iter){
		int bz_ID = Power_network_inform.points.bidding_zone(point_iter);
		double value = nominal_demand.mu(point_iter) * Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
		Constraint_imbalance_Trip.push_back(Eigen::TripletXd(point_iter, bz_ID, value));
	}
	Constraint_imbalance.setFromTriplets(Constraint_imbalance_Trip.begin(), Constraint_imbalance_Trip.end());

	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Estimate the annual average of the imbalance field
	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Initialization of parameters
	estimation_inform imbalance;
	imbalance.mu_mean = Imbalance_ts.colwise().sum();
	imbalance.mu_mean /= Time;
	mu_0_scale = imbalance.mu_mean.sum() / Constraint_imbalance.sum();
	imbalance.mu_scale = Eigen::VectorXd::Constant(point_num, mu_0_scale);
	//std::cout << imbalance.mu_scale.transpose() * Constraint_imbalance << "\n\n";
	imbalance.Conversion_Mat_1 = Constraint_imbalance.transpose() * Power_network_inform.points.covariance * Constraint_imbalance;
	imbalance.Conversion_Mat_2 = Power_network_inform.points.covariance * Constraint_imbalance;

	// estimation step
	BME_linear(imbalance, Constraint_imbalance);
	//std::cout << imbalance.mu_mean.transpose() << "\n";
	std::cout << imbalance.mu.transpose() * Constraint_imbalance << "\n\n";

	// Output the annual average of normalized mean demand field
	fout_name = "csv/processed/spatial_field/imbalance_field_10km_annual_mean.csv";
	col_name.clear();
	col_name.push_back("imbalance");
	basic::write_file(imbalance.mu, fout_name, col_name);

	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Estimate the normalized mean demand / imbalance field
	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Initialization of parameters
	nominal_demand.alpha_iteration = 1.;
	nominal_demand.mu_scale = nominal_demand.mu * 2. / pow(pi, .5);
	nominal_demand.x_scale = Eigen::VectorXd(point_num);
	for(int item = 0; item < point_num; ++ item){
		nominal_demand.x_scale(item) = quantile(nominal_demand.norm_dist, 1. - exp(-pow(1., 2.)));
	}
	nominal_demand.mu = nominal_demand.mu_scale;
	nominal_demand.x = nominal_demand.x_scale;
	imbalance.mu_scale = imbalance.mu;

	// Run the algorithm for each time slice
	for(int tick = process_par.time_boundary[0]; tick < process_par.time_boundary[0] + process_par.time_boundary[1]; ++ tick){
		nominal_demand.mu_mean = Demand_ts.row(tick);

		// Estimation step
		BME_copula(nominal_demand, Power_network_inform, Constraint_demand, 1E-3);
		std::cout << nominal_demand.mu.transpose() * Constraint_demand << "\n\n";

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
		col_name.clear();
		col_name.push_back("nominal_mean_demand");
		fout_name = "csv/processed/spatial_field/nominal_mean_demand_field_10km_ts_" + digit_zeros + std::to_string(tick) + ".csv";
		basic::write_file(nominal_demand.mu, fout_name, col_name);

		// Set sparse matrix for equality constraints for imbalance
		// must come after nominal demand!!
		std::vector<Eigen::TripletXd> Constraint_imbalance_temp_Trip;
		Constraint_imbalance_temp_Trip.reserve(point_num);
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int bz_ID = Power_network_inform.points.bidding_zone(point_iter);
			double value = nominal_demand.mu(point_iter) * Power_network_inform.points.population_density(point_iter) * Power_network_inform.points.point_area / 1000.;
			Constraint_imbalance_temp_Trip.push_back(Eigen::TripletXd(point_iter, bz_ID, value));
		}
		Constraint_imbalance.setFromTriplets(Constraint_imbalance_temp_Trip.begin(), Constraint_imbalance_temp_Trip.end());
		imbalance.Conversion_Mat_1 = Constraint_imbalance.transpose() * Power_network_inform.points.covariance * Constraint_imbalance;
		imbalance.Conversion_Mat_2 = Power_network_inform.points.covariance * Constraint_imbalance;
		//std::cout << imbalance.mu_scale.transpose() * Constraint_imbalance << "\n\n";

		imbalance.mu_mean = Imbalance_ts.row(tick);

		// Estimation step
		BME_linear(imbalance, Constraint_imbalance);
		//std::cout << imbalance.mu_mean.transpose() << "\n";
		std::cout << imbalance.mu.transpose() * Constraint_imbalance << "\n\n";

		// Output the annual average of imbalance field
		fout_name = "csv/processed/spatial_field/imbalance_field_10km_ts_" + digit_zeros + std::to_string(tick) + ".csv";
		col_name.clear();
		col_name.push_back("imbalance");
		basic::write_file(imbalance.mu, fout_name, col_name);
		digit_zeros.clear();
	}
}

// Function that calculates onshore wind capacity factor field
void spatial_field::wind_on_cf_estimation(power_network::network_inform &Power_network_inform, configuration::process_config &process_par){
	int bz_num = Power_network_inform.points.bidding_zone.maxCoeff() + 1;
	int point_num = Power_network_inform.points.bidding_zone.size();
	int Time = configuration::parameters::Time();
	double pi = boost::math::constants::pi<double>();

	// Read onshore wind power data
	bool fin_wind_on_row_name = 1;
	std::string fin_wind_on = "csv/input/power_market/generation_wind_onshore_forecast_2021.csv";
	auto fin_wind_on_dim = basic::get_file_dim(fin_wind_on);
	auto Wind_on_ts = basic::read_file(fin_wind_on_dim[0], fin_wind_on_dim[1], fin_wind_on);

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
	// Estimate the annual average of onshore wind capacity factor field
	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Initialization of parameters
	estimation_inform wind_on_cf;
	wind_on_cf.alpha_iteration = 1.;
	wind_on_cf.mu_mean = (Wind_on_ts.colwise().sum()).segment(fin_wind_on_row_name, bz_num);
	wind_on_cf.mu_mean /= Time;

	// Column number 4 is redundant
	int redunant_num = 1;
	std::map<int, int> corres_col;
	corres_col.insert(std::pair<int, int>(0, 0));
	corres_col.insert(std::pair<int, int>(1, 1));
	corres_col.insert(std::pair<int, int>(2, 2));
	corres_col.insert(std::pair<int, int>(3, 3));
	corres_col.insert(std::pair<int, int>(4, -1));
	Eigen::SparseMatrix <double> Redundant_col(bz_num - redunant_num, bz_num);
	std::vector <Eigen::TripletXd> Redundant_col_trip;
	Redundant_col_trip.reserve(bz_num - redunant_num);
	for(int zone_iter = 0; zone_iter < bz_num; ++ zone_iter){
		if(corres_col[zone_iter] == -1){
			continue;
		}
		else{
			Redundant_col_trip.push_back(Eigen::TripletXd(corres_col[zone_iter], zone_iter, 1.));
		}
	}
	Redundant_col.setFromTriplets(Redundant_col_trip.begin(), Redundant_col_trip.end());
	wind_on_cf.mu_mean = Redundant_col * wind_on_cf.mu_mean;
	Constraint_wind_on = Constraint_wind_on * Redundant_col.transpose();

	// Mean of weibull distribution = mu_0_scale * gamma(1 + 1 / k), here k = 2
	double mu_0_scale = wind_on_cf.mu_mean.sum() / Constraint_wind_on.sum() * 2. / pow(pi, .5);
	wind_on_cf.mu_scale = Eigen::VectorXd::Constant(point_num, mu_0_scale);
	wind_on_cf.mu = Eigen::VectorXd::Constant(point_num, mu_0_scale);
	wind_on_cf.x = Eigen::VectorXd(point_num);
	for(int item = 0; item < point_num; ++ item){
		wind_on_cf.x(item) = quantile(wind_on_cf.norm_dist, 1. - exp(-pow(wind_on_cf.mu(item) / wind_on_cf.mu_scale(item), 2.)));
	}

	// Estimation step
	BME_copula(wind_on_cf, Power_network_inform, Constraint_wind_on, 1E-12);
	std::cout << wind_on_cf.mu.transpose() * Constraint_wind_on << "\n\n";

	// Output the annual average of onshore wind capacity factor field
	std::string fout_name = "csv/processed/spatial_field/wind_onshore_cf_field_10km_annual_mean.csv";
	std::vector <std::string> col_name;
	col_name.push_back("wind_onshore_cf");
	basic::write_file(wind_on_cf.mu, fout_name, col_name);

	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Estimate onshore wind capacity factor field
	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Initialization of parameters
	wind_on_cf.alpha_iteration = .1;
	wind_on_cf.mu_scale = wind_on_cf.mu * 2. / pow(pi, .5);
	wind_on_cf.x_scale = Eigen::VectorXd(point_num);
	for(int item = 0; item < point_num; ++ item){
		wind_on_cf.x_scale(item) = quantile(wind_on_cf.norm_dist, 1. - exp(-pow(1., 2.)));
	}
	wind_on_cf.mu = wind_on_cf.mu_scale;
	wind_on_cf.x = wind_on_cf.x_scale;

	// Run the algorithm for each time slice
	for(int tick = process_par.time_boundary[0]; tick < process_par.time_boundary[0] + process_par.time_boundary[1]; ++ tick){
		wind_on_cf.mu_mean = (Wind_on_ts.row(tick)).segment(fin_wind_on_row_name, bz_num);
		wind_on_cf.mu_mean = Redundant_col * wind_on_cf.mu_mean;

		// estimation step
		BME_copula(wind_on_cf, Power_network_inform, Constraint_wind_on, 1E-3);
		std::cout << tick << ":\t" << wind_on_cf.mu.transpose() * Constraint_wind_on << "\n\n";

		// Output onshore wind capacity factor
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
		col_name.clear();
		col_name.push_back("wind_onshore_cf");
		fout_name = "csv/processed/spatial_field/wind_onshore_cf_field_10km_ts_" + digit_zeros + std::to_string(tick) + ".csv";
		basic::write_file(wind_on_cf.mu, fout_name, col_name);
	}
}

// Function that calculates solar radiation field
void spatial_field::solar_radiation_estimation(power_network::network_inform &Power_network_inform, configuration::process_config &process_par){
	int point_num = Power_network_inform.points.bidding_zone.size();
	int Time = configuration::parameters::Time();
	double pi = boost::math::constants::pi<double>();

	// Read solar radiation data
	std::string fin_solar = "csv/input/power_market/solar_radiation_2021.csv";
	auto fin_solar_dim = basic::get_file_dim(fin_solar);
	auto Solar_ts = basic::read_file(fin_solar_dim[0], fin_solar_dim[1], fin_solar);
	bool fin_solar_row_name = 1;
	int station_num = fin_solar_dim[1] - fin_solar_row_name;

	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Estimate the annual average of solar radiation field
	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Initialization of parameters
	Eigen::VectorXd solar_radiation_mean = Eigen::VectorXd::Zero(station_num);

	for(int station_iter = 0; station_iter < station_num; ++ station_iter){
		int col_ID = station_iter + fin_solar_row_name;
		int tick_count = 0;
		for(int tick = 0; tick < Time; ++ tick){
			if(Solar_ts(tick, col_ID) < 0.){
				continue;
			}

			solar_radiation_mean(station_iter) += Solar_ts(tick, col_ID);
			tick_count += 1;
		}
		solar_radiation_mean(station_iter) /= tick_count;
	}
	//std::cout << solar_radiation_mean.transpose() << "\n\n";

	// Set sparse matrix for equality constraints for solar radiation
	Eigen::VectorXi station_freq = Eigen::VectorXi::Zero(point_num);
	Eigen::VectorXd average_field = -Eigen::VectorXd::Ones(point_num);
	for(int station_iter = 0; station_iter < station_num; ++ station_iter){
		int point_ID = Power_network_inform.weather_stations.point(station_iter);
		if(point_ID == -1){
			continue;
		}
		average_field(point_ID) *= station_freq(point_ID);
		station_freq(point_ID) += 1;
		average_field(point_ID) += solar_radiation_mean(station_iter);
		average_field(point_ID) /= station_freq(point_ID);
	}

	estimation_inform solar_radiation;
	solar_radiation.alpha_iteration = 1.;
	solar_radiation.mu_mean = Eigen::VectorXd(station_num);
	std::vector<Eigen::TripletXd> Constraint_solar_Trip;
	Constraint_solar_Trip.reserve(station_num);
	int constraint_count = 0;
	for(int point_iter = 0; point_iter < point_num; ++ point_iter){
		if(average_field(point_iter) < 0.){
			continue;
		}
		Constraint_solar_Trip.push_back(Eigen::TripletXd(point_iter, constraint_count, 1.));
		solar_radiation.mu_mean(constraint_count) = average_field(point_iter);
		constraint_count += 1;
	}
	solar_radiation.mu_mean = solar_radiation.mu_mean.head(constraint_count);
	//std::cout << solar_radiation.mu_mean.transpose() << "\n\n";
	Eigen::SparseMatrix <double> Constraint_solar(point_num, Constraint_solar_Trip.size());
	Constraint_solar.setFromTriplets(Constraint_solar_Trip.begin(), Constraint_solar_Trip.end());

	// Mean of weibull distribution = mu_0_scale * gamma(1 + 1 / k), here k = 2
	double mu_0_scale = solar_radiation.mu_mean.sum() / Constraint_solar.sum() * 2. / pow(pi, .5);
	solar_radiation.mu_scale = Eigen::VectorXd::Constant(point_num, mu_0_scale);
	solar_radiation.mu = Eigen::VectorXd::Constant(point_num, mu_0_scale);
	solar_radiation.x = Eigen::VectorXd(point_num);
	for(int item = 0; item < point_num; ++ item){
		solar_radiation.x(item) = quantile(solar_radiation.norm_dist, 1. - exp(-pow(solar_radiation.mu(item) / solar_radiation.mu_scale(item), 2.)));
	}

	// Estimation step
	BME_copula(solar_radiation, Power_network_inform, Constraint_solar, 1E-12);
	std::cout << tick << ":\t" << solar_radiation.mu.transpose() * Constraint_solar << "\n\n";

	// Output the annual average of normalized mean demand field
	std::string fout_name;
	fout_name = "csv/processed/spatial_field/solar_radiation_field_10km_annual_mean.csv";
	std::vector <std::string> col_name;
	col_name.push_back("solar_radiation");
	basic::write_file(solar_radiation.mu, fout_name, col_name);

	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Estimate solar radiation field
	// ------------------------------------------------------------------------------------------------------------------------------------------
	// Initialization of parameters
	solar_radiation.alpha_iteration = .5;
	solar_radiation.mu_scale = solar_radiation.mu * 2. / pow(pi, .5);
	solar_radiation.x_scale = Eigen::VectorXd(point_num);
	for(int item = 0; item < point_num; ++ item){
		solar_radiation.x_scale(item) = quantile(solar_radiation.norm_dist, 1. - exp(-pow(1., 2.)));
	}
	solar_radiation.mu = solar_radiation.mu_scale;
	solar_radiation.x = solar_radiation.x_scale;

	// Run the algorithm for each time slice
	for(int tick = process_par.time_boundary[0]; tick < process_par.time_boundary[0] + process_par.time_boundary[1]; ++ tick){
		Eigen::VectorXd solar_radiation_temp = -Eigen::VectorXd::Ones(station_num);
		for(int station_iter = 0; station_iter < station_num; ++ station_iter){
			int col_ID = station_iter + fin_solar_row_name;
			if(Solar_ts(tick, col_ID) < 0.){
				continue;
			}
			solar_radiation_temp(station_iter) = Solar_ts(tick, col_ID);
		}
		//std::cout << solar_radiation_temp.transpose() << "\n\n";

		// Set sparse matrix for equality constraints for solar radiation
		Eigen::VectorXi station_freq_temp = Eigen::VectorXi::Zero(point_num);
		Eigen::VectorXd average_field_temp = -Eigen::VectorXd::Ones(point_num);
		for(int station_iter = 0; station_iter < station_num; ++ station_iter){
			if(solar_radiation_temp(station_iter) < 0.){
				continue;
			}

			int point_ID = Power_network_inform.weather_stations.point(station_iter);
			if(point_ID == -1){
				continue;
			}
			average_field_temp(point_ID) *= station_freq_temp(point_ID);
			station_freq_temp(point_ID) += 1;
			average_field_temp(point_ID) += solar_radiation_temp(station_iter);
			average_field_temp(point_ID) /= station_freq_temp(point_ID);
		}

		std::vector<Eigen::TripletXd> Constraint_solar_Trip_temp;
		Constraint_solar_Trip_temp.reserve(station_num);
		solar_radiation.mu_mean = Eigen::VectorXd(station_num);
		int constraint_count_temp = 0;
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			if(average_field_temp(point_iter) < 0.){
				continue;
			}
			average_field_temp(point_iter) += 1.;
			Constraint_solar_Trip_temp.push_back(Eigen::TripletXd(point_iter, constraint_count_temp, 1.));
			solar_radiation.mu_mean(constraint_count_temp) = average_field_temp(point_iter);
			constraint_count_temp += 1;
		}
		solar_radiation.mu_mean = solar_radiation.mu_mean.head(constraint_count_temp);
		//std::cout << solar_radiation.mu_mean.transpose() << "\n\n";
		Eigen::SparseMatrix <double> Constraint_solar_temp(point_num, Constraint_solar_Trip_temp.size());
		Constraint_solar_temp.setFromTriplets(Constraint_solar_Trip_temp.begin(), Constraint_solar_Trip_temp.end());

		// Estimation step
		BME_copula(solar_radiation, Power_network_inform, Constraint_solar_temp, 1E-3);
		solar_radiation.mu = solar_radiation.mu.array() - 1.;
		std::cout << solar_radiation.mu.transpose() * Constraint_solar_temp << "\n\n";

		// Output onshore wind capacity factor
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
		col_name.clear();
		col_name.push_back("solar_radiation");
		fout_name = "csv/processed/spatial_field/solar_radiation_field_10km_ts_" + digit_zeros + std::to_string(tick) + ".csv";
		basic::write_file(solar_radiation.mu, fout_name, col_name);
	}
}

// Function that stores processed mean demand field
void spatial_field::spatial_field_store(power_network::network_inform &Power_network_inform, fin_field fin_field_processed, configuration::process_config &process_par, int Time){
	int row_num = Power_network_inform.points.bidding_zone.rows();
	Power_network_inform.points.nominal_mean_demand_field = Eigen::MatrixXd(row_num, Time);
	Power_network_inform.points.imbalance_field = Eigen::MatrixXd(row_num, Time);
	Power_network_inform.points.wind_on_cf = Eigen::MatrixXd(row_num, Time);
	Power_network_inform.points.solar_cf = Eigen::MatrixXd(row_num, Time);

	//for(int tick = 0; tick < Time; ++ tick){
	for(int tick = process_par.time_boundary[0]; tick < process_par.time_boundary[0] + process_par.time_boundary[1]; ++ tick){
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
		std::string fin_demand_temp = fin_field_processed.demand + digit_zeros + std::to_string(tick) + ".csv";
		std::string fin_imbalance_temp = fin_field_processed.imbalance + digit_zeros + std::to_string(tick) + ".csv";
		std::string fin_wind_on_temp = fin_field_processed.wind_on + digit_zeros + std::to_string(tick) + ".csv";
		std::string fin_solar_temp = fin_field_processed.solar + digit_zeros + std::to_string(tick) + ".csv";

		Power_network_inform.points.nominal_mean_demand_field.col(tick) = basic::read_file(row_num, 1, fin_demand_temp);
		Power_network_inform.points.imbalance_field.col(tick) = basic::read_file(row_num, 1, fin_imbalance_temp);
		Power_network_inform.points.wind_on_cf.col(tick) = basic::read_file(row_num, 1, fin_wind_on_temp);
		Eigen::VectorXd solar_radiation = basic::read_file(row_num, 1, fin_solar_temp);

		for(int point_iter = 0; point_iter < row_num; ++ point_iter){
			Power_network_inform.points.solar_cf(point_iter, tick) = solar_cf_calculation(solar_radiation(point_iter));
		}
	}
}
