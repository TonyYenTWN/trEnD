// Main Source File that disaggregate electricity demand of bidding zones to demand density field at each spatial point
#include <iostream>
#include <chrono>
#include <boost/math/distributions/normal.hpp>
#include "Geostat.h"
#include "../basic/rw_csv.cpp"

boost::math::normal norm_dist(0.0, 1.0);

int main(){
	auto start = std::chrono::high_resolution_clock::now();
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast <std::chrono::microseconds> (stop - start);
	
	std::string fin_name;
	// Read population density data; can change to other types of spatial point data in the future
	fin_name = "input/population_density_10km.csv";
	int num_row = 4041; 
	int num_col = 6;
	Eigen::MatrixXd sample_inform = read_file(num_row, num_col, fin_name);

	// Read electricity demand data; can change to other types of spatial point data in the future
	fin_name = "input/demand_actual_2021.csv";
	int num_time = 8760; 
	int num_category = sample_inform.col(0).maxCoeff() + 1;
	Eigen::MatrixXd Demand_ts = read_file(num_time, num_category + 1, fin_name);
	
	// Compute distance of the spatial points, and the corresponding covariance matrix of the random variable
	Eigen::MatrixXd Distance_Points(num_row, num_row);
	Eigen::MatrixXd Covariance_Points(num_row, num_row);
	
	#pragma omp parallel
	{
		#pragma omp for
		for(int row_ID = 0; row_ID < num_row; ++ row_ID){
			for(int col_ID = row_ID; col_ID < num_row; ++ col_ID){
				if(row_ID == col_ID){
					Distance_Points(row_ID, col_ID) = 0;
				}
				else{
					Distance_Points(row_ID, col_ID) = geodist(sample_inform.row(row_ID).tail(2) * pi / 180, sample_inform.row(col_ID).tail(2) * pi / 180);
				}
				Covariance_Points(row_ID, col_ID) = exp(-pow(Distance_Points(row_ID, col_ID), 1) / 2 / pow(10, 4));
				Distance_Points(col_ID, row_ID) = Distance_Points(row_ID, col_ID);
				Covariance_Points(col_ID, row_ID) = Covariance_Points(row_ID, col_ID);
			}
		}
	}
	
	// Read Demand File of Each Demand Pool
	// Inititalization of Constraint Matrix
	Eigen::MatrixXd Constraint = Eigen::MatrixXd::Zero(num_row, num_category);

	#pragma omp parallel
	{
		#pragma omp for
		for(int row_ID = 0; row_ID < num_row; ++ row_ID){
			Constraint(row_ID, int(sample_inform(row_ID, 0))) = sample_inform(row_ID, 1);
		}
	}

	Eigen::VectorXd Constraint_sum = Constraint.colwise().sum();
	
	// Initialization of Demand Data
	Eigen::VectorXd Demand_0(num_category);
	Demand_0 = Demand_ts.colwise().sum().tail(num_category);
	Demand_0 /= num_time;
	double multiply_factor = 1;
	Demand_0 *= multiply_factor;
	double mu_0_mean = Demand_0.sum() / Constraint.sum() * 2 / pow(pi, .5);
	
	// ---------------------------------------------------------------------
	// Infer temporal average of the normalized mean demand field
	// ---------------------------------------------------------------------
	// MaxEnt, known covariance and mean = 0, reference equiprobable space = x
	// COV^(-1) * x* - Constraint * lambda = 0
	// x* = COV * Constraint * lambda
	// Constraint^T * x* = Constraint^T * COV * Constraint * lambda = demand
	// lambda = solve(Constraint^T * COV * Constraint, demand)
	// ---------------------------------------------------------------------
	/*
	// Uncomment this section to run the code for linear case
	MatrixXd Conversion_Mat_1 = Constraint.transpose() * Covariance_Points * Constraint;
	MatrixXd Conversion_Mat_2 = Covariance_Points * Constraint;
	VectorXd lambda = Conversion_Mat_1.colPivHouseholderQr().solve(Demand_0 - mu_0_mean * Constraint_sum);
	VectorXd mu_0 = Conversion_Mat_2 * lambda;
	mu_0.array() += mu_0_mean;
	cout << Constraint.transpose() * mu_0 << "\n" << endl;
	cout << mu_0.minCoeff() << endl;
	cout << mu_0.maxCoeff() << endl;
	cout << "\n" << endl;
	*/
	// ---------------------------------------------------------------------
	// Weibull Distribution using Gaussian Copula
	
	// Uncomment this section to run the code for non-linear case
	// Initialization of constants, vectors, and Matrixs
	int count;
	double alpha_iteration = 1;
	Eigen::VectorXd mu_0 = Eigen::VectorXd::LinSpaced(num_row, mu_0_mean, mu_0_mean); 
	Eigen::VectorXd x_0(num_row);
	#pragma omp parallel
	{
		#pragma omp for
		for(int item = 0; item < num_row; ++ item){
			x_0(item) = quantile(norm_dist, 1 - exp(-pow(mu_0(item) / mu_0_mean, 2)));
		}
	}
	Eigen::VectorXd dx_0 = Eigen::VectorXd::LinSpaced(num_row, 1, 1);
	Eigen::VectorXd dmu_0 = Eigen::VectorXd::LinSpaced(num_row, 1, 1);
	Eigen::VectorXd mu_inv_0(num_row);
	Eigen::VectorXd lambda(num_category);
	Eigen::VectorXd Conversion_vec(num_row);
	Eigen::MatrixXd Conversion_Mat_1(num_row, num_row);
	Eigen::MatrixXd Conversion_Mat_2(num_row, num_category);

	// Iterations
	count = 0; 
	while(count < 5000 && dmu_0.lpNorm<Eigen::Infinity>()> pow(10, -12)){
		#pragma omp parallel
		{
			#pragma omp for
			for(int item = 0; item < num_row; ++ item){
				Conversion_vec(item) = .5 * mu_0_mean * pow((-log(1 - cdf(norm_dist, x_0(item)))), -.5) * pow(1 - cdf(norm_dist, x_0(item)), -1) * pdf(norm_dist, x_0(item));
			}
		} 
		Conversion_Mat_1 = Conversion_vec.asDiagonal();
		Conversion_Mat_2 = Conversion_Mat_1 * Constraint;
		
		lambda = (Conversion_Mat_2.transpose() * Covariance_Points * Conversion_Mat_2).colPivHouseholderQr().solve(Demand_0 - Constraint.transpose() * mu_0 + Conversion_Mat_2.transpose() * x_0);
		dx_0 = Covariance_Points * Conversion_Mat_2 * lambda - x_0;
		x_0 += alpha_iteration * dx_0;

		#pragma omp parallel
		{
			#pragma omp for
			for(int item = 0; item < num_row; ++ item){
				dmu_0(item) = mu_0_mean * pow(-log(1 - cdf(norm_dist, x_0(item))), .5) - mu_0(item);
				mu_0(item) += dmu_0(item);
			}
		}
		mu_inv_0 = pow(mu_0.array(), -1);
		dmu_0 *= mu_inv_0;
		
		count += 1;
	}
	
	std::cout << "--------------------------------------------------------------------------" << std::endl;
	std::cout << "Annual Mean" << std::endl;
	std::cout << "Iterations: " << count << std::endl;
	std::cout << "--------------------------------------------------------------------------" << std::endl;
	std::cout << "Aggregated Demand at Each Bidding Zone:" << std::endl;
	std::cout << Constraint.transpose() * mu_0 << std::endl;
	std::cout << "\n" << std::endl;
	std::cout << "Range of Nominal Demand:" << std::endl;
	std::cout << mu_0.minCoeff() << std::endl;
	std::cout << mu_0.maxCoeff() << std::endl;
	std::cout << "\n" << std::endl;
	std::cout << "Relative Error:" << std::endl;
	std::cout << dmu_0.lpNorm<Eigen::Infinity>() << std::endl;		
	std::cout << "\n" << std::endl;
	stop = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast <std::chrono::microseconds> (stop - start);
	std::cout << "Timelapse of Code:\n";
	std::cout << duration.count() / pow(10, 6) << " seconds" << std::endl;
	std::cout << "\n";
	
	// ---------------------------------------------------------------------
	// Output Annual Average Normalized Mean Demand Field Data
	std::string fout_name;
	fout_name = "output/nominal_mean_demand_field_10km_annual_mean.csv";
	std::vector<std::string> col_name;
	col_name.push_back("nominal_mean_demand");
	write_file(mu_0, fout_name, col_name);
	
	// ---------------------------------------------------------------------
	// Infer the normalized mean demand field
	// ---------------------------------------------------------------------
	// Initialization of Constants, Vectors and Matrixs
	int tick_temp;
	int count_zeros;
	alpha_iteration = .01;
	std::string digit_zeros;
	Eigen::VectorXd Demand(num_category);
	Eigen::VectorXd mu_scale = mu_0 * 2 / pow(pi, .5);
	Eigen::VectorXd x_scale(num_row);
	#pragma omp parallel
	{
		#pragma omp for
		for(int item = 0; item < num_row; ++ item){
			x_scale(item) = quantile(norm_dist, 1 - exp(-pow(1, 2)));
		}
	}
	Eigen::VectorXd mu = mu_scale; 
	Eigen::VectorXd x = x_scale;
	Eigen::VectorXd dx = Eigen::VectorXd::LinSpaced(num_row, 1, 1);
	Eigen::VectorXd dmu = Eigen::VectorXd::LinSpaced(num_row, 1, 1);
	Eigen::VectorXd mu_inv(num_row);

	for(int tick = 0; tick < num_time; ++ tick){
		Demand = Demand_ts.row(tick).tail(num_category);
		
		// Iterated Optimization
		count = 0;
		dmu = Eigen::VectorXd::LinSpaced(num_row, 1, 1);
		while(count < 5000 && dmu.lpNorm<Eigen::Infinity>()> pow(10, -3)){
		#pragma omp parallel
		{
			#pragma omp for			
			for(int item = 0; item < num_row; ++ item){
				Conversion_vec(item) = .5 * mu_scale(item) * pow((-log(1 - cdf(norm_dist, x(item)))), -.5) * pow((-log(1 - cdf(norm_dist, x(item)))), -1) * pdf(norm_dist, x(item));
			}
		}	
			Conversion_Mat_1 = Conversion_vec.asDiagonal();
			Conversion_Mat_2 = Conversion_Mat_1 * Constraint;
			
			lambda = (Conversion_Mat_2.transpose() * Covariance_Points * Conversion_Mat_2).colPivHouseholderQr().solve(Demand - Constraint.transpose() * mu + Conversion_Mat_2.transpose() * x);
			dx = Covariance_Points * Conversion_Mat_2 * lambda - x;
			x += alpha_iteration * dx;

		#pragma omp parallel
		{
			#pragma omp for
			for(int item = 0; item < num_row; ++ item){
				dmu(item) = mu_scale(item) * pow(-log(1 - cdf(norm_dist, x(item))), .5) - mu(item);
				mu(item) += dmu(item);
			}
		}
			mu_inv = pow(mu.array(), -1);
			dmu *= mu_inv;
			
			count += 1;
		}
	
		// Output Normalized Mean Demand Field Data
		count_zeros = 0;
		tick_temp = tick;
		while(int (tick_temp / 10) != 0){
			count_zeros += 1;
			tick_temp /= 10;
		}
		for(int item = 0; item < 5 - count_zeros; ++item){
			digit_zeros += std::to_string(0);
		}
		fout_name = "output/nominal_mean_demand_field_10km_ts_" + digit_zeros + std::to_string(tick) + ".csv";
		write_file(mu, fout_name, col_name);
		digit_zeros.clear();
		
		if((tick + 1) % 24 == 0){
			std::cout << "--------------------------------------------------------------------------" << std::endl;
			std::cout << "Tick: " << tick << std::endl;
			std::cout << "Iterations: " << count << std::endl;
			std::cout << "--------------------------------------------------------------------------" << std::endl;
			std::cout << "Range of Nominal Demand:" << std::endl;
			std::cout << mu.minCoeff() << std::endl;
			std::cout << mu.maxCoeff() << std::endl;
			std::cout << "\n" << std::endl;
			std::cout << "Relative Error:" << std::endl;
			std::cout << dmu.lpNorm<Eigen::Infinity>() << std::endl;		
			std::cout << "\n" << std::endl;
			stop = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast <std::chrono::microseconds> (stop - start);
			std::cout << "Timelapse of Code:\n";
			std::cout << duration.count() / pow(10, 6) << " seconds" << std::endl;
			std::cout << "\n";
		}
	}
}   