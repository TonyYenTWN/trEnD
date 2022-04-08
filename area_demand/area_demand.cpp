// Main Source File
#include <iostream>
#include <chrono>
#include <omp.h>
#include <boost/math/distributions/normal.hpp>
#include "Basic_Definitions.h"
#include "Geostat.cpp"
#include "rw_csv.cpp"

boost::math::normal norm_dist(0.0, 1.0);

int main(){
	auto start = chrono::high_resolution_clock::now();
	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast <chrono::microseconds> (stop - start);
	
	string fin_name;
	// Read population density data; can change to other types of spatial point data in the future
	fin_name = "input/population_density_10km.csv";
	int num_row = 4041; 
	int num_col = 4;
	MatrixXd sample_inform = read_file(num_row, num_col, fin_name);

	// Read population density data; can change to other types of spatial point data in the future
	fin_name = "input/demand_actual_2021.csv";
	int num_time = 8760; 
	int num_category = sample_inform.col(0).maxCoeff() + 1;
	MatrixXd Demand_ts = read_file(num_time, num_category + 1, fin_name);
	
	// Compute distance of the spatial points, and the corresponding covariance matrix of the random variable
	MatrixXd Distance_Points(num_row, num_row);
	MatrixXd Covariance_Points(num_row, num_row);
	
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
	// int num_category = sample_inform.col(0).maxCoeff() + 1;
	MatrixXd Constraint = MatrixXd::Zero(num_row, num_category);

	#pragma omp parallel
	{
		#pragma omp for
		for(int row_ID = 0; row_ID < num_row; ++ row_ID){
			Constraint(row_ID, int(sample_inform(row_ID, 0))) = sample_inform(row_ID, 1);
		}
	}

	VectorXd Constraint_sum = Constraint.colwise().sum();
	
	// Initialization of Demand Data
	VectorXd Demand_0(num_category);
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
	VectorXd mu_0 = VectorXd::LinSpaced(num_row, mu_0_mean, mu_0_mean); 
	VectorXd x_0(num_row);
	#pragma omp parallel
	{
		#pragma omp for
		for(int item = 0; item < num_row; ++ item){
			x_0(item) = quantile(norm_dist, 1 - exp(-pow(mu_0(item) / mu_0_mean, 2)));
		}
	}
	VectorXd dx_0 = VectorXd::LinSpaced(num_row, 1, 1);
	VectorXd dmu_0 = VectorXd::LinSpaced(num_row, 1, 1);
	VectorXd mu_inv_0(num_row);
	VectorXd lambda(num_category);
	VectorXd Conversion_vec(num_row);
	MatrixXd Conversion_Mat_1(num_row, num_row);
	MatrixXd Conversion_Mat_2(num_row, num_category);

	// Iterations
	count = 0; 
	while(count < 5000 && dmu_0.lpNorm<Infinity>() > pow(10, -12)){
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
	
	cout << "--------------------------------------------------------------------------" << endl;
	cout << "Annual Mean" << endl;
	cout << "Iterations: " << count << endl;
	cout << "--------------------------------------------------------------------------" << endl;
	cout << "Aggregated Demand at Each Bidding Zone:" << endl;
	cout << Constraint.transpose() * mu_0 << endl;
	cout << "\n" << endl;
	cout << "Range of Nominal Demand:" << endl;
	cout << mu_0.minCoeff() << endl;
	cout << mu_0.maxCoeff() << endl;
	cout << "\n" << endl;
	cout << "Relative Error:" << endl;
	cout << dmu_0.lpNorm<Infinity>() << endl;		
	cout << "\n" << endl;
	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast <chrono::microseconds> (stop - start);
	cout << "Timelapse of Code:\n";
	cout << duration.count() / pow(10, 6) << " seconds" << endl;
	cout << "\n";
	
	// ---------------------------------------------------------------------
	// Output Annual Average Normalized Mean Demand Field Data
	string fout_name;
	fout_name = "output/nominal_mean_demand_field_10km_annual_mean.csv";
	vector<string> col_name;
	col_name.push_back("nominal_mean_demand");
	write_file(mu_0, fout_name, col_name);
	
	// ---------------------------------------------------------------------
	// Infer the normalized mean demand field
	// ---------------------------------------------------------------------
	// Initialization of Constants, Vectors and Matrixs
	int tick_temp;
	int count_zeros;
	alpha_iteration = .01;
	string digit_zeros;
	VectorXd Demand(num_category);
	VectorXd mu_scale = mu_0 * 2 / pow(pi, .5);
	VectorXd x_scale(num_row);
	#pragma omp parallel
	{
		#pragma omp for
		for(int item = 0; item < num_row; ++ item){
			x_scale(item) = quantile(norm_dist, 1 - exp(-pow(1, 2)));
		}
	}
	VectorXd mu = mu_scale; 
	VectorXd x = x_scale;
	VectorXd dx = VectorXd::LinSpaced(num_row, 1, 1);
	VectorXd dmu = VectorXd::LinSpaced(num_row, 1, 1);
	VectorXd mu_inv(num_row);

	for(int tick = 0; tick < num_time; ++ tick){
		Demand = Demand_ts.row(tick).tail(num_category);
		
		// Iterated Optimization
		count = 0;
		dmu = VectorXd::LinSpaced(num_row, 1, 1);
		while(count < 5000 && dmu.lpNorm<Infinity>() > pow(10, -3)){
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
			digit_zeros += to_string(0);
		}
		fout_name = "output/nominal_mean_demand_field_10km_ts_" + digit_zeros + to_string(tick) + ".csv";
		write_file(mu, fout_name, col_name);
		digit_zeros.clear();
		
		if((tick + 1) % 24 == 0){
			cout << "--------------------------------------------------------------------------" << endl;
			cout << "Tick: " << tick << endl;
			cout << "Iterations: " << count << endl;
			cout << "--------------------------------------------------------------------------" << endl;
			cout << "Range of Nominal Demand:" << endl;
			cout << mu.minCoeff() << endl;
			cout << mu.maxCoeff() << endl;
			cout << "\n" << endl;
			cout << "Relative Error:" << endl;
			cout << dmu.lpNorm<Infinity>() << endl;		
			cout << "\n" << endl;
			stop = chrono::high_resolution_clock::now();
			duration = chrono::duration_cast <chrono::microseconds> (stop - start);
			cout << "Timelapse of Code:\n";
			cout << duration.count() / pow(10, 6) << " seconds" << endl;
			cout << "\n";
		}
	}
}   