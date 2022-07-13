// Main Source File
#include <iostream>
#include <chrono>
#include <omp.h>
#include <boost/math/distributions/normal.hpp>
#include "../basic/rw_csv.cpp"

boost::math::normal norm_dist(0.0, 1.0);

Eigen::VectorXd Neighbor_Search(Eigen::VectorXd par, Eigen::VectorXd p, Eigen::VectorXd q_0, Eigen::VectorXd weight){
	
	// Neighborhood distance of each parameters
	double diff_K = q_0.maxCoeff() / 100;
  	double diff_p = (p.maxCoeff() - p.minCoeff()) / 100;
  	double diff_var = 50;
  	
  	// Initialization of parameters
  	double K_min = par(0);
  	double K_max = par(1);
  	double p_0 = par(2);
  	double Var = par(3);
  	//std::cout << K_min << " " << K_max << " " << p_0 << " " << Var << std::endl;
  	
  	int same = 0;
  	Eigen::Vector3d error;
  	Eigen::VectorXd q(q_0.size());
  	
  	int max_iter = 250;
  	int count = 0;
  	
  	while(same != 4 && count <= max_iter){
  		count += 1;
  		same = 0;
  		
  		// Scanning neighboring minimum capacity
  		error << 0, 0, 0;
		for(int item = 0; item < q.size(); ++ item){
			if(q(item) != -9999 && p(item) != -9999){
				q(item) = (K_min - diff_K) + (K_max + diff_K - K_min) * cdf(norm_dist, (p(item) - p_0) / pow(Var, .5));	
				error(0) += weight(item) * abs(q(item) - q_0(item)) * pdf(norm_dist, (p(item) - p_0 + diff_p) / pow(Var, .5));
				//std::cout << abs(q(item) - q_0(item)) << std::endl;
				q(item) = K_min + (K_max - K_min) * cdf(norm_dist, (p(item) - p_0) / pow(Var, .5));
				error(1) += weight(item) * abs(q(item) - q_0(item)) * pdf(norm_dist, (p(item) - p_0 + diff_p) / pow(Var, .5));
				q(item) = (K_min + diff_K) + (K_max - diff_K - K_min) * cdf(norm_dist, (p(item) - p_0) / pow(Var, .5));
				error(2) += weight(item) * abs(q(item) - q_0(item)) * pdf(norm_dist, (p(item) - p_0 + diff_p) / pow(Var, .5));	
			}
		}
		K_min += (error(2) <= error(0)) * (error(2) < error(1)) * diff_K; 
		K_min -= (error(0) < error(1)) * (error(0) < error(2)) * diff_K;
  		same += (error(1) <= error(0)) * (error(1) <= error(2));
  		
  		// Scanning neighboring maximum capacity
  		error << 0, 0, 0;
		for(int item = 0; item < q.size(); ++ item){
			if(q(item) != -9999 && p(item) != -9999){
				q(item) = K_min + (K_max - diff_K - K_min) * cdf(norm_dist, (p(item) - p_0) / pow(Var, .5));	
				error(0) += weight(item) * abs(q(item) - q_0(item)) * pdf(norm_dist, (p(item) - p_0 + diff_p) / pow(Var, .5));
				q(item) = K_min + (K_max - K_min) * cdf(norm_dist, (p(item) - p_0) / pow(Var, .5));
				error(1) += weight(item) * abs(q(item) - q_0(item)) * pdf(norm_dist, (p(item) - p_0 + diff_p) / pow(Var, .5));
				q(item) = K_min + (K_max + diff_K - K_min) *  cdf(norm_dist, (p(item) - p_0) / pow(Var, .5));
				error(2) += weight(item) * abs(q(item) - q_0(item)) * pdf(norm_dist, (p(item) - p_0 + diff_p) / pow(Var, .5));	
			}
		}
		K_max += (error(2) <= error(0)) * (error(2) < error(1)) * diff_K; 
		K_max -= (error(0) < error(1)) * (error(0) < error(2)) * diff_K;
  		same += (error(1) <= error(0)) * (error(1) <= error(2));
  		
  		// Scanning neighboring critical price
  		error << 0, 0, 0;
		for(int item = 0; item < q.size(); ++ item){
			if(q(item) != -9999 && p(item) != -9999){
				q(item) = K_min + (K_max - K_min) * cdf(norm_dist, (p(item) - p_0 + diff_p) / pow(Var, .5));
				error(0) += weight(item) * abs(q(item) - q_0(item)) * pdf(norm_dist, (p(item) - p_0 + diff_p) / pow(Var, .5));
				q(item) = K_min + (K_max - K_min) * cdf(norm_dist, (p(item) - p_0) / pow(Var, .5));
				error(1) += weight(item) * abs(q(item) - q_0(item)) * pdf(norm_dist, (p(item) - p_0 + diff_p) / pow(Var, .5));
				q(item) = K_min + (K_max - K_min) * cdf(norm_dist, (p(item) - p_0 - diff_p) / pow(Var, .5));
				error(2) += weight(item) * abs(q(item) - q_0(item)) * pdf(norm_dist, (p(item) - p_0 + diff_p) / pow(Var, .5));
			}
		}
		p_0 += (error(2) <= error(0)) * (error(2) < error(1)) * diff_p;
		p_0 -= (error(0) < error(1)) * (error(0) < error(2)) * diff_p;
  		same += (error(1) <= error(0)) * (error(1) <= error(2));
  		
  		// Scanning neighboring variance
  		error << 0, 0, 0;
		for(int item = 0; item < q.size(); ++ item){
			if(q(item) != -9999 && p(item) != -9999){
				q(item) = K_min + (K_max - K_min) * cdf(norm_dist, (p(item) - p_0) / pow(Var - diff_var, .5));
				error(0) += weight(item) * abs(q(item) - q_0(item)) * pdf(norm_dist, (p(item) - p_0 + diff_p) / pow(Var, .5));
				q(item) = K_min + (K_max - K_min) * cdf(norm_dist, (p(item) - p_0) / pow(Var, .5));
				error(1) += weight(item) * abs(q(item) - q_0(item)) * pdf(norm_dist, (p(item) - p_0 + diff_p) / pow(Var, .5));
				q(item) = K_min + (K_max - K_min) * cdf(norm_dist, (p(item) - p_0) / pow(Var + diff_var, .5));
				error(2) += weight(item) * abs(q(item) - q_0(item)) * pdf(norm_dist, (p(item) - p_0 + diff_p) / pow(Var, .5));
			}
		}
		Var += (error(2) <= error(0)) * (error(2) < error(1)) * diff_var;
		Var -= (error(0) < error(1)) * (error(0) < error(2)) * diff_var;
		if(Var <= 0){
			break;
		}
		else{
			same += (error(1) <= error(0)) * (error(1) <= error(2));
		}
	}
	
	Eigen::VectorXd par_new(4);
	if(count <= max_iter && K_max > K_min && Var > 0){
		par_new << K_min, K_max, p_0, Var;
	}
	else{
		par_new << -9999, -9999, -9999, -9999;
	}
	
	return par_new;
}

int main(){
	auto start = std::chrono::high_resolution_clock::now();
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast <std::chrono::microseconds> (stop - start);
	
	std::string fin_name;
	// Read day ahead market prices data
	fin_name = "input/da_price_2021.csv";
	int num_row = 8760; 
	int num_col = 14;
	Eigen::MatrixXd da_price = read_file(num_row, num_col, fin_name);
	
	// Read day ahead market prices data
	fin_name = "input/residual_load_forecast_2021.csv";
	Eigen::MatrixXd residual_load_forecast = read_file(num_row, num_col, fin_name);
	
	// Initialization of parameters
	int iteration_num = 0;
	double decay_time = 1440;
	std::string col_names[num_col - 1] = {"NO1","NO2","NO3","NO4","NO5","DE-LU","DK1","FI","GB","NL","SE1","SE2","SE3"};
	std::vector<std::string> col_output = {"Minimum Capacity", "Maximum Capacity", "Price_0", "Variance"};
	std::string fout_name;
	Eigen::VectorXd par(4);
	Eigen::MatrixXd par_ts(num_row, 4);
	
	Eigen::MatrixXd weight(num_row, num_row);
	for(int tick = 0; tick < num_row; ++ tick){
		weight.col(tick) = exp(-abs(tick - da_price.col(0).array()) / decay_time);
	}
	
	for(int item = 0; item < num_col - 1; ++ item){
		
		std::cout << "-------------------------------------------" << std::endl;
		std::cout << col_names[item] << std::endl;
		std::cout << "-------------------------------------------" << std::endl;
		par << 0, residual_load_forecast.col(item + 1).maxCoeff(), da_price.col(item + 1).mean(), 1000;
		
		for(int tick = 0; tick < num_row; ++ tick){
			iteration_num = 0;	
			if(par(3) != -9999){
				par = Neighbor_Search(par, da_price.col(item + 1), residual_load_forecast.col(item + 1), weight.col(tick));
			}	
			while(par(3) == -9999){
				if(iteration_num <= 1){
					iteration_num += 1;
					par << 0, (rand() % 100) / 100 * residual_load_forecast.col(item + 1).maxCoeff(), da_price.col(item + 1).mean(), 1000;
					par = Neighbor_Search(par, da_price.col(item + 1), residual_load_forecast.col(item + 1), weight.col(tick));					
				}
				else{
					std::cout << "Tick: " << tick << "; An error occurred" << std::endl;
					break;
					//par = par_ts.row(tick - 1);
				}		
			}
			par_ts.row(tick) = par;
			
			if(tick % 1000 == 0){
				stop = std::chrono::high_resolution_clock::now();
				duration = std::chrono::duration_cast <std::chrono::microseconds> (stop - start);				
				std::cout << "Tick: " << tick << "; Time Laspe: " << duration.count() / pow(10, 6) << " seconds" << std::endl;
			}
		}
		
		fout_name = "output/merit_order_curve_" + col_names[item] + ".csv";
		write_file(par_ts, fout_name, col_output);
		std::cout << "\n" << std::endl;
	}

}