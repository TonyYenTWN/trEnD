// Main Source File
#include <iostream>
#include <chrono>
#include <omp.h>
#include <boost/math/distributions/normal.hpp>
#include "../basic/rw_csv.cpp"

boost::math::normal norm_dist(0.0, 1.0);

Eigen::Vector3d Neighbor_Search(Eigen::Vector3d par, Eigen::VectorXd p, Eigen::VectorXd q_0, Eigen::VectorXd weight){
	
	// Neighborhood distance of each parameters
	double diff_K = 50;
  	double diff_p = 1;
  	double diff_var = 10;
  	
  	// Initialization of parameters
  	double K = par(0);
  	double p_0 = par(1);
  	double Var = par(2);
  	
  	int same = 0;
  	Eigen::Vector3d error;
  	Eigen::VectorXd q(q_0.size());
  	
  	while(same != 3){
  		same = 0;
  		error << 0, 0, 0;
  		
  		// Scanning neighboring capacity
		for(int item = 0; item < q.size(); ++ item){
			if(q(item) != -9999){
				q(item) = (K - diff_K) * cdf(norm_dist, (p(item) - p_0) / pow(Var, .5));	
				error(0) += weight(item) * pow(q(item) - q_0(item), 2);
				q(item) = K * cdf(norm_dist, (p(item) - p_0) / pow(Var, .5));
				error(1) += weight(item) * pow(q(item) - q_0(item), 2);
				q(item) = (K + diff_K) * cdf(norm_dist, (p(item) - p_0) / pow(Var, .5));
				error(2) += weight(item) * pow(q(item) - q_0(item), 2);	
			}
		}
		K += (error(2) < error(0)) * (error(2) < error(1)) * diff_K;
		K -= (error(0) < error(1)) * (error(0) < error(2)) * diff_K;
  		same += (error(1) <= error(0)) * (error(1) <= error(2));
  		
  		// Scanning neighboring critical price
		for(int item = 0; item < q.size(); ++ item){
			if(q(item) != -9999){
				q(item) = K * cdf(norm_dist, (p(item) - p_0 + diff_p) / pow(Var, .5));
				error(0) += weight(item) * pow(q(item) - q_0(item), 2);
				q(item) = K * cdf(norm_dist, (p(item) - p_0) / pow(Var, .5));
				error(1) += weight(item) * pow(q(item) - q_0(item), 2);
				q(item) = K * cdf(norm_dist, (p(item) - p_0 - diff_p) / pow(Var, .5));
				error(2) += weight(item) * pow(q(item) - q_0(item), 2);
			}
		}
		p_0 += (error(2) < error(0)) * (error(2) < error(1)) * diff_p;
		p_0 -= (error(0) < error(1)) * (error(0) < error(2)) * diff_p;
  		same += (error(1) <= error(0)) * (error(1) <= error(2));
  		
  		// Scanning neighboring variance
		for(int item = 0; item < q.size(); ++ item){
			if(q(item) != -9999){
				q(item) = K * cdf(norm_dist, (p(item) - p_0) / pow(Var - diff_var, .5));
				error(0) += weight(item) * pow(q(item) - q_0(item), 2);
				q(item) = K * cdf(norm_dist, (p(item) - p_0) / pow(Var, .5));
				error(1) += weight(item) * pow(q(item) - q_0(item), 2);
				q(item) = K * cdf(norm_dist, (p(item) - p_0) / pow(Var + diff_var, .5));
				error(2) += weight(item) * pow(q(item) - q_0(item), 2);
			}
		}
		Var += (error(2) < error(0)) * (error(2) < error(1)) * diff_var;
		Var -= (error(0) < error(1)) * (error(0) < error(2)) * diff_var;
  		same += (error(1) <= error(0)) * (error(1) <= error(2));
	}
	
	Eigen::Vector3d par_new;
	par_new << K, p_0, Var;
	
	return par_new;
}

int main(){
	
	std::string fin_name;
	// Read day ahead market prices data
	fin_name = "input/da_price_2021.csv";
	int num_row = 8760; 
	int num_col = 14;
	Eigen::MatrixXd da_price = read_file(num_row, num_col, fin_name);
	
	// Read day ahead market prices data
	fin_name = "input/residual_load_forecast_2021.csv";
	num_row = 8760; 
	num_col = 14;
	Eigen::MatrixXd residual_load_forecast = read_file(num_row, num_col, fin_name);
	
	// Initialization of parameters
	Eigen::Vector3d par;
	par << 10000, 100, 100;
	
	Eigen::VectorXd weight = Eigen::VectorXd::Ones(num_row);
	par = Neighbor_Search(par, da_price.col(4), residual_load_forecast.col(4), weight);	
	
	std::cout << par << std::endl;
}