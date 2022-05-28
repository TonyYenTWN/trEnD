// Source file for linear programming using gradient projection algorithm
#include <iostream>
#include <omp.h>
#include "../basic/Basic_Definitions.h"
#include <Eigen/Sparse>

typedef Eigen::Triplet <double> Trip;    	// Define a triplet object

struct LP_object{
	int num_eq_constraints;
	int num_ie_constraints;
	int num_variables;
	Eigen::VectorXd Objective_vector;
	Eigen::VectorXd Ortho_proj_length_vector;
	Eigen::VectorXd Proj_grad_vector;
	Eigen::VectorXd Boundary_vector;
	Eigen::VectorXd Solution_vector;
	Eigen::SparseMatrix <double> Constraint_matrix;
};

void LP_optimized(LP_object &Problem){
	
}

int main(){
	LP_object test;
	test.num_eq_constraints = 0;
	test.num_ie_constraints = 7;
	test.num_variables = 4;
	LP_optimized(test);
}