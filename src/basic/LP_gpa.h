// Header file for linear programming using gradient projection algorithm
#pragma once

#ifndef LP_GPA
#define LP_GPA

#include <iostream>
#include <iomanip>
#include <omp.h>
#include "../basic/Basic_Definitions.h"
#include "../basic/Eigen_Sparse.h"

// Constraint object
struct LP_constraint{
	Eigen::VectorXd norm;
	Eigen::VectorXd ie_reduced_increment;				// Increment per unit variation in the gradient direction
	Eigen::SparseMatrix <double> eq_orig_matrix;		// Coefficients for original equality constraints
	Eigen::SparseMatrix <double> ie_orig_matrix;		// Coefficients for original inequality constraints
	Eigen::SparseMatrix <double> ie_reduced_matrix;		// Coefficients for reduced inequality constraints
	Eigen::SparseMatrix <double> permutation_matrix; 	// Record of permutation order
	Eigen::SparseMatrix <double> ie_reduced_cov_matrix;	// Dot product of reduced inequality constraints
};

// Boundary object
struct LP_boundary{
	Eigen::VectorXd eq_vector;				// Value of original equality constraints
	Eigen::MatrixXd ie_orig_matrix;			// Boundary of original inequality constraints; 0th column: lower bound value; 1st column: upper bound value
	Eigen::MatrixXd ie_reduced_matrix;		// Boundary of reduced inequality constraints
};

// Objective object
struct LP_objective{
	Eigen::VectorXd orig_vector;
	Eigen::VectorXd reduced_vector;
	Eigen::VectorXd ie_reduced_cov_vector;
	Eigen::VectorXd varying_vector;			// Indicates whether the coefficients of the objective is stepwise-linear for a variable 
	Eigen::VectorXd update_coeff;			// Indicates the stepwise-linear coefficient to be updated next loop
	double orig_value;
	double orig_value_sum;			// Sum of stepwise-linear optimization iterations
	double reduced_value;
};

// Solver object
struct LP_matrix_solver{
	Eigen::SimplicialLDLT <Eigen::SparseMatrix <double>> ldlt; 								// LDLT solver for a symmetric positive definite matrix
	Eigen::SparseLU <Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering <int>> lu;			// LU solver for an arbitrary matrix
	Eigen::SparseQR <Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering <int>> qr;			// QR solver for the transpose of an arbitrary matrix
};

// Solution object
struct LP_solution{
	Eigen::VectorXd orig_vector;
	Eigen::VectorXd reduced_vector;
};

struct LP_object{
	// Input parameters
	int Constraints_eq_num;
	int Constraints_ie_num;
	int Variables_num;
	LP_objective Objective;
	LP_constraint Constraint;
	LP_boundary Boundary;
	LP_matrix_solver Solver;
	
	// Process and output variables
	Eigen::VectorXd Proj_grad_vector;
	LP_solution Solution;
};

// Functions
void LP_constraint_eq_redundant_deletion(LP_object&);
void LP_variables_permutation(LP_object&, bool stepwise_obj = 0);
void LP_solution_permutation(LP_object&, bool inverse = 0);
void LP_constraint_redundant_matrix_solver(LP_object&);
void LP_constraint_ie_reduced_generation(LP_object&);
void LP_boundary_ie_reduced_generation(LP_object&);
void LP_objective_reduced_generation(LP_object&);
void LP_feasible_solution_reduced_generation(LP_object&);
void LP_constraint_ie_reduced_normalization(LP_object&);
void LP_boundary_ie_reduced_normalization(LP_object&);
void LP_constraint_ie_reduced_cov_matrix_generation(LP_object&);
void LP_optimization(LP_object&, bool stepwise_obj = 0);
void LP_result_print(LP_object&, std::string);
void LP_process(LP_object &Problem, std::string Problem_name = "Linear Problem", bool result_output = 1, bool find_sol = 1, bool stepwise_obj = 0, bool constraint_update = 1, bool boundary_update = 1, bool objective_update = 1);

#endif