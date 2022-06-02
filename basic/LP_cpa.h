// Geostat functions header file
#pragma once

#ifndef LP_CPA
#define LP_CPA

#include <iostream>
#include <omp.h>
#include "../basic/Basic_Definitions.h"
#include "../basic/Eigen_Sparse.h"

// Structures
// Constraint object
struct LP_constraint{
	Eigen::VectorXd eq_norm;
	Eigen::VectorXd ie_increment;										// Increment per unit variation in the gradient direction
	Eigen::SparseMatrix <double> eq_matrix;								// Coefficients for equality constraints
	Eigen::SparseMatrix <double> ie_matrix;								// Coefficients for inequality constraints
	Eigen::SparseMatrix <double> eq_cov_matrix;							// Dot product of the equality constraints
	Eigen::SimplicialLDLT <Eigen::SparseMatrix <double>> eq_cov_solver; // Solver for the covariance (dot product) matrix
};

// Boundary object
struct LP_boundary{
	Eigen::VectorXd eq_vector;			// Value of original equality constraints
	Eigen::MatrixXd ie_matrix;			// Boundary of original inequality constraints; 0th column: lower bound value; 1st column: upper bound value
};

// Objective object
struct LP_objective{
	Eigen::VectorXd vector;
	double value;
};

// Solver object
struct LP_solver{
	Eigen::SimplicialLDLT <Eigen::SparseMatrix <double>> Spd; 	// Solver for a symmetric positive definite matrix
	Eigen::SparseQR <Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering <int>> qr;
};

// Projection gradient object
struct LP_proj_gradient{
	Eigen::VectorXd vector_default; 	// The projection of objective vector onto the equality constraints
};

// Solution object
struct LP_solution{
	Eigen::VectorXd vector;
	Eigen::VectorXd vector_ref;
	Eigen::VectorXd vector_ini;
};

struct LP_object{
	// Input parameters
	int Constraints_eq_num;
	int Constraints_ie_num;
	int Variables_num;
	
	// Mixed substructure
	LP_objective Objective;
	LP_constraint Constraint;
	LP_boundary Boundary;
	LP_solver Solver;
	LP_proj_gradient Proj_grad;
	
	// Process and output variables
	LP_solution Solution;
};

// Functions
void LP_constraint_eq_normalization(LP_object&);
void LP_boundary_eq_normalization(LP_object&);
void LP_constraint_cov_eq_matrix_generation(LP_object&);
void LP_default_gradient(LP_object&);
void LP_optimization(LP_object&);
void LP_result_print(LP_object&, std::string);
void LP_process(LP_object&, std::string Problem_name = "Linear Problem", bool print_result = 1, bool constraint_flag = 1, bool boundary_flag = 1, bool objective_flag = 1);

#endif