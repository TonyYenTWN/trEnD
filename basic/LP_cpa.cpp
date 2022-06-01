// Source file for linear programming using contraction and gradient projection algorithm
// Applicable when an initial feasible point in the interior is trivial

#include <iostream>
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
	double orig_value;
	double reduced_value;
};

// Solver object
struct LP_matrix_solver{
	Eigen::SimplicialLDLT <Eigen::SparseMatrix <double>> Spd; 									// Solver for a symmetric positive definite matrix
	Eigen::SparseLU <Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering <int>> Normal;			// Solver for an arbitrary matrix
	Eigen::SparseLU <Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering <int>> Normal_Trans;	// Solver for the transpose of an arbitrary matrix
};

// Solution object
struct LP_solution{
	Eigen::VectorXd orig_vector;
	Eigen::VectorXd reduced_vector;
	Eigen::VectorXd reduced_vector_ref;
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

// Delete redundant equality constraints
void LP_constraint_eq_redundant_deletion(LP_object &Problem){
	// Set precision of zero detection
	double tol = pow(10, -12);
	
	// Covariance (dot product) of the equality constraints in the original space
	Eigen::SparseMatrix <double> Constraint_eq_orig_cov = Problem.Constraint.eq_orig_matrix * Problem.Constraint.eq_orig_matrix.transpose();
	Problem.Solver.Spd.compute(Constraint_eq_orig_cov);
	
	// Check if the coavriance matrix is not of full rank (determinant = 0)
	if(abs(Problem.Solver.Spd.determinant()) < tol){
		// Declare a dynamic vector of list of active (non-redeundant) equality constraints
		std::vector<int> active_const_seq;
		active_const_seq.reserve(Problem.Constraints_eq_num);
		Eigen::VectorXd Constraint_eq_cov_D = Problem.Solver.Spd.vectorD();
		for(int item_iter = 0; item_iter < Problem.Constraints_eq_num; ++ item_iter){
			if(abs(Constraint_eq_cov_D(item_iter)) >= tol){
				active_const_seq.push_back(item_iter);
			}
		}
		
		std::vector<Trip> Constraint_eq_sub_trip;
		Constraint_eq_sub_trip.reserve(active_const_seq.size());
		Eigen::SparseMatrix <double> Span_reduction(active_const_seq.size(), Problem.Constraints_eq_num);
		for(int item_iter = 0; item_iter < active_const_seq.size(); ++ item_iter){
			Constraint_eq_sub_trip.push_back(Trip(active_const_seq[item_iter], active_const_seq[item_iter], 1));
		}
		Span_reduction.setFromTriplets(Constraint_eq_sub_trip.begin(), Constraint_eq_sub_trip.end());

		// Update the information about equality constraints after removing the redundant ones
		Problem.Constraint.eq_orig_matrix = Span_reduction * Problem.Constraint.eq_orig_matrix;
		Problem.Boundary.eq_vector = Span_reduction * Problem.Boundary.eq_vector;
		Problem.Constraints_eq_num = active_const_seq.size();
	}
}

// Find redundant variables forming full rank square matrix from the equality constraints
void LP_variables_permutation(LP_object &Problem){
	// Set precision of zero detection
	double tol = pow(10, -12);
	
	// Check which of the original variables should be chosen as redundant
	Problem.Solver.Spd.compute(Problem.Constraint.eq_orig_matrix.transpose() * Problem.Constraint.eq_orig_matrix);
	std::vector<int> redundant_var_seq;
	redundant_var_seq.reserve(Problem.Constraints_eq_num);
	Eigen::VectorXd Constraint_eq_cov_D = Problem.Solver.Spd.vectorD();
	for(int var_id = Problem.Variables_num - 1; var_id > -1; -- var_id){
		if(abs(Constraint_eq_cov_D(var_id)) >= tol){
			redundant_var_seq.push_back(var_id);
			
			// Exit for loop once sufficient redundant variables are gathered
			if(redundant_var_seq.size() == Problem.Constraints_eq_num){
				break;
			}
		}
	}
	std::reverse(redundant_var_seq.begin(), redundant_var_seq.end());
	
	// Set a permutation matrix so that the redundant variables are put to the end
	std::vector<Trip> Var_permut_trip;
	Var_permut_trip.reserve(redundant_var_seq.size());
}

// Linear system for reducing the redundant variables
void LP_constraint_redundant_matrix_solver(LP_object &Problem){
	Problem.Solver.Normal.compute(Problem.Constraint.eq_orig_matrix.rightCols(Problem.Constraints_eq_num));
	Problem.Solver.Normal_Trans.compute(Problem.Constraint.eq_orig_matrix.rightCols(Problem.Constraints_eq_num).transpose());
}

void LP_constraint_ie_reduced_generation(LP_object &Problem){
	// Set default value of reduced constraint matrix
	Eigen::MatrixXd reduced_matrix = Eigen::MatrixXd(Problem.Constraint.ie_orig_matrix.leftCols(Problem.Variables_num - Problem.Constraints_eq_num));
	
	// Check if there are eqaulity constraints in the original problem
	if(Problem.Constraints_eq_num != 0){
		// If equality constraints exist, replacing the redundant variables with them
		// A_e * x_e + A_r * x_r = b
		// A_e * x_e - b = -A_r * x_r
		// l_r <= -A_r^(-1) * (A_e * x_e - b) = x_r <= u_r
		// l_r - A_r^(-1) * b <= -A_r^(-1) * A_e * x_e <= u_r - A_r^(-1) * b
		Eigen::VectorXd col_orig;
		for(int var_id = 0; var_id < Problem.Variables_num - Problem.Constraints_eq_num; ++ var_id){
			col_orig = Problem.Constraint.eq_orig_matrix.col(var_id);
			reduced_matrix.middleRows(Problem.Variables_num - Problem.Constraints_eq_num, Problem.Constraints_eq_num).col(var_id) = -Problem.Solver.Normal.solve(col_orig);
		}
						
		// Check if there are ineqaulity constraints in the original problem
		if(Problem.Constraints_ie_num != 0){
			// If inequality constraints exist in the original, renew the coefficients
			// l_i <= M_e * x_e + M_r * x_r <= u_i
			// l_i <= M_e * x_e - M_r * A_r^(-1) * (A_e * x_e - b) <= u_i
			// l_i - M_r * A_r^(-1) * b <= (M_e - M_r * A_r^(-1) * A_e) * x_e <= u_i - M_r * A_r^(-1) * b			
		}
	}
	
	double tol = pow(10, -8);
	Problem.Constraint.ie_reduced_matrix = reduced_matrix.sparseView(tol, 1);
}

void LP_boundary_ie_reduced_generation(LP_object &Problem){
	// Set default value of reduced boundary vector
	Problem.Boundary.ie_reduced_matrix = Problem.Boundary.ie_orig_matrix;
	
	// Check if there are eqaulity constraints in the original problem
	if(Problem.Constraints_eq_num != 0){
		// If equality constraints exist, replacing the redundant variables with them
		// A_e * x_e + A_r * x_r = b
		// A_e * x_e - b = -A_r * x_r
		// l_r <= -A_r^(-1) * (A_e * x_e - b) = x_r <= u_r
		// l_r - A_r^(-1) * b <= -A_r^(-1) * A_e * x_e <= u_r - A_r^(-1) * b
		Problem.Boundary.ie_reduced_matrix.middleRows(Problem.Variables_num - Problem.Constraints_eq_num, Problem.Constraints_eq_num).col(0) -= Problem.Solver.Normal.solve(Problem.Boundary.eq_vector);
		Problem.Boundary.ie_reduced_matrix.middleRows(Problem.Variables_num - Problem.Constraints_eq_num, Problem.Constraints_eq_num).col(1) -= Problem.Solver.Normal.solve(Problem.Boundary.eq_vector);
						
		// Check if there are ineqaulity constraints in the original problem
		if(Problem.Constraints_ie_num != 0){			
		}	
	}
}

void LP_objective_reduced_generation(LP_object &Problem){
	// Set default value of reduced objective vector
	Problem.Objective.reduced_vector = Problem.Objective.orig_vector.head(Problem.Variables_num - Problem.Constraints_eq_num);
	
	// Check if there are eqaulity constraints in the original problem
	if(Problem.Constraints_eq_num != 0){
		// Z = c_e^T * x_e + c_r^T * x_r 
		//   = c_e^T * x_e - c_r^T * A_r^(-1) * (A_e * x_e - b)
		//   = (c_e^T - c_r^T * A_r^(-1) * A_e) * x_e + c_r^T * A_r^(-1) * b
		Problem.Objective.reduced_vector = Problem.Objective.orig_vector.head(Problem.Variables_num - Problem.Constraints_eq_num);
		Eigen::VectorXd coeff_redundant = Problem.Objective.orig_vector.tail(Problem.Constraints_eq_num);
		Problem.Objective.reduced_vector -= (Problem.Constraint.eq_orig_matrix.leftCols(Problem.Variables_num - Problem.Constraints_eq_num)).transpose() * Problem.Solver.Normal_Trans.solve(coeff_redundant);		
	}
}

// Normalization of reduced inequality constraints
void LP_constraint_ie_reduced_normalization(LP_object &Problem){
	Problem.Constraint.norm = Eigen::VectorXd::Ones(Problem.Variables_num + Problem.Constraints_ie_num);
	
	for(int row_id = Problem.Variables_num - Problem.Constraints_eq_num; row_id < Problem.Constraint.ie_reduced_matrix.rows(); ++ row_id){
		Problem.Constraint.norm(row_id) = Problem.Constraint.ie_reduced_matrix.row(row_id).norm();
		Problem.Constraint.ie_reduced_matrix.row(row_id) /= Problem.Constraint.norm(row_id);
	}
}

// Normalization of reduced inequality boundaries
void LP_boundary_ie_reduced_normalization(LP_object &Problem){
	for(int row_id = Problem.Variables_num - Problem.Constraints_eq_num; row_id < Problem.Constraint.ie_reduced_matrix.rows(); ++ row_id){
		Problem.Boundary.ie_reduced_matrix.row(row_id) /= Problem.Constraint.norm(row_id);
	}
}

// Function for the optimization algorithm
void LP_optimization(LP_object &Problem){
	// Parameters declaration
	double tol = pow(10, -12);
	double eps = pow(10, -6); 
	
	// Process variables (in the loop) declaration
	double min_increment;
	double Previous_obj_value;
	double Obj_value;
	Eigen::MatrixXd Boundary_gap(Problem.Variables_num + Problem.Constraints_ie_num, 2); 
	Eigen::VectorXd Projected_increment(Problem.Variables_num + Problem.Constraints_ie_num);	
	Eigen::VectorXd Projected_gradient;
	Eigen::VectorXd Previous_solution;
	Eigen::VectorXd Contracted_solution;
	Eigen::Vector2d current_increment;

	// Initialization of starting point and other process variables
	Problem.Solution.reduced_vector = Eigen::VectorXd::Zero(Problem.Variables_num - Problem.Constraints_eq_num);
	Problem.Objective.reduced_value = Problem.Solution.reduced_vector.dot(Problem.Objective.reduced_vector);
	
	// Main loop for searching the optimal solution
	while(1){
		// Update contracted solution and previous solution for later use in the loop
		Previous_obj_value = Problem.Objective.reduced_value;
		Previous_solution = Problem.Solution.reduced_vector;
		Contracted_solution = Problem.Solution.reduced_vector_ref;
		Contracted_solution += (1 - eps) * (Problem.Solution.reduced_vector - Problem.Solution.reduced_vector_ref);
		
		// Calculate the allowed increment along the objective direction before reaching a constraint
		Boundary_gap.col(0) = Problem.Constraint.ie_reduced_matrix * Contracted_solution - Problem.Boundary.ie_reduced_matrix.col(0);
		Boundary_gap.col(1) = Problem.Boundary.ie_reduced_matrix.col(1) - Problem.Constraint.ie_reduced_matrix * Contracted_solution;
		min_increment = std::numeric_limits<double>::infinity();
		for(int constraint_id = 0; constraint_id < Boundary_gap.rows(); ++ constraint_id){
			current_increment(0) = -Boundary_gap(constraint_id, 0) / Problem.Constraint.ie_reduced_increment(constraint_id);
			if(current_increment(0) > tol && current_increment(0) < min_increment){
				min_increment = current_increment(0);
			}
			current_increment(1) = Boundary_gap(constraint_id, 1) / Problem.Constraint.ie_reduced_increment(constraint_id);
			if(current_increment(1) > tol && current_increment(1) < min_increment){
				min_increment = current_increment(1);
			}		
		}
		
		// Update the solution point and projected gradient on / near the active constraint(s)
		Problem.Solution.reduced_vector = Contracted_solution + min_increment * Problem.Objective.reduced_vector;
		Projected_gradient = Problem.Solution.reduced_vector - Previous_solution;
		Projected_gradient /= Projected_gradient.norm();
		
		// Calculate the allowed increment along the projected gradient before reaching another set of constraint(s)
		Boundary_gap.col(0) = Problem.Constraint.ie_reduced_matrix * Problem.Solution.reduced_vector - Problem.Boundary.ie_reduced_matrix.col(0);
		Boundary_gap.col(1) = Problem.Boundary.ie_reduced_matrix.col(1) - Problem.Constraint.ie_reduced_matrix * Problem.Solution.reduced_vector;
		Projected_increment = Problem.Constraint.ie_reduced_matrix * Projected_gradient;
		min_increment = std::numeric_limits<double>::infinity();
		for(int constraint_id = 0; constraint_id < Boundary_gap.rows(); ++ constraint_id){
			if(abs(Projected_increment(constraint_id)) > tol){
				current_increment(0) = -Boundary_gap(constraint_id, 0) / Projected_increment(constraint_id);
				if(current_increment(0) > tol && current_increment(0) < min_increment){
					min_increment = current_increment(0);
				}
				current_increment(1) = Boundary_gap(constraint_id, 1) / Projected_increment(constraint_id);
				if(current_increment(1) > tol && current_increment(1) < min_increment){
					min_increment = current_increment(1);
				}			
			}			
		}
		Problem.Solution.reduced_vector += min_increment * Projected_gradient;
		Obj_value = Problem.Solution.reduced_vector.dot(Problem.Objective.reduced_vector);
		
		// Check if objective function can be improved or optimality has been reached
		if(Obj_value >= Previous_obj_value){
			Problem.Objective.reduced_value = Obj_value;
		}
		else{
			// If current objective value is smaller than previous, optimality occurs at previous solution
			Problem.Solution.reduced_vector = Previous_solution;
			break;
		}		
	}
	
	// Reconstruct the original solution and objective value from the reduced problem, if necessary
	if(Problem.Constraints_eq_num != 0){
		// Reconstruction of the solution
		// x_r = -A_r^(-1) * (A_e * x_e - b)
		// Problem.Constraint.eq_orig_matrix
		Eigen::VectorXd redundant_variables = -Problem.Solver.Normal.solve(Problem.Constraint.eq_orig_matrix.leftCols(Problem.Variables_num - Problem.Constraints_eq_num) * Problem.Solution.reduced_vector - Problem.Boundary.eq_vector);
		Problem.Solution.orig_vector = Eigen::VectorXd(Problem.Variables_num);
		Problem.Solution.orig_vector << Problem.Solution.reduced_vector, redundant_variables;
		
		// Reconstruction of the objective value
		// Z = c_e^T * x_e + c_r^T * x_r 
		//   = c_e^T * x_e - c_r^T * A_r^(-1) * (A_e * x_e - b)
		//   = (c_e^T - c_r^T * A_r^(-1) * A_e) * x_e + c_r^T * A_r^(-1) * b
		Problem.Objective.orig_value = 	Problem.Objective.reduced_value + Problem.Objective.orig_vector.tail(Problem.Constraints_eq_num).transpose() * Problem.Solver.Normal.solve(Problem.Boundary.eq_vector);
	}
	else{
		Problem.Solution.orig_vector = Problem.Solution.reduced_vector;
	}
}

// Function for the optimization algorithm
void LP_result_print(LP_object &Problem, std::string Problem_name = "Linear Problem"){
	std::cout << "---------------------------------------------------------------" << std::endl;
	std::cout << "| " << Problem_name << std::endl;
	std::cout << "---------------------------------------------------------------" << std::endl;
	std::cout << "| Reduced Problem  | " << std::endl;	
	std::cout << "| Solution         | " << Problem.Solution.reduced_vector.transpose() << std::endl;
	std::cout << "| Objective value  | " << Problem.Objective.reduced_value << std::endl;
	std::cout << "---------------------------------------------------------------" << std::endl;
	std::cout << "| Original Problem | " << std::endl;
	std::cout << "| Solution         | " << Problem.Solution.orig_vector.transpose() << std::endl;
	std::cout << "| Objective value  | " << Problem.Objective.orig_value << std::endl;
	std::cout << "---------------------------------------------------------------" << std::endl;
	std::cout << "\n" << std::endl;
}

// Wrapping up all the process
void LP_process(LP_object &Problem, std::string Problem_name = "Linear Problem", bool constraint_update = 1, bool boundary_update = 1, bool objective_update = 1){
	// Generate sparse matrix for reduced inequality constraints, if needed
	if(constraint_update){
		LP_constraint_eq_redundant_deletion(Problem);
		LP_variables_permutation(Problem);
		LP_constraint_redundant_matrix_solver(Problem);
		LP_constraint_ie_reduced_generation(Problem);
		LP_boundary_ie_reduced_generation(Problem);
		LP_objective_reduced_generation(Problem);
		LP_constraint_ie_reduced_normalization(Problem);
		LP_boundary_ie_reduced_normalization(Problem);
		Problem.Constraint.ie_reduced_increment = Problem.Constraint.ie_reduced_matrix * Problem.Objective.reduced_vector;
	}
	else{
		// Update reduced inequality boudaries, if needed
		if(boundary_update){
			LP_boundary_ie_reduced_generation(Problem);
			LP_boundary_ie_reduced_normalization(Problem);
		}
		
		// Update reduced inequality increments along the gradient, if needed
		if(objective_update){
			Problem.Constraint.ie_reduced_increment = Problem.Constraint.ie_reduced_matrix * Problem.Objective.reduced_vector;
		}		
	} 
	
	// Solve the LP
	LP_optimization(Problem);
	
	// Print results
	LP_result_print(Problem, Problem_name);	
}

// Test case 1, with only equality constraints:
// maximize z: x_1 + 2 * x_2 + x_3 + x_4;
// s.t.
// 0 <= x_1 <= 1;
// 0 <= x_2 <= 1;
// 0 <= x_3 <= 1;
// 0 <= x_4 <= 1;
// 0 <= x_2 + x_3 <= 1; -> this is redundant
// x_1 + x_2 + 2 * x_3 + 3 * x_4 == 2;
// 2 * x_1 + x_2 + x_3 + x_4 == 1;
// Optimal solution is (0, .5, 0, .5) and optimal value is 1.5
// A trivial initial feasible point is (0, 0, 1, 0) or (.2, 0, 0, .6)
void test_problem_1_set(LP_object &Problem){
	// Set dimension of the problem
	Problem.Constraints_eq_num = 3;
	Problem.Constraints_ie_num = 0;
	Problem.Variables_num = 4;
	
	// Set objective vector
	Problem.Objective.orig_vector = Eigen::VectorXd(Problem.Variables_num);
	Problem.Objective.orig_vector << 1, 2, 1, 1;
	
	// Set boudary values for original equality and inequality constraints
	Problem.Boundary.eq_vector = Eigen::VectorXd(Problem.Constraints_eq_num);
	Problem.Boundary.eq_vector << 2, 1;
	Problem.Boundary.ie_orig_matrix = Eigen::MatrixXd(Problem.Variables_num + Problem.Constraints_ie_num, 2);
	Problem.Boundary.ie_orig_matrix.col(0) = Eigen::VectorXd::Zero(Problem.Variables_num + Problem.Constraints_ie_num);
	Problem.Boundary.ie_orig_matrix.col(1) = Eigen::VectorXd::Ones(Problem.Variables_num + Problem.Constraints_ie_num);
	
	// Set sparse matrix for original equality constraints
	std::vector<Trip> Constraint_eq_trip;
	Constraint_eq_trip.reserve(Problem.Constraints_eq_num * Problem.Variables_num);
	Constraint_eq_trip.push_back(Trip(0, 0, 1));
	Constraint_eq_trip.push_back(Trip(0, 1, 1));
	Constraint_eq_trip.push_back(Trip(0, 2, 2));
	Constraint_eq_trip.push_back(Trip(0, 3, 3));
	Constraint_eq_trip.push_back(Trip(1, 0, 2));
	Constraint_eq_trip.push_back(Trip(1, 1, 1));
	Constraint_eq_trip.push_back(Trip(1, 2, 1));
	Constraint_eq_trip.push_back(Trip(1, 3, 1));
	Constraint_eq_trip.push_back(Trip(2, 0, 2));
	Constraint_eq_trip.push_back(Trip(2, 1, 1));
	Constraint_eq_trip.push_back(Trip(2, 2, 1));
	Constraint_eq_trip.push_back(Trip(2, 3, 1));	
	Problem.Constraint.eq_orig_matrix = Eigen::SparseMatrix <double> (Problem.Constraints_eq_num, Problem.Variables_num);
	Problem.Constraint.eq_orig_matrix.setFromTriplets(Constraint_eq_trip.begin(), Constraint_eq_trip.end());
	
	// Set sparse matrix for original inequality constraints
	Problem.Constraint.ie_orig_matrix = Eigen::SparseMatrix <double> (Problem.Variables_num + Problem.Constraints_ie_num, Problem.Variables_num);
	std::vector<Trip> Constraint_ie_orig_trip;
	Constraint_ie_orig_trip.reserve((Problem.Constraints_ie_num + 1) * Problem.Variables_num);
	for(int var_id = 0; var_id < Problem.Variables_num; ++ var_id){
		Constraint_ie_orig_trip.push_back(Trip(var_id, var_id, 1));
	}
	Problem.Constraint.ie_orig_matrix.setFromTriplets(Constraint_ie_orig_trip.begin(), Constraint_ie_orig_trip.end());
	// Put additional inequality constraints here
	
	// Set reference point for contraction
	Problem.Solution.reduced_vector_ref = Eigen::VectorXd(Problem.Variables_num - Problem.Constraints_eq_num);
	Problem.Solution.reduced_vector_ref << 1. / 6., 1. / 15.;
}

// Test case 2, with only equality constraints:
// maximize z: 100 * x_1 + 13 * x_2 + 10 * x_3 + x_4 + x_5 + x_6;
// s.t.
// 0 <= x_1 <= 1;
// 0 <= x_2 <= 1;
// 0 <= x_3 <= 1;
// 0 <= x_4 <= 1;
// 0 <= x_5 <= 1;
// 0 <= x_6 <= 1;
// x_1 + 2 * x_2 + 2 * x_3 + 3 * x_4 == 2;
// 2 * x_1 + x_2 + x_3 + x_4 == 1;
// x_1 + x_2 + x_5 + x_6 == 1;
// Optimal solution is (.2, 0, 0, .6, 0, .8) and optimal value is 21.4
// 3 trivial feasible points is (0, 1, 0, 0, 0, 0), (.2, 0, 0, .6, .4, .4), or (0, 0, 1, 0, .5, .5)
// A feasible interior point for reference of contraction can be found by averaging the 3 trivial points
void test_problem_2_set(LP_object &Problem){
	// Set dimension of the problem
	Problem.Constraints_eq_num = 3;
	Problem.Constraints_ie_num = 0;
	Problem.Variables_num = 6;
	
	// Set objective vector
	Problem.Objective.orig_vector = Eigen::VectorXd(Problem.Variables_num);
	Problem.Objective.orig_vector << 100, 13, 10, 1, 1, 1;

	// Set boudary values for original equality and inequality constraints
	Problem.Boundary.eq_vector = Eigen::VectorXd(Problem.Constraints_eq_num);
	Problem.Boundary.eq_vector << 2, 1, 1;
	Problem.Boundary.ie_orig_matrix = Eigen::MatrixXd(Problem.Variables_num + Problem.Constraints_ie_num, 2);
	Problem.Boundary.ie_orig_matrix.col(0) = Eigen::VectorXd::Zero(Problem.Variables_num + Problem.Constraints_ie_num);
	Problem.Boundary.ie_orig_matrix.col(1) = Eigen::VectorXd::Ones(Problem.Variables_num + Problem.Constraints_ie_num);

	// Set sparse matrix for original equality constraints
	std::vector<Trip> Constraint_eq_trip;
	Constraint_eq_trip.reserve(Problem.Constraints_eq_num * Problem.Variables_num);
	Constraint_eq_trip.push_back(Trip(0, 0, 1));
	Constraint_eq_trip.push_back(Trip(0, 1, 2));
	Constraint_eq_trip.push_back(Trip(0, 2, 2));
	Constraint_eq_trip.push_back(Trip(0, 3, 3));
	Constraint_eq_trip.push_back(Trip(1, 0, 2));
	Constraint_eq_trip.push_back(Trip(1, 1, 1));
	Constraint_eq_trip.push_back(Trip(1, 2, 1));
	Constraint_eq_trip.push_back(Trip(1, 3, 1));
	Constraint_eq_trip.push_back(Trip(2, 0, 1));
	Constraint_eq_trip.push_back(Trip(2, 1, 1));
	Constraint_eq_trip.push_back(Trip(2, 4, 1));
	Constraint_eq_trip.push_back(Trip(2, 5, 1));
	Problem.Constraint.eq_orig_matrix = Eigen::SparseMatrix <double> (Problem.Constraints_eq_num, Problem.Variables_num);
	Problem.Constraint.eq_orig_matrix.setFromTriplets(Constraint_eq_trip.begin(), Constraint_eq_trip.end());

	// Set sparse matrix for original inequality constraints
	Problem.Constraint.ie_orig_matrix = Eigen::SparseMatrix <double> (Problem.Variables_num + Problem.Constraints_ie_num, Problem.Variables_num);
	std::vector<Trip> Constraint_ie_orig_trip;
	Constraint_ie_orig_trip.reserve((Problem.Constraints_ie_num + 1) * Problem.Variables_num);
	for(int var_id = 0; var_id < Problem.Variables_num; ++ var_id){
		Constraint_ie_orig_trip.push_back(Trip(var_id, var_id, 1));
	}
	Problem.Constraint.ie_orig_matrix.setFromTriplets(Constraint_ie_orig_trip.begin(), Constraint_ie_orig_trip.end());
	
	// Set reference point for contraction
	Eigen::VectorXd feasible_point_1(Problem.Variables_num - Problem.Constraints_eq_num);
	Eigen::VectorXd feasible_point_2(Problem.Variables_num - Problem.Constraints_eq_num);
	Eigen::VectorXd feasible_point_3(Problem.Variables_num - Problem.Constraints_eq_num);
	feasible_point_1 << 0, 1, 0;
	feasible_point_2 << .2, 0, 0;
	feasible_point_3 << 0, 0, 1;
	Problem.Solution.reduced_vector_ref = 1. / 3. * (feasible_point_1 + feasible_point_2 + feasible_point_3);
}

int main(){
	LP_object* Problem = new LP_object;
	test_problem_1_set(*Problem);
	LP_process(*Problem);
	delete Problem;
}
