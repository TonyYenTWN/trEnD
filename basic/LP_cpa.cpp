// Source file for linear programming using contraction and gradient projection algorithm
// Applicable when an initial feasible point in the interior is trivial

#include <iostream>
#include <omp.h>
#include "../basic/Basic_Definitions.h"
#include "../basic/Eigen_Sparse.h"

// Constraint object
struct LP_constraint{
	bool eq_normalized_flag = 0;
	bool eq_solved_flag = 0;
	Eigen::VectorXd eq_norm;
	Eigen::VectorXd ie_increment;										// Increment per unit variation in the gradient direction
	Eigen::SparseMatrix <double> eq_matrix;								// Coefficients for equality constraints
	Eigen::SparseMatrix <double> ie_matrix;								// Coefficients for inequality constraints
	Eigen::SparseMatrix <double> eq_cov_matrix;							// Dot product of the equality constraints
	Eigen::SimplicialLDLT <Eigen::SparseMatrix <double>> eq_cov_solver; // Solver for the covariance (dot product) matrix
};

// Boundary object
struct LP_boundary{
	bool eq_normalized_flag = 0;		// Whenever the boundary of equality constraints is updated, this boolean variable should be set to 0
	Eigen::VectorXd eq_vector;			// Value of original equality constraints
	Eigen::MatrixXd ie_matrix;			// Boundary of original inequality constraints; 0th column: lower bound value; 1st column: upper bound value
};

// Objective object
struct LP_objective{
	Eigen::VectorXd vector;
	double value;
};

// Solver object
struct LP_cov_eq_solver{
	bool solved_flag = 0;										// Whenever the equality constraints are updated, this boolean variable should be set to 0
	Eigen::SimplicialLDLT <Eigen::SparseMatrix <double>> Spd; 	// Solver for a symmetric positive definite matrix
};

// Projection gradient object
struct LP_proj_gradient{
	bool uodated_flag = 0;				// Whenever the objective function is updated, this boolean variable should be set to 0
	Eigen::VectorXd vector_default; 	// The projection of objective vector onto the equality constraints
	Eigen::VectorXd vector_current;
};

// Solution object
struct LP_solution{
	Eigen::VectorXd vector;
	Eigen::VectorXd vector_ref;
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
	
	// Process and output variables
	LP_solution Solution;
};

void LP_constraint_eq_normalization(LP_object &Problem){
	// Check if equality constraints has been normalized
	if(!Problem.Constraint.eq_normalized_flag){
		// Update both the constraints and the boundaries
		Problem.Constraint.eq_norm = Eigen::VectorXd(Problem.Constraints_eq_num);
		
		for(int item_iter = 0; item_iter < Problem.Constraints_eq_num; ++ item_iter){
			Problem.Constraint.eq_norm(item_iter) = Problem.Constraint.eq_matrix.row(item_iter).norm();
			Problem.Constraint.eq_matrix.row(item_iter) /= Problem.Constraint.eq_norm(item_iter);
			Problem.Boundary.eq_vector.row(item_iter) /= Problem.Constraint.eq_norm(item_iter);
		}
		
		// Update boolean variables indicating normalization of the constraint and the boundary have been made for the current problem
		Problem.Constraint.eq_normalized_flag = 1;
		Problem.Boundary.eq_normalized_flag = 1;
	}
	else{
		// Normalize only the boundaries if not yet done
		if(!Problem.Boundary.eq_normalized_flag){
			for(int item_iter = 0; item_iter < Problem.Constraints_eq_num; ++ item_iter){
				Problem.Boundary.eq_vector.row(item_iter) /= Problem.Constraint.eq_norm(item_iter);
			}			
			
			// Update boolean variables indicating normalization of the boundary has been made for the current problem
			Problem.Boundary.eq_normalized_flag = 1;
		}
		else{
			std::cout << "WARNING: The equality constraints and boundaries have already been normalized!" << std::endl;
		}
	}
}

// Generate sparse matrix for covariance equality constraints
void LP_constraint_cov_eq_matrix_generation(LP_object &Problem){
	// Check if equality constraints has been normalized; normalize if not
	if(!Problem.Constraint.eq_normalized_flag){
		LP_constraint_eq_normalization(Problem);
	}
	
	// Check if covariance matrix has been generated
	if(!Problem.Constraint.eq_solved_flag){
		// Set precision of zero detection
		double tol = pow(10, -12);
		
		// Compute the default covariance matrix and its solver
		Problem.Constraint.eq_cov_matrix = Problem.Constraint.eq_matrix * Problem.Constraint.eq_matrix.transpose();
		Problem.Constraint.eq_cov_solver.compute(Problem.Constraint.eq_cov_matrix);
		
		// Check if the coavriance matrix is not of full rank (determinant = 0)
		if(abs(Problem.Constraint.eq_cov_solver.determinant()) < tol){
			// Declare a dynamic vector of list of active (non-redeundant) equality constraints
			std::vector<int> active_const_seq;
			active_const_seq.reserve(Problem.Constraints_eq_num);
			Eigen::VectorXd Constraint_eq_cov_D = Problem.Constraint.eq_cov_solver.vectorD();
			for(int item_iter = 0; item_iter < Problem.Constraints_eq_num; ++ item_iter){
				if(abs(Constraint_eq_cov_D(item_iter)) >= tol){
					active_const_seq.push_back(item_iter);
				}
			}
			
			std::vector<Trip> Constraint_cov_sub_trip;
			Constraint_cov_sub_trip.reserve(active_const_seq.size());
			Eigen::SparseMatrix <double> Span_reduction(active_const_seq.size(), Problem.Constraints_eq_num);
			for(int item_iter = 0; item_iter < active_const_seq.size(); ++ item_iter){
				Constraint_cov_sub_trip.push_back(Trip(active_const_seq[item_iter], active_const_seq[item_iter], 1));
			}
			Span_reduction.setFromTriplets(Constraint_cov_sub_trip.begin(), Constraint_cov_sub_trip.end());
			
			// Update the information about equality constraints after removing the redundant ones
			Problem.Constraint.eq_cov_matrix = Span_reduction * Problem.Constraint.eq_cov_matrix * Span_reduction.transpose();
			Problem.Constraint.eq_matrix = Span_reduction * Problem.Constraint.eq_matrix;
			Problem.Constraints_eq_num = active_const_seq.size();
		}
		
		// Update boolean variables indicating covariance matrix has been generated for the current problem
		Problem.Constraint.eq_solved_flag = 1;
	}
	else{
		std::cout << "WARNING: Covariance matrix of the equality constraints has already been generated!" << std::endl;
	}
}

// Compute the default gradient vector projected on the equality constraints
void LP_default_gradient(LP_object &Problem){
	// Check if covariance matrix has been generated
	if(!Problem.Constraint.eq_solved_flag){
		LP_constraint_cov_eq_matrix_generation(Problem);
	}	
	
	Eigen::VectorXd Obj_cov_vector = Problem.Constraint.eq_matrix * Problem.Objective.vector;
	Problem.Proj_grad_vector_default = Problem.Objective.vector - Problem.Constraint.eq_cov_solver.solve(Obj_cov_vector);
}

void LP_optimization(LP_object &Problem, bool print_result = 1){
	// Parameters declaration
	double tol = pow(10, -12);
	double eps = pow(10, -6);
	
	// Pre-processing steps
	// Generate the covariance matrix for equality constraints if not yet done
	if(!Problem.Constraint.eq_solved_flag){
		LP_constraint_cov_eq_matrix_generation(Problem);	// For test of overwriting prevention, run the function twice
		LP_default_gradient(Problem);
	}
	else{
		// Normalize boundaries if not yet done
		if(!Problem.Boundary.eq_normalized_flag){
			LP_constraint_eq_normalization(Problem);
		}
	}
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
// A trivial initial feasible point is (0, 0, 1, 0), (0, .5, 0, .5) or (.2, 0, 0, .6)
// A feasible interior point for reference of contraction can be found by averaging the 3 trivial points
void test_problem_1(){
	// Set dimension of the problem
	LP_object Problem;
	Problem.Constraints_eq_num = 3;		// Repeat the second equality constraint twice to test for redundancy
	Problem.Constraints_ie_num = 0;
	Problem.Variables_num = 4;

	// Set objective vector
	Problem.Objective.vector = Eigen::VectorXd(Problem.Variables_num);
	Problem.Objective.vector << 1, 2, 1, 1;
	
	// Set boudary values for original equality and inequality constraints
	Problem.Boundary.eq_vector = Eigen::VectorXd(Problem.Constraints_eq_num);
	Problem.Boundary.eq_vector << 2, 1;
	Problem.Boundary.ie_matrix = Eigen::MatrixXd(Problem.Variables_num + Problem.Constraints_ie_num, 2);
	Problem.Boundary.ie_matrix.col(0) = Eigen::VectorXd::Zero(Problem.Variables_num + Problem.Constraints_ie_num);
	Problem.Boundary.ie_matrix.col(1) = Eigen::VectorXd::Ones(Problem.Variables_num + Problem.Constraints_ie_num);
	
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
	Problem.Constraint.eq_matrix = Eigen::SparseMatrix <double> (Problem.Constraints_eq_num, Problem.Variables_num);
	Problem.Constraint.eq_matrix.setFromTriplets(Constraint_eq_trip.begin(), Constraint_eq_trip.end());
	
	// Set sparse matrix for original inequality constraints
	Problem.Constraint.ie_matrix = Eigen::SparseMatrix <double> (Problem.Variables_num + Problem.Constraints_ie_num, Problem.Variables_num);
	std::vector<Trip> Constraint_ie_trip;
	Constraint_ie_trip.reserve((Problem.Constraints_ie_num + 1) * Problem.Variables_num);
	for(int var_id = 0; var_id < Problem.Variables_num; ++ var_id){
		Constraint_ie_trip.push_back(Trip(var_id, var_id, 1));
	}
	Problem.Constraint.ie_matrix.setFromTriplets(Constraint_ie_trip.begin(), Constraint_ie_trip.end());
	// Put additional inequality constraints here
	
	// Set reference point for contraction
	Eigen::VectorXd feasible_point_1(Problem.Variables_num);
	Eigen::VectorXd feasible_point_2(Problem.Variables_num);
	Eigen::VectorXd feasible_point_3(Problem.Variables_num);
	feasible_point_1 << 0, 0, 1, 0;
	feasible_point_2 << 0, .5, 0, .5;
	feasible_point_3 << .2, 0, 0, .6;	
	Problem.Solution.vector_ref = 1. / 3. * (feasible_point_1 + feasible_point_2 + feasible_point_3);
	
	// Optimization of the Problem
	LP_optimization(Problem);
}

int main(){
	test_problem_1();
}