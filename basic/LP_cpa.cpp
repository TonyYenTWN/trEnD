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
	bool eq_normalized_flag = 0;
	Eigen::VectorXd eq_vector;			// Value of original equality constraints
	Eigen::MatrixXd ie_matrix;			// Boundary of original inequality constraints; 0th column: lower bound value; 1st column: upper bound value
};

// Objective object
struct LP_objective{
	Eigen::VectorXd vector;
	Eigen::VectorXd cov_obj_vector;		// Dot product between the reduced objective and inequality constraints
	double value;
};

// Solver object
struct LP_cov_eq_solver{
	bool solved_flag = 0;
	Eigen::SimplicialLDLT <Eigen::SparseMatrix <double>> Spd; 	// Solver for a symmetric positive definite matrix
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
	Eigen::VectorXd Proj_grad_vector;
	LP_solution Solution;
};

void LP_constraint_eq_normalization(LP_object &Problem){
	// Check if equality constraints has been updated
	if(!Problem.Constraint.eq_normalized_flag){
		// Update both the constraints and the boundaries
		Problem.Constraint.eq_norm = Eigen::VectorXd(Problem.Constraints_eq_num);
		
		for(int item_iter = 0; item_iter < Problem.Constraints_eq_num; ++ item_iter){
			Problem.Constraint.eq_norm(item_iter) = Problem.Constraint.eq_matrix.row(item_iter).norm();
			Problem.Constraint.eq_matrix.row(item_iter) /= Problem.Constraint.eq_norm(item_iter);
			Problem.Boundary.ie_matrix.row(item_iter) /= Problem.Constraint.eq_norm(item_iter);
		}
		
		// Update boolean variables indicating normalization of the constraint and the boundary have been made for the current problem
		Problem.Constraint.eq_normalized_flag = 1;
		Problem.Boundary.eq_normalized_flag = 1;
	}
	else{
		// Update only the boundaries
		if(!Problem.Boundary.eq_normalized_flag){
			for(int item_iter = 0; item_iter < Problem.Constraints_eq_num; ++ item_iter){
				Problem.Boundary.ie_matrix.row(item_iter) /= Problem.Constraint.eq_norm(item_iter);
			}			
			
			// Update boolean variables indicating normalization of the boundary has been made for the current problem
			Problem.Boundary.eq_normalized_flag = 1;
		}
		else{
			std::cout << "The equality constraints and boundaries have already been normalized!" << std::endl;
		}
	}
}

void LP_constraint_cov_eq_matrix_generation(LP_object &Problem){
	if(!Problem.Constraint.eq_solved_flag){
		Problem.Constraint.eq_cov_matrix = Problem.Constraint.eq_matrix * Problem.Constraint.eq_matrix.transpose();
		Problem.Constraint.eq_cov_solver.compute(Problem.Constraint.eq_cov_matrix);
		Problem.Constraint.eq_solved_flag = 1;
		
		std::cout << Problem.Constraint.eq_cov_matrix << std::endl;
	}
	else{
		std::cout << "Covariance matrix of the equality constraints has already been generated!" << std::endl;
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
	Problem.Constraints_eq_num = 2;
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
	Eigen::VectorXd feasible_point_1(Problem.Variables_num - Problem.Constraints_eq_num);
	Eigen::VectorXd feasible_point_2(Problem.Variables_num - Problem.Constraints_eq_num);
	Eigen::VectorXd feasible_point_3(Problem.Variables_num - Problem.Constraints_eq_num);
	feasible_point_1 << 0, 0, 1, 0;
	feasible_point_2 << 0, .5, 0, .5;
	feasible_point_3 << .2, 0, 0, .6;	
	Problem.Solution.vector_ref = 1. / 3. * (feasible_point_1 + feasible_point_2 + feasible_point_3);
	
	// Generate sparse matrix for covariance equality constraints
	LP_constraint_cov_eq_matrix_generation(Problem);
}

int main(){
	//test_problem_1();
}