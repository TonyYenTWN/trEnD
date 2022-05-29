// Source file for linear programming using gradient projection algorithm
// Applicable when an initial feasible point is trivial
// Test case:
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

#include <iostream>
#include <omp.h>
#include "../basic/Basic_Definitions.h"
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>							// Use only when the matric is symmetric positive definite (ex. covariance or laplacian matrix)

typedef Eigen::Triplet <double> Trip;    				// Define a triplet object
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;

// Constraint object
struct LP_constraint{
	bool normalized = 0;
	Eigen::VectorXd norm;
	MatrixXb Status_ie_matrix;							// 1 = currently active; 0 = not yet called; -1 = previously called and currently inactive
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

// Solver object
struct LP_matrix_solver{
	bool Spd_flag = 0;
	bool Normal_flag = 0;
	Eigen::SimplicialLDLT <Eigen::SparseMatrix <double>> Spd; 							// Solver for a symmetric positive definite matrix
	Eigen::SparseLU <Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering <int>> Normal;	// Solver for an arbitrary matrix
};

struct LP_object{
	// Input parameters
	int Constraints_eq_num;
	int Constraints_ie_num;
	int Variables_num;
	Eigen::VectorXd Objective_vector;
	LP_constraint Constraint;
	LP_boundary Boundary;
	LP_matrix_solver Solver;
	
	// Process and output variables
	Eigen::VectorXd Proj_grad_vector;
	Eigen::VectorXd Solution_vector;
};

void LP_constraint_ie_reduced_generation(LP_object &Problem){
	// Warm up the solver if not yet done for the current matrix
	if(!Problem.Solver.Normal_flag){
		Problem.Solver.Normal.compute(Problem.Constraint.eq_orig_matrix.rightCols(Problem.Constraints_eq_num));
		Problem.Solver.Normal_flag = 1;
	}
	
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
		}	
	}
	
	double tol = pow(10, -6);
	Problem.Constraint.ie_reduced_matrix = reduced_matrix.sparseView(tol, 1);
}

void LP_boundary_ie_reduced_generation(LP_object &Problem){
	// Warm up the solver if not yet done for the current matrix
	if(!Problem.Solver.Normal_flag){
		Problem.Solver.Normal.compute(Problem.Constraint.eq_orig_matrix.rightCols(Problem.Constraints_eq_num));
		Problem.Solver.Normal_flag = 1;
	}
	
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

// Normalization of reduced inequality constraints
void LP_constraint_ie_reduced_normalization(LP_object &Problem, bool constraint_update = 1){
	Problem.Constraint.norm = Eigen::VectorXd::Ones(Problem.Variables_num + Problem.Constraints_ie_num);
	
	for(int row_id = Problem.Variables_num - Problem.Constraints_eq_num; row_id < Problem.Constraint.ie_reduced_matrix.rows(); ++ row_id){
		// Check if the constraint coefficients need to be normalized
		if(constraint_update){
			Problem.Constraint.norm(row_id) = Problem.Constraint.ie_reduced_matrix.row(row_id).norm();
			Problem.Constraint.ie_reduced_matrix.row(row_id) /= Problem.Constraint.norm(row_id);
		}
		// If not, only boundary values will be normalized
		Problem.Boundary.ie_reduced_matrix.row(row_id) /= Problem.Constraint.norm(row_id);
	}
	Problem.Constraint.normalized = 1;
	
	std::cout << Problem.Constraint.ie_reduced_matrix << "\n" << std::endl;
	std::cout << Problem.Boundary.ie_reduced_matrix << "\n" << std::endl;
	std::cout << Problem.Constraint.norm.transpose() << "\n" << std::endl;
}

void test_problem(){
	// Set dimension of the problem
	LP_object Problem;
	Problem.Constraints_eq_num = 2;
	Problem.Constraints_ie_num = 0;
	Problem.Variables_num = 4;
	
	// Set objective vector
	Problem.Objective_vector = Eigen::VectorXd(Problem.Variables_num);
	Problem.Objective_vector << 1, 2, 1, 1;
	
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
	Problem.Constraint.eq_orig_matrix = Eigen::SparseMatrix <double> (Problem.Constraints_eq_num, Problem.Variables_num);
	Problem.Constraint.eq_orig_matrix.setFromTriplets(Constraint_eq_trip.begin(), Constraint_eq_trip.end());
	
	// Set sparse matrix for original inequality constraints
	Problem.Constraint.ie_orig_matrix = Eigen::SparseMatrix <double> (Problem.Variables_num + Problem.Constraints_ie_num, Problem.Variables_num);
	std::vector<Trip> Constraint_ie_orig_trip;
	Constraint_ie_orig_trip.reserve((Problem.Constraints_eq_num + 1) * (Problem.Variables_num - Problem.Constraints_eq_num));
	for(int var_id = 0; var_id < Problem.Variables_num; ++ var_id){
		Constraint_ie_orig_trip.push_back(Trip(var_id, var_id, 1));
	}
	Problem.Constraint.ie_orig_matrix.setFromTriplets(Constraint_ie_orig_trip.begin(), Constraint_ie_orig_trip.end());
	// Put additional inequality constraints here
	
	// Generate sparse matrix for reduced inequality constraints
	Problem.Solver.Normal.compute(Problem.Constraint.eq_orig_matrix.rightCols(Problem.Constraints_eq_num));
	Problem.Solver.Normal_flag = 1;
	LP_constraint_ie_reduced_generation(Problem);
	LP_boundary_ie_reduced_generation(Problem);
	LP_constraint_ie_reduced_normalization(Problem);
}

int main(){
	test_problem();
}