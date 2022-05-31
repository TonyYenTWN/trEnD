// Source file for linear programming using contraction and gradient projection algorithm
// Applicable when an initial feasible point in the interior is trivial
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
#include "../basic/Eigen_Sparse.h"

// Constraint object
struct LP_constraint{
	bool normalized = 0;
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
	double value;
};

// Solver object
struct LP_matrix_solver{
	bool Normal_flag = 0;
	//Eigen::SimplicialLDLT <Eigen::SparseMatrix <double>> Spd; 								// Solver for a symmetric positive definite matrix
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

// Linear system for reducing the redundant variables
void LP_constraint_redundant_matrix_solved(LP_object &Problem){
	Problem.Solver.Normal.compute(Problem.Constraint.eq_orig_matrix.rightCols(Problem.Constraints_eq_num));
	Problem.Solver.Normal_Trans.compute(Problem.Constraint.eq_orig_matrix.rightCols(Problem.Constraints_eq_num).transpose());
	Problem.Solver.Normal_flag = 1;
}

void LP_constraint_ie_reduced_generation(LP_object &Problem){
	// Set default value of reduced constraint matrix
	Eigen::MatrixXd reduced_matrix = Eigen::MatrixXd(Problem.Constraint.ie_orig_matrix.leftCols(Problem.Variables_num - Problem.Constraints_eq_num));
	
	// Check if there are eqaulity constraints in the original problem
	if(Problem.Constraints_eq_num != 0){
		// Warm up the solver if not yet done for the current matrix
		if(!Problem.Solver.Normal_flag){
			LP_constraint_redundant_matrix_solved(Problem);
		}
		
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
	
	double tol = pow(10, -8);
	Problem.Constraint.ie_reduced_matrix = reduced_matrix.sparseView(tol, 1);
}

void LP_boundary_ie_reduced_generation(LP_object &Problem){
	// Set default value of reduced boundary vector
	Problem.Boundary.ie_reduced_matrix = Problem.Boundary.ie_orig_matrix;
	
	// Check if there are eqaulity constraints in the original problem
	if(Problem.Constraints_eq_num != 0){
		// Warm up the solver if not yet done for the current matrix
		if(!Problem.Solver.Normal_flag){
			LP_constraint_redundant_matrix_solved(Problem);
		}

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
		// Warm up the solver if not yet done for the current matrix
		if(!Problem.Solver.Normal_flag){
			LP_constraint_redundant_matrix_solved(Problem);
		}

		// Z = c_e^T * x_e + c_r^T * x_r 
		//   = c_e^T * x_e - c_r^T * A_r^(-1) * (A_e * x_e - b)
		//   = (c_e^T - c_r^T * A_r^(-1) * A_e) * x_e + c_r^T * A_r^(-1) * b
		Problem.Objective.reduced_vector = Problem.Objective.orig_vector.head(Problem.Variables_num - Problem.Constraints_eq_num);
		Eigen::VectorXd coeff_redundant = Problem.Objective.orig_vector.tail(Problem.Constraints_eq_num);
		Problem.Objective.reduced_vector -= (Problem.Constraint.eq_orig_matrix.leftCols(Problem.Variables_num - Problem.Constraints_eq_num)).transpose() * Problem.Solver.Normal_Trans.solve(coeff_redundant);		
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
	
	if(constraint_update){
		Problem.Constraint.ie_reduced_increment = Problem.Constraint.ie_reduced_matrix * Problem.Objective.reduced_vector;
	}
	Problem.Constraint.normalized = 1;
}

// Function for the optimization algorithm
void LP_optimized(LP_object &Problem){
	// Parameters declaration
	double tol = pow(10, -8);
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
	Problem.Objective.value = Problem.Solution.reduced_vector.dot(Problem.Objective.reduced_vector);
	
	// Main loop for searching the optimal solution
	while(1){
		// Update contracted solution and previous solution for later use in the loop
		Previous_obj_value = Problem.Objective.value;
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
			Problem.Objective.value = Obj_value;
		}
		else{
			break;
		}		
	}
	
	// Reconstruct the original solution and objective value from the reduced problem, if necessary
	if(Problem.Constraints_eq_num != 0){
		// x_r = -A_r^(-1) * (A_e * x_e - b)
		//Problem.Constraint.eq_orig_matrix
		Eigen::VectorXd redundant_variables = -Problem.Solver.Normal.solve(Problem.Solution.reduced_vector - Problem.Boundary.eq_vector);
	}
	else{
		Problem.Solution.orig_vector = Problem.Solution.reduced_vector;
	}
	
	std::cout << "Solution to the reduced problem:" << std::endl;
	std::cout << Problem.Solution.reduced_vector.transpose() << "\n" << std::endl;
	std::cout << "Objective value of the reduced problem:" << std::endl;
	std::cout << Problem.Objective.value << "\n" << std::endl;
}

void test_problem(){
	// Set dimension of the problem
	LP_object Problem;
	Problem.Constraints_eq_num = 2;
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
	Problem.Constraint.eq_orig_matrix = Eigen::SparseMatrix <double> (Problem.Constraints_eq_num, Problem.Variables_num);
	Problem.Constraint.eq_orig_matrix.setFromTriplets(Constraint_eq_trip.begin(), Constraint_eq_trip.end());
	
	// Set sparse matrix for original inequality constraints
	Problem.Constraint.ie_orig_matrix = Eigen::SparseMatrix <double> (Problem.Variables_num + Problem.Constraints_ie_num, Problem.Variables_num);
	std::vector<Trip> Constraint_ie_orig_trip;
	Constraint_ie_orig_trip.reserve((Problem.Constraints_eq_num + 1) * (Problem.Variables_num - Problem.Constraints_eq_num) + Problem.Constraints_ie_num * Problem.Variables_num);
	for(int var_id = 0; var_id < Problem.Variables_num; ++ var_id){
		Constraint_ie_orig_trip.push_back(Trip(var_id, var_id, 1));
	}
	Problem.Constraint.ie_orig_matrix.setFromTriplets(Constraint_ie_orig_trip.begin(), Constraint_ie_orig_trip.end());
	// Put additional inequality constraints here
	
	// Set reference point for contraction
	Problem.Solution.reduced_vector_ref = Eigen::VectorXd(Problem.Variables_num - Problem.Constraints_eq_num);
	Problem.Solution.reduced_vector_ref << 1. / 6., 1. / 15.;
	
	// Generate sparse matrix for reduced inequality constraints
	LP_constraint_redundant_matrix_solved(Problem);
	LP_constraint_ie_reduced_generation(Problem);
	LP_boundary_ie_reduced_generation(Problem);
	LP_objective_reduced_generation(Problem);
	LP_constraint_ie_reduced_normalization(Problem);
	
	// Solve the LP
	LP_optimized(Problem);
}

int main(){
	test_problem();
}
