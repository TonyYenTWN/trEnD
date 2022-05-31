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
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>							// Use only when the matric is symmetric positive definite (ex. covariance or laplacian matrix)

typedef Eigen::Triplet <double> Trip;    				// Define a triplet object

// Constraint object
struct LP_constraint{
	bool normalized = 0;
	Eigen::VectorXd norm;
	Eigen::MatrixXi status_ie_reduced_matrix;			// 1 = currently active; 0 = not yet called; -1 = previously called and currently inactive
	Eigen::SparseMatrix <double> eq_orig_matrix;		// Coefficients for original equality constraints
	Eigen::SparseMatrix <double> ie_orig_matrix;		// Coefficients for original inequality constraints
	Eigen::SparseMatrix <double> ie_reduced_matrix;		// Coefficients for reduced inequality constraints			
	Eigen::SparseMatrix <double> cov_ie_reduced_matrix; // Dot product between the reduced inequality constraints
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
	Eigen::VectorXd cov_ie_obj_reduced_vector;	// Dot product between the reduced objective and inequality constraints
};

// Solver object
struct LP_matrix_solver{
	bool Normal_flag = 0;
	Eigen::SimplicialLDLT <Eigen::SparseMatrix <double>> Spd; 									// Solver for a symmetric positive definite matrix
	Eigen::SparseLU <Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering <int>> Normal;			// Solver for an arbitrary matrix
	Eigen::SparseLU <Eigen::SparseMatrix <double>, Eigen::COLAMDOrdering <int>> Normal_Trans;	// Solver for the transpose of an arbitrary matrix
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

// Linear system for reducing the redundant variables
void LP_constraint_redundant_matrix_solved(LP_object &Problem){
	Problem.Solver.Normal.compute(Problem.Constraint.eq_orig_matrix.rightCols(Problem.Constraints_eq_num));
	Problem.Solver.Normal_Trans.compute(Problem.Constraint.eq_orig_matrix.rightCols(Problem.Constraints_eq_num).transpose());
	Problem.Solver.Normal_flag = 1;
}

void LP_constraint_ie_reduced_generation(LP_object &Problem){
	// Warm up the solver if not yet done for the current matrix
	if(!Problem.Solver.Normal_flag){
		LP_constraint_redundant_matrix_solved(Problem);
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
		LP_constraint_redundant_matrix_solved(Problem);
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

void LP_objective_reduced_generation(LP_object &Problem){
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
}

// Covariance (dot product) of reduced inequality constraints
void LP_constraint_ie_reduced_covariance(LP_object &Problem){
	Problem.Constraint.cov_ie_reduced_matrix = Problem.Constraint.ie_reduced_matrix * Problem.Constraint.ie_reduced_matrix.transpose();
	Problem.Objective.cov_ie_obj_reduced_vector = Problem.Constraint.ie_reduced_matrix * Problem.Objective.reduced_vector;
}

// Function for the optimization algorithm
void LP_optimized(LP_object &Problem){
	// Parameters declaration
	double tol = pow(10, -8);
	double eps = pow(10, -6); 
	
	// Variables declaration before loop
	bool relax_flag = 0;
	int active_const_rank;
	double cov_value_temp;
	std::vector<Eigen::Vector2i> active_const_seq;
	std::vector<Eigen::Vector2i> active_const_seq_prev;
	std::vector<Eigen::Vector2i> active_const_reduced_seq;
	std::vector<Trip> Constraint_cov_sub_trip;
	Eigen::VectorXd Constraint_cov_D;
	Eigen::VectorXd Cov_ie_obj_reduced_sub;
	Eigen::SparseMatrix <double> Constraint_ie_reduced_sub;
	Eigen::SparseMatrix <double> Constraint_ie_reduced_covariance_sub;
	Eigen::SparseMatrix <double> Span_reduction;
	Eigen::MatrixXd Boundary_gap(Problem.Variables_num + Problem.Constraints_ie_num, 2); 
	Eigen::VectorXd Boundary_increment(Problem.Variables_num + Problem.Constraints_ie_num);
	
	// Initialize the starting point and status of activeness of each constraints
	Problem.Constraint.status_ie_reduced_matrix = Eigen::MatrixXi::Zero(Problem.Variables_num + Problem.Constraints_ie_num, 2);
	Problem.Solution.reduced_vector = Eigen::VectorXd::Zero(Problem.Variables_num - Problem.Constraints_eq_num);
	
	int count = 0;
	while(count < 2){
		active_const_reduced_seq.clear();
		active_const_seq.clear();
		count += 1;
		
		if(!relax_flag){
			// Check whether any of the constraints are active at current position
			Boundary_gap.col(0) = Problem.Constraint.ie_reduced_matrix * Problem.Solution.reduced_vector - Problem.Boundary.ie_reduced_matrix.col(0);
			Boundary_gap.col(1) = Problem.Boundary.ie_reduced_matrix.col(1) - Problem.Constraint.ie_reduced_matrix * Problem.Solution.reduced_vector;
			active_const_seq.reserve(Problem.Variables_num + Problem.Constraints_ie_num);
			for(int row_id = 0; row_id < Boundary_gap.rows(); ++ row_id){
				if(Boundary_gap(row_id, 0) < tol){
					Problem.Constraint.status_ie_reduced_matrix(row_id, 0) = 1;
					active_const_seq.push_back(Eigen::Vector2i(row_id, 0));
				}
				else if(Boundary_gap(row_id, 1) < tol){
					Problem.Constraint.status_ie_reduced_matrix(row_id, 1) = 1;
					active_const_seq.push_back(Eigen::Vector2i(row_id, 1));
				}
			}			
		}	
		
		// Check if there are active constraints
		if(active_const_seq.size() > 0){
			// Construct sub_covariance_matrix of all the active constraints, should the set be not empty
			Constraint_cov_sub_trip.reserve(pow(active_const_seq.size(), 2));
			Constraint_ie_reduced_covariance_sub = Eigen::SparseMatrix <double>(active_const_seq.size(), active_const_seq.size());
			Span_reduction = Eigen::SparseMatrix <double>(active_const_seq.size(), Problem.Variables_num + Problem.Constraints_ie_num);
			for(int row_id = 0; row_id < active_const_seq.size(); ++ row_id){
				Constraint_cov_sub_trip.push_back(Trip(active_const_seq[row_id](0), active_const_seq[row_id](0), 1));
			}
			Span_reduction.setFromTriplets(Constraint_cov_sub_trip.begin(), Constraint_cov_sub_trip.end());
			Constraint_ie_reduced_covariance_sub = Span_reduction * Problem.Constraint.cov_ie_reduced_matrix * Span_reduction.transpose();
			Constraint_cov_sub_trip.clear();
			Problem.Solver.Spd.compute(Constraint_ie_reduced_covariance_sub);
			
			// Check if some of the active constraints are redundant
			if(abs(Problem.Solver.Spd.determinant()) < tol){
				active_const_rank = active_const_seq.size();
				active_const_reduced_seq.reserve(active_const_rank);	
				Constraint_cov_D = Problem.Solver.Spd.vectorD();
				for(int row_id = 0; row_id < active_const_seq.size(); ++ row_id){
					if(abs(Constraint_cov_D(row_id)) < tol){
						active_const_rank -= 1;
					}
					else{
						active_const_reduced_seq.push_back(active_const_seq[row_id]);
					}
				}
				
				// Check if the active constraints form a degenerate extreme point
				if(active_const_rank == Problem.Variables_num - Problem.Constraints_eq_num){
					// Relaxation of the constraints; all active constraints thus become inactive
					for(int row_id = 0; row_id < active_const_seq.size(); ++ row_id){
						if(active_const_seq[row_id](1) == 0){
							Problem.Boundary.ie_reduced_matrix(active_const_seq[row_id](0), 0) -= eps;
							Problem.Constraint.status_ie_reduced_matrix(active_const_seq[row_id](0), 0) = 0;
						}
						else{
							Problem.Boundary.ie_reduced_matrix(active_const_seq[row_id](0), 1) += eps;
							Problem.Constraint.status_ie_reduced_matrix(active_const_seq[row_id](0), 1) = 0;
						}
					}
					
					// Exit the current loop; now lies in the interior and update variables accordingly
					relax_flag = 1;
					active_const_seq_prev.clear();
					continue;
				}
				// If not, find the span of the normal vectors of the active constraints
				else{
					Constraint_cov_sub_trip.reserve(pow(active_const_rank, 2));
					Constraint_ie_reduced_covariance_sub = Eigen::SparseMatrix <double>(active_const_rank, active_const_rank);
					Span_reduction = Eigen::SparseMatrix <double>(active_const_rank, Problem.Variables_num + Problem.Constraints_ie_num);
					for(int row_id = 0; row_id < active_const_rank; ++ row_id){
						Constraint_cov_sub_trip.push_back(Trip(active_const_reduced_seq[row_id](0), active_const_reduced_seq[row_id](0), 1));
					}
					Span_reduction.setFromTriplets(Constraint_cov_sub_trip.begin(), Constraint_cov_sub_trip.end());
					Constraint_ie_reduced_covariance_sub = Span_reduction * Problem.Constraint.cov_ie_reduced_matrix * Span_reduction.transpose();					
					Constraint_cov_sub_trip.clear();
					Problem.Solver.Spd.compute(Constraint_ie_reduced_covariance_sub);			
				}				
			}
			Cov_ie_obj_reduced_sub = Span_reduction * Problem.Objective.cov_ie_obj_reduced_vector;
			Constraint_ie_reduced_sub = Span_reduction * Problem.Constraint.ie_reduced_matrix;
			Problem.Proj_grad_vector = Problem.Objective.reduced_vector;
			Problem.Proj_grad_vector -= Constraint_ie_reduced_sub.transpose() * Problem.Solver.Spd.solve(Cov_ie_obj_reduced_sub);					
		}
		// If there are no active constraints, projected gradient is the objective vector
		else{
			Problem.Proj_grad_vector = Problem.Objective.reduced_vector;
		}
		
		// Check maximum allowed increment along the projected gradient before reaching a constraint
		Boundary_increment = Problem.Constraint.ie_reduced_matrix * Problem.Proj_grad_vector;
		
		// Update status of constraints
		if(active_const_seq_prev.size() > 0){
			for(int row_id = 0; row_id < active_const_seq_prev.size(); ++ row_id){
				if(Problem.Constraint.status_ie_reduced_matrix(active_const_seq_prev[row_id](0), active_const_seq_prev[row_id](1)) == 0){
					Problem.Constraint.status_ie_reduced_matrix(active_const_seq_prev[row_id](0), active_const_seq_prev[row_id](1)) = -1;
				}
			}			
		}
		
		// Update the list of previous active constraints for the next loop
		active_const_seq_prev = active_const_seq;
		
		// Code execution reaches here only if relaxation of constraints does not occur at this loop 
		if(relax_flag){
			relax_flag = 0;
		}
		std::cout << Problem.Boundary.ie_reduced_matrix << "\n" << std::endl;		
	}

	//std::cout << Boundary_increment.transpose() << "\n" << std::endl;
	//std::cout << Problem.Constraint.status_ie_reduced_matrix << "\n" << std::endl;
	//std::cout << Constraint_ie_reduced_covariance_sub << "\n" << std::endl;
	//std::cout << Problem.Constraint.cov_ie_reduced_matrix << "\n" << std::endl;
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
	Constraint_ie_orig_trip.reserve((Problem.Constraints_eq_num + 1) * (Problem.Variables_num - Problem.Constraints_eq_num));
	for(int var_id = 0; var_id < Problem.Variables_num; ++ var_id){
		Constraint_ie_orig_trip.push_back(Trip(var_id, var_id, 1));
	}
	Problem.Constraint.ie_orig_matrix.setFromTriplets(Constraint_ie_orig_trip.begin(), Constraint_ie_orig_trip.end());
	// Put additional inequality constraints here
	
	// Generate sparse matrix for reduced inequality constraints
	LP_constraint_redundant_matrix_solved(Problem);
	LP_constraint_ie_reduced_generation(Problem);
	LP_boundary_ie_reduced_generation(Problem);
	LP_objective_reduced_generation(Problem);
	LP_constraint_ie_reduced_normalization(Problem);
	LP_constraint_ie_reduced_covariance(Problem);
	
	// Solve the LP
	LP_optimized(Problem);
}

int main(){
	test_problem();
}
