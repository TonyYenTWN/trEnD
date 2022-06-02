// Source file for linear programming using contraction and gradient projection algorithm
// Applicable when an initial feasible point in the interior is trivial

#include "../basic/LP_cpa.h"

// Normalize the equalitiy constraints
void LP_constraint_eq_normalization(LP_object &Problem){
	// Update both the constraints and the boundaries
	Problem.Constraint.eq_norm = Eigen::VectorXd(Problem.Constraints_eq_num);
	
	for(int item_iter = 0; item_iter < Problem.Constraints_eq_num; ++ item_iter){
		Problem.Constraint.eq_norm(item_iter) = Problem.Constraint.eq_matrix.row(item_iter).norm();
		Problem.Constraint.eq_matrix.row(item_iter) /= Problem.Constraint.eq_norm(item_iter);
		Problem.Boundary.eq_vector.row(item_iter) /= Problem.Constraint.eq_norm(item_iter);
	}
}

// Normalize the equalitiy boundaries
void LP_boundary_eq_normalization(LP_object &Problem){
	for(int item_iter = 0; item_iter < Problem.Constraints_eq_num; ++ item_iter){
		Problem.Boundary.eq_vector.row(item_iter) /= Problem.Constraint.eq_norm(item_iter);
	}
}

// Generate sparse matrix for covariance equality constraints
void LP_constraint_cov_eq_matrix_generation(LP_object &Problem){
	
	Problem.Solver.qr.compute(Problem.Constraint.eq_matrix.transpose());
	
	// Check if the equality constraints has maximum possible rank
	if(Problem.Solver.qr.rank() < Problem.Constraints_eq_num){
		// If redundant eqaulity constraint occurs, check feasibility of the system
		// Because an initial feasible point is needed for the algorithm, feasibility should be checked earlier
		// However, redundant eqaulity constraints can still occur so this step should still be performed
		
		// Update the information about equality constraints after removing the redundant ones
		Problem.Constraint.eq_matrix = Problem.Solver.qr.colsPermutation() * Problem.Constraint.eq_matrix;
		Problem.Constraint.eq_matrix = Problem.Constraint.eq_matrix.topRows(Problem.Solver.qr.rank());
		Problem.Constraints_eq_num = Problem.Solver.qr.rank();
		Problem.Boundary.eq_vector = Problem.Solver.qr.colsPermutation() * Problem.Boundary.eq_vector;
		Problem.Boundary.eq_vector = Problem.Boundary.eq_vector.head(Problem.Solver.qr.rank());
	}
	
	// Compute the default covariance matrix and its solver
	Problem.Constraint.eq_cov_matrix = Problem.Constraint.eq_matrix * Problem.Constraint.eq_matrix.transpose();
	Problem.Solver.ldlt.compute(Problem.Constraint.eq_cov_matrix);
}

// Compute the default gradient vector projected on the equality constraints
void LP_default_gradient(LP_object &Problem){
	Eigen::VectorXd Obj_cov_vector = Problem.Constraint.eq_matrix * Problem.Objective.vector;
	Problem.Proj_grad.vector_default = Problem.Objective.vector - Problem.Constraint.eq_matrix.transpose() * Problem.Solver.ldlt.solve(Obj_cov_vector);
	Problem.Constraint.ie_increment = Problem.Constraint.ie_matrix * Problem.Proj_grad.vector_default;
}

void LP_optimization(LP_object &Problem){
	// Parameters declaration
	double tol = pow(10, -12);
	double eps = pow(10, -8);
	
	// Main process
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
	Problem.Solution.vector = Problem.Solution.vector_ini;
	Problem.Objective.value = Problem.Solution.vector.dot(Problem.Objective.vector);
	
	// Main loop for searching the optimal solution
	while(1){
		// Update contracted solution and previous solution for later use in the loop
		Previous_obj_value = Problem.Objective.value;
		Previous_solution = Problem.Solution.vector;
		Contracted_solution = Problem.Solution.vector_ref;
		Contracted_solution += (1 - eps) * (Problem.Solution.vector - Problem.Solution.vector_ref);
		
		// Calculate the allowed increment along the objective direction before reaching a constraint
		Boundary_gap.col(0) = Problem.Constraint.ie_matrix * Contracted_solution - Problem.Boundary.ie_matrix.col(0);
		Boundary_gap.col(1) = Problem.Boundary.ie_matrix.col(1) - Problem.Constraint.ie_matrix * Contracted_solution;
		min_increment = std::numeric_limits<double>::infinity();
		for(int constraint_id = 0; constraint_id < Boundary_gap.rows(); ++ constraint_id){
			current_increment(0) = -Boundary_gap(constraint_id, 0) / Problem.Constraint.ie_increment(constraint_id);
			if(current_increment(0) > 0 && current_increment(0) < min_increment){
				min_increment = current_increment(0);
			}
			current_increment(1) = Boundary_gap(constraint_id, 1) / Problem.Constraint.ie_increment(constraint_id);
			if(current_increment(1) > 0 && current_increment(1) < min_increment){
				min_increment = current_increment(1);
			}
		}
		
		// Update the solution point and projected gradient on / near the active constraint(s)
		Problem.Solution.vector = Contracted_solution + min_increment * Problem.Proj_grad.vector_default;
		Projected_gradient = Problem.Solution.vector - Previous_solution;
		Projected_gradient /= Projected_gradient.norm();
		
		// Calculate the allowed increment along the projected gradient before reaching another set of constraint(s)
		Boundary_gap.col(0) = Problem.Constraint.ie_matrix * Problem.Solution.vector - Problem.Boundary.ie_matrix.col(0);
		Boundary_gap.col(1) = Problem.Boundary.ie_matrix.col(1) - Problem.Constraint.ie_matrix * Problem.Solution.vector;
		Projected_increment = Problem.Constraint.ie_matrix * Projected_gradient;
		min_increment = std::numeric_limits<double>::infinity();
		for(int constraint_id = 0; constraint_id < Boundary_gap.rows(); ++ constraint_id){
			if(abs(Projected_increment(constraint_id)) > tol){
				current_increment(0) = -Boundary_gap(constraint_id, 0) / Projected_increment(constraint_id);
				if(current_increment(0) > 0 && current_increment(0) < min_increment){
					min_increment = current_increment(0);
				}
				current_increment(1) = Boundary_gap(constraint_id, 1) / Projected_increment(constraint_id);
				if(current_increment(1) > 0 && current_increment(1) < min_increment){
					min_increment = current_increment(1);
				}			
			}			
		}
		Problem.Solution.vector += min_increment * Projected_gradient;
		Obj_value = Problem.Solution.vector.dot(Problem.Objective.vector);
		
		// Check if objective function can be improved or optimality has been reached
		if(Obj_value >= Previous_obj_value){
			Problem.Objective.value = Obj_value;
		}
		else{
			// If current objective value is smaller than previous, optimality occurs at previous solution
			Problem.Solution.vector = Previous_solution;
			break;
		}					
	}
}

// Function for the optimization algorithm
void LP_result_print(LP_object &Problem, std::string Problem_name){
	std::cout << "----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << "| " << Problem_name << std::endl;
	std::cout << "----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << "| Solution         | " << Problem.Solution.vector.transpose() << std::endl;
	std::cout << "| Objective value  | " << Problem.Objective.value << std::endl;
	std::cout << "----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << "\n" << std::endl;
}

void LP_process(LP_object &Problem, std::string Problem_name, bool print_result, bool constraint_flag, bool boundary_flag, bool objective_flag){
	// Pre-processing steps
	// Generate the covariance matrix for equality constraints if not yet done
	if(constraint_flag){
		LP_constraint_eq_normalization(Problem);
		LP_boundary_eq_normalization(Problem);
		LP_constraint_cov_eq_matrix_generation(Problem);	
		LP_default_gradient(Problem);
	}
	else{
		// Normalize boundaries if not yet done
		if(boundary_flag){
			LP_boundary_eq_normalization(Problem);
		}
		
		// Update the default gradient vector projected on the equality constraints if not yet done
		if(objective_flag){
			LP_default_gradient(Problem);
		}
	}
	
	// Optimization
	LP_optimization(Problem);
	
	// Output result
	if(print_result){
		LP_result_print(Problem, Problem_name);
	}
}