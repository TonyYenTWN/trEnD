// Source file for linear programming using gradient projection algorithm
// Applicable when an initial feasible solution is trivial

#include "../basic/LP_gpa.h"

// Delete redundant equality constraints
void LP_constraint_eq_redundant_deletion(LP_object &Problem){
	// Because an initial feasible point is needed for the algorithm, feasibility should be checked earlier
	// However, redundant eqaulity constraints can still occur so this step should still be performed
	
	Problem.Solver.qr.compute(Problem.Constraint.eq_orig_matrix.transpose());
	
	// Check if the equality constraints has maximum possible rank
	if(Problem.Solver.qr.rank() < Problem.Constraints_eq_num){
		// Update the information about equality constraints after removing the redundant ones
		Problem.Constraint.eq_orig_matrix = Problem.Solver.qr.colsPermutation() * Problem.Constraint.eq_orig_matrix;
		Problem.Constraint.eq_orig_matrix = Problem.Constraint.eq_orig_matrix.topRows(Problem.Solver.qr.rank());
		Problem.Constraints_eq_num = Problem.Solver.qr.rank();
		Problem.Boundary.eq_vector = Problem.Solver.qr.colsPermutation() * Problem.Boundary.eq_vector;
		Problem.Boundary.eq_vector = Problem.Boundary.eq_vector.head(Problem.Solver.qr.rank());
	}
}

// Choose redundant variables forming full rank square matrix from the equality constraints
void LP_variables_permutation(LP_object &Problem, bool stepwise_obj){
	// Set precision of zero detection
	double tol = pow(10, -12);
	
	std::vector<Trip> Permutation_trip;
	Permutation_trip.reserve(Problem.Variables_num);	
	int redundant_size = 0;
	int remaining_size = 0;
	int start_id;
	Eigen::SparseMatrix <double> Redundant_matrix;
	Eigen::SparseMatrix <double> Subspan_matrix(Problem.Variables_num, Problem.Constraints_eq_num);
	Subspan_matrix.reserve(Eigen::VectorXi::Constant(Problem.Constraints_eq_num, 2));
	Problem.Constraint.permutation_matrix = Eigen::SparseMatrix <double> (Problem.Variables_num, Problem.Variables_num);
	
	// The initial column should be non-zero
	for(int var_iter = 0; var_iter < Problem.Variables_num; ++ var_iter){
		if(Problem.Constraint.eq_orig_matrix.col(var_iter).norm() > tol){
			Subspan_matrix.insert(var_iter, redundant_size) = 1;
			Permutation_trip.push_back(Trip(var_iter, Problem.Variables_num - Problem.Constraints_eq_num, 1));
			start_id = var_iter + 1;
			redundant_size += 1;
			break;
		}
		else{
			Permutation_trip.push_back(Trip(var_iter, remaining_size, 1));
			remaining_size += 1;	
		}
	}
	for(int var_iter = start_id; var_iter < Problem.Variables_num; ++ var_iter){
		if(redundant_size < Problem.Constraints_eq_num){
			Subspan_matrix.insert(var_iter, redundant_size) = 1;
			Redundant_matrix = Problem.Constraint.eq_orig_matrix * Subspan_matrix;
			Problem.Solver.qr.compute(Redundant_matrix.leftCols(redundant_size + 1));
			
			// Check if the new column is linearly independent to other ones
			if(Problem.Solver.qr.rank() == redundant_size + 1){
				Permutation_trip.push_back(Trip(var_iter, Problem.Variables_num - Problem.Constraints_eq_num + redundant_size, 1));			
				redundant_size += 1;
			}
			else{
				// If not indepedent, remove the current column
				Subspan_matrix.coeffRef(var_iter, redundant_size) = 0;
				Permutation_trip.push_back(Trip(var_iter, remaining_size, 1));
				remaining_size += 1;
			}
		}
		else{
			Permutation_trip.push_back(Trip(var_iter, remaining_size, 1));
			remaining_size += 1;
		}
	}
	Problem.Constraint.permutation_matrix.setFromTriplets(Permutation_trip.begin(), Permutation_trip.end());
	
	// Update the constraints, objective function, and initial feasible solution with the permutation
	Problem.Constraint.eq_orig_matrix = Problem.Constraint.eq_orig_matrix * Problem.Constraint.permutation_matrix;
	Problem.Boundary.ie_orig_matrix = Problem.Constraint.permutation_matrix.transpose() * Problem.Boundary.ie_orig_matrix;
	Problem.Objective.orig_vector = Problem.Constraint.permutation_matrix.transpose() * Problem.Objective.orig_vector;
	if(stepwise_obj){
		Problem.Objective.varying_vector = Problem.Constraint.permutation_matrix.transpose() * Problem.Objective.varying_vector;	
	}
}

// Permute solution to correct order
void LP_solution_permutation(LP_object &Problem, bool inverse){
	if(!inverse){
		Problem.Solution.orig_vector = Problem.Constraint.permutation_matrix.transpose() * Problem.Solution.orig_vector;
	}
	else{
		Problem.Solution.orig_vector = Problem.Constraint.permutation_matrix * Problem.Solution.orig_vector;
	}
}

// Linear system for reducing the redundant variables
void LP_constraint_redundant_matrix_solver(LP_object &Problem){
	Problem.Solver.lu.compute(Problem.Constraint.eq_orig_matrix.rightCols(Problem.Constraints_eq_num));
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
		// Do not use openmp parallel thread to run the loop below; may yield wrong solution
		for(int var_id = 0; var_id < Problem.Variables_num - Problem.Constraints_eq_num; ++ var_id){
			col_orig = Problem.Constraint.eq_orig_matrix.col(var_id);
			reduced_matrix.middleRows(Problem.Variables_num - Problem.Constraints_eq_num, Problem.Constraints_eq_num).col(var_id) = -Problem.Solver.lu.solve(col_orig);
		}
						
		// Check if there are ineqaulity constraints in the original problem
		if(Problem.Constraints_ie_num != 0){
			// If inequality constraints exist in the original, renew the coefficients
			// l_i <= M_e * x_e + M_r * x_r <= u_i
			// l_i <= M_e * x_e - M_r * A_r^(-1) * (A_e * x_e - b) <= u_i
			// l_i - M_r * A_r^(-1) * b <= (M_e - M_r * A_r^(-1) * A_e) * x_e <= u_i - M_r * A_r^(-1) * b			
		}
	}
	
	double tol = pow(10, -12);
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
		Problem.Boundary.ie_reduced_matrix.middleRows(Problem.Variables_num - Problem.Constraints_eq_num, Problem.Constraints_eq_num).col(0) -= Problem.Solver.lu.solve(Problem.Boundary.eq_vector);
		Problem.Boundary.ie_reduced_matrix.middleRows(Problem.Variables_num - Problem.Constraints_eq_num, Problem.Constraints_eq_num).col(1) -= Problem.Solver.lu.solve(Problem.Boundary.eq_vector);
						
		// Check if there are ineqaulity constraints in the original problem
		if(Problem.Constraints_ie_num != 0){
			// If inequality constraints exist in the original, renew the coefficients
			// l_i <= M_e * x_e + M_r * x_r <= u_i
			// l_i <= M_e * x_e - M_r * A_r^(-1) * (A_e * x_e - b) <= u_i
			// l_i - M_r * A_r^(-1) * b <= (M_e - M_r * A_r^(-1) * A_e) * x_e <= u_i - M_r * A_r^(-1) * b
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
		Problem.Objective.reduced_vector -= (Problem.Constraint.eq_orig_matrix.leftCols(Problem.Variables_num - Problem.Constraints_eq_num)).transpose() * Problem.Solver.lu.transpose().solve(coeff_redundant);		
	}
}

void LP_feasible_solution_reduced_generation(LP_object &Problem){
	Problem.Solution.reduced_vector = Problem.Solution.orig_vector.head(Problem.Variables_num - Problem.Constraints_eq_num);
}

// Normalization of reduced inequality constraints
void LP_constraint_ie_reduced_normalization(LP_object &Problem){
	Problem.Constraint.norm = Eigen::VectorXd::Ones(Problem.Variables_num + Problem.Constraints_ie_num);

	#pragma omp parallel
	{
		#pragma omp for
		for(int row_id = Problem.Variables_num - Problem.Constraints_eq_num; row_id < Problem.Constraint.ie_reduced_matrix.rows(); ++ row_id){
			Problem.Constraint.norm(row_id) = Problem.Constraint.ie_reduced_matrix.row(row_id).norm();
			Problem.Constraint.ie_reduced_matrix.row(row_id) /= Problem.Constraint.norm(row_id);
		}
	}
}

// Normalization of reduced inequality boundaries
void LP_boundary_ie_reduced_normalization(LP_object &Problem){
	#pragma omp parallel
	{
		#pragma omp for
		for(int row_id = Problem.Variables_num - Problem.Constraints_eq_num; row_id < Problem.Constraint.ie_reduced_matrix.rows(); ++ row_id){
			Problem.Boundary.ie_reduced_matrix.row(row_id) /= Problem.Constraint.norm(row_id);
		}
	}
}

// Covariance between the reduced inequality constraints and the reduced gradient
void LP_constraint_ie_reduced_cov_matrix_generation(LP_object &Problem){
	Problem.Constraint.ie_reduced_cov_matrix = Problem.Constraint.ie_reduced_matrix * Problem.Constraint.ie_reduced_matrix.transpose();
}

// Covariance of the reduced inequality constraints and the reduced gradient
void LP_objective_ie_reduced_cov_matrix_generation(LP_object &Problem){
	Problem.Objective.ie_reduced_cov_vector = Problem.Constraint.ie_reduced_matrix * Problem.Objective.reduced_vector;
}

// Main function for the optimization algorithm
void LP_optimization(LP_object &Problem, bool stepwise_obj){
	// Set precision for 0 detection
	double tol = pow(10, -9);
	double eps = pow(10, -7);
	
	// Declare variables for the main loop
	bool ldlt_flag;
	bool coeff_update_flag;
	int active_constraint_num;
	int active_constraint_max;
	int active_constraint_max_temp;
	double min_increment;
	double current_increment;
	double Previous_Obj = Problem.Solution.reduced_vector.dot(Problem.Objective.reduced_vector);
	std::vector <Eigen::Vector2i> Active_constraint_now;
	Active_constraint_now.reserve(Problem.Variables_num + Problem.Constraints_ie_num);
	std::vector <Eigen::Vector2i> Active_constraint_prior;
	Active_constraint_prior.reserve(Problem.Variables_num + Problem.Constraints_ie_num);
	std::vector <Eigen::Vector2i> Active_constraint_later;
	Active_constraint_later.reserve(Problem.Variables_num + Problem.Constraints_ie_num);
	std::vector <Trip> Subspan_trip;
	Eigen::VectorXd Projected_grad;
	Eigen::VectorXd Projected_increment;
	Eigen::VectorXd redundant_variables;
	Eigen::VectorXi Previous_active_constraint = Eigen::VectorXi::Ones(Problem.Variables_num + Problem.Constraints_ie_num);
	Eigen::MatrixXd Boundary_gap(Problem.Variables_num + Problem.Constraints_ie_num, 2);
	Eigen::SparseMatrix <double> Subspan_matrix;
	Eigen::SparseMatrix <double> Subcov_matrix;
	Eigen::SparseMatrix <double> Subconstraint_matrix;
	
	// Initialize previous active constraint
	Projected_grad = Problem.Objective.reduced_vector;
	if(Projected_grad.norm() != 0){
		Projected_grad /= Projected_grad.norm();
	}
	Projected_increment = Problem.Constraint.ie_reduced_matrix * Projected_grad;
	Boundary_gap.col(0) = Problem.Constraint.ie_reduced_matrix * Problem.Solution.reduced_vector - Problem.Boundary.ie_reduced_matrix.col(0);
	Boundary_gap.col(1) = Problem.Boundary.ie_reduced_matrix.col(1) - Problem.Constraint.ie_reduced_matrix * Problem.Solution.reduced_vector;		
	for(int constraint_iter = 0; constraint_iter < Boundary_gap.rows(); ++ constraint_iter){
		if(Boundary_gap(constraint_iter, 0) < tol && Projected_increment(constraint_iter) < 0){
			Previous_active_constraint(constraint_iter) = 0;
		}
		else if(Boundary_gap(constraint_iter, 1) < tol && Projected_increment(constraint_iter) > 0){
			Previous_active_constraint(constraint_iter) = 0;
		}
	}
	
	int loop_count = 0;
	//while(loop_count < Problem.Variables_num){
	while(1){
		loop_count += 1;
		std::cout << "---------------------------------------------------------------------------" << std::endl;
		std::cout << "New loop" << std::endl;
		std::cout << "---------------------------------------------------------------------------" << std::endl;
		// Clear list of current active constraints
		Active_constraint_now.clear();
		Active_constraint_prior.clear();
		Active_constraint_later.clear();
		
		// Check which constraints are active currently
		std::cout << "\nCheck Active Constraints\n\n";
		Boundary_gap.col(0) = Problem.Constraint.ie_reduced_matrix * Problem.Solution.reduced_vector - Problem.Boundary.ie_reduced_matrix.col(0);
		Boundary_gap.col(1) = Problem.Boundary.ie_reduced_matrix.col(1) - Problem.Constraint.ie_reduced_matrix * Problem.Solution.reduced_vector;
		for(int constraint_iter = 0; constraint_iter < Boundary_gap.rows(); ++ constraint_iter){
			if(Boundary_gap(constraint_iter, 0) < tol){
				Problem.Boundary.ie_reduced_matrix(constraint_iter, 0) += Boundary_gap(constraint_iter, 0);
				Boundary_gap(constraint_iter, 0) = 0.;
				
				if(Previous_active_constraint(constraint_iter) == 0){
					Active_constraint_prior.push_back(Eigen::Vector2i(constraint_iter, 0));
				}
				else{
					Active_constraint_later.push_back(Eigen::Vector2i(constraint_iter, 0));
				}
			}
			else if(Boundary_gap(constraint_iter, 1) < tol){
				Problem.Boundary.ie_reduced_matrix(constraint_iter, 1) -= Boundary_gap(constraint_iter, 1);
				Boundary_gap(constraint_iter, 1) = 0.;
				
				if(Previous_active_constraint(constraint_iter) == 0){
					Active_constraint_prior.push_back(Eigen::Vector2i(constraint_iter, 1));
				}
				else{
					Active_constraint_later.push_back(Eigen::Vector2i(constraint_iter, 1));
				}
			}
		}
		Active_constraint_now.insert(Active_constraint_now.begin(), Active_constraint_prior.begin(), Active_constraint_prior.end());
		Active_constraint_now.insert(Active_constraint_now.end(), Active_constraint_later.begin(), Active_constraint_later.end());

		// Check number of active constraints
		// Active constraints more than dimension of free variables: the active constraints form a degenerate extreme point
		if(Active_constraint_now.size() > Problem.Variables_num - Problem.Constraints_eq_num){
			std::cout << "\nBoundary Relaxed\n\n";
			#pragma omp parallel
			{
				#pragma omp for
				for(int constraint_iter = 0; constraint_iter < Active_constraint_now.size(); ++ constraint_iter){
					if(Active_constraint_now[constraint_iter](1) == 0){
						Problem.Boundary.ie_reduced_matrix(Active_constraint_now[constraint_iter](0), 0) -= eps;
					}
					else{
						Problem.Boundary.ie_reduced_matrix(Active_constraint_now[constraint_iter](0), 1) += eps;
					}
				}
			}
			Previous_active_constraint = Eigen::VectorXi::Zero(Problem.Variables_num + Problem.Constraints_ie_num);
			continue;
		}
		// Active constraints less than dimension of free variables: project on the null space of these constraints
		else if(Active_constraint_now.size() != 0){
			std::cout << "\nAt Boundary\n\n";
			ldlt_flag = 1;
			active_constraint_max = std::min(int (Active_constraint_now.size()), Problem.Variables_num - Problem.Constraints_eq_num - 1);
			Previous_active_constraint = Eigen::VectorXi::Zero(Problem.Variables_num + Problem.Constraints_ie_num);
			Subspan_matrix = Eigen::SparseMatrix <double> (active_constraint_max, Problem.Variables_num + Problem.Constraints_ie_num);
			Subspan_matrix.reserve(Eigen::VectorXi::Constant(Problem.Variables_num + Problem.Constraints_ie_num, 2));

			for(int constraint_iter = 0; constraint_iter < active_constraint_max; ++ constraint_iter){
				Previous_active_constraint(Active_constraint_now[constraint_iter](0)) = 1;
				Subspan_matrix.insert(constraint_iter, Active_constraint_now[constraint_iter](0)) = 1;	
			}
			Subcov_matrix = Subspan_matrix.topRows(active_constraint_max) * Problem.Constraint.ie_reduced_cov_matrix * Subspan_matrix.topRows(active_constraint_max).transpose();
			
			// Check if subspan of covariance matrix is full rank
			std::cout << "\nCheck Full Rank LDLT\n\n";
			Problem.Solver.ldlt.compute(Subcov_matrix);
			if(abs(Problem.Solver.ldlt.determinant()) == 0.){
				// LDLT has numerical stability issues so use qr solver to check for full rank again
				std::cout << "\nCheck Full Rank QR\n\n";
				ldlt_flag = 0;
				Subconstraint_matrix = Subspan_matrix.topRows(active_constraint_max) * Problem.Constraint.ie_reduced_matrix;
				Problem.Solver.qr.compute(Subconstraint_matrix.transpose());
				
				// If the constraints are degenerate, need to pick only part of the constraints for projection
				if(Problem.Solver.qr.rank() < active_constraint_max){
					std::cout << "\nFind subspan of the constraint\n\n";
					Subspan_matrix = Subspan_matrix * Problem.Solver.qr.colsPermutation().transpose();				
					Subcov_matrix = Subspan_matrix.topRows(Problem.Solver.qr.rank()) * Problem.Constraint.ie_reduced_cov_matrix * Subspan_matrix.topRows(Problem.Solver.qr.rank()).transpose();

//					Subconstraint_matrix = Problem.Solver.qr.colsPermutation() * Subconstraint_matrix;
//					Subconstraint_matrix = Subconstraint_matrix.topRows(Problem.Solver.qr.rank());
//					Subcov_matrix = Subconstraint_matrix * Subconstraint_matrix.transpose();
				}
				
				Problem.Solver.qr.compute(Subcov_matrix);					
			}
			
			if(ldlt_flag){
				std::cout << "\nProject with LDLT\n\n";
				Projected_grad = Problem.Objective.reduced_vector;
				Projected_grad -= (Subspan_matrix.topRows(active_constraint_max) * Problem.Constraint.ie_reduced_matrix).transpose() * Problem.Solver.ldlt.solve(Subspan_matrix.topRows(active_constraint_max) * Problem.Objective.ie_reduced_cov_vector);
				if(Projected_grad.norm() != 0){
					Projected_grad /= Projected_grad.norm();
				}
				else{
					// The objective function is degenerate so no further improvement is possible
					break;
				}				
			}
			else{
				std::cout << "\nProject with QR\n\n";
				Projected_grad = Problem.Objective.reduced_vector;
				Projected_grad -= (Subspan_matrix.topRows(active_constraint_max) * Problem.Constraint.ie_reduced_matrix).transpose() * Problem.Solver.qr.solve(Subspan_matrix.topRows(active_constraint_max) * Problem.Objective.ie_reduced_cov_vector);
				if(Projected_grad.norm() != 0){
					Projected_grad /= Projected_grad.norm();
				}
				else{
					// The objective function is degenerate so no further improvement is possible
					break;
				}				
			}			
		}
		// No active constraints; interior point
		else{
			std::cout << "\nAt Interior\n\n";
			Projected_grad = Problem.Objective.reduced_vector;
			if(Projected_grad.norm() != 0){
				Projected_grad /= Projected_grad.norm();
			}
			else{
				// The objective function is degenerate so no further improvement is possible
				break;
			}
		}
		
		// Check the maximum allow increment along the projected gradient greater than 0
		std::cout << "\nCheck Maximum Allowed Increment\n\n";
		min_increment = std::numeric_limits<double>::infinity();
		Projected_increment = Problem.Constraint.ie_reduced_matrix * Projected_grad;
		#pragma omp parallel
		{
			#pragma omp for reduction(min: min_increment) private(current_increment)
			for(int constraint_iter = 0; constraint_iter < Boundary_gap.rows(); ++ constraint_iter){
				if(abs(Projected_increment(constraint_iter)) > tol){
					current_increment = std::max(-Boundary_gap(constraint_iter, 0) / Projected_increment(constraint_iter), Boundary_gap(constraint_iter, 1) / Projected_increment(constraint_iter));
					min_increment = std::min(current_increment, min_increment);												
				}
			}
		}
				
		// Check if there are feasible directions for improvement
		if(min_increment != std::numeric_limits<double>::infinity()){
			Problem.Solution.reduced_vector += min_increment * Projected_grad;
		}
		else{
			break;
		}

		// If the objective is step-wise linear, check if the coefficient should be updated
		std::cout << "\nCheck Coefficient Update\n\n";
		if(stepwise_obj){
			coeff_update_flag = 0;
			Problem.Objective.update_coeff = Eigen::VectorXd::Zero(Problem.Variables_num);
			Boundary_gap.col(0) = Problem.Constraint.ie_reduced_matrix * Problem.Solution.reduced_vector - Problem.Boundary.ie_reduced_matrix.col(0);
			Boundary_gap.col(1) = Problem.Boundary.ie_reduced_matrix.col(1) - Problem.Constraint.ie_reduced_matrix * Problem.Solution.reduced_vector;
			
			for(int constraint_iter = 0; constraint_iter < Problem.Objective.varying_vector.size(); ++ constraint_iter){
				if(Problem.Objective.varying_vector(constraint_iter) == 1. && abs(Projected_increment(constraint_iter)) > tol){
					if(Boundary_gap(constraint_iter, 0) < tol){
						Problem.Objective.update_coeff(constraint_iter) = -1.;
						coeff_update_flag = 1;
					}
					else if(Boundary_gap(constraint_iter, 1) < tol){
						Problem.Objective.update_coeff(constraint_iter) = 1.;
						coeff_update_flag = 1;
					}				
				}
			}
			
			// Exit the while loop if there exist some step-wise obkective coefficients to be updated
			if(coeff_update_flag){
				Previous_Obj = Problem.Solution.reduced_vector.dot(Problem.Objective.reduced_vector);
				Problem.Objective.update_coeff = Problem.Constraint.permutation_matrix * Problem.Objective.update_coeff;
				break;
			}			
		}
		
		// Check if objective value actually improved significantly
		std::cout << "\nCheck Objective Improvement\n\n";
		if(Problem.Solution.reduced_vector.dot(Problem.Objective.reduced_vector) - Previous_Obj > tol){
			// If improved, update the previous solution
			Previous_Obj = Problem.Solution.reduced_vector.dot(Problem.Objective.reduced_vector);
		}
		else{
			Problem.Solution.reduced_vector -= min_increment * Projected_grad;
			break;
		}
	}
	
	// Calculate the objective value of the optimal solution
	Problem.Objective.reduced_value = Previous_Obj;

	// Reconstruct the original solution and objective value from the reduced problem, if necessary
	if(Problem.Constraints_eq_num != 0){
		// Reconstruction of the solution
		// x_r = -A_r^(-1) * (A_e * x_e - b)
		redundant_variables = -Problem.Solver.lu.solve(Problem.Constraint.eq_orig_matrix.leftCols(Problem.Variables_num - Problem.Constraints_eq_num) * Problem.Solution.reduced_vector - Problem.Boundary.eq_vector);
		Problem.Solution.orig_vector = Eigen::VectorXd(Problem.Variables_num);
		Problem.Solution.orig_vector << Problem.Solution.reduced_vector, redundant_variables;

		// Reconstruction of the objective value
		// Z = c_e^T * x_e + c_r^T * x_r 
		//   = c_e^T * x_e - c_r^T * A_r^(-1) * (A_e * x_e - b)
		//   = (c_e^T - c_r^T * A_r^(-1) * A_e) * x_e + c_r^T * A_r^(-1) * b
		Problem.Objective.orig_value = 	Problem.Objective.reduced_value + Problem.Objective.orig_vector.tail(Problem.Constraints_eq_num).transpose() * Problem.Solver.lu.solve(Problem.Boundary.eq_vector);
	}
	else{
		Problem.Solution.orig_vector = Problem.Solution.reduced_vector;
		Problem.Objective.orig_value = Problem.Objective.reduced_value;
	}
	
	// Permute the original solution to the correct order
	Problem.Solution.orig_vector = Problem.Constraint.permutation_matrix * Problem.Solution.orig_vector;	 
}

// Function for the optimization algorithm
void LP_result_print(LP_object &Problem, std::string Problem_name = "Linear Problem"){
	std::cout << std::fixed << std::setprecision(8);
	std::cout << "---------------------------------------------------------------------------" << std::endl;
	std::cout << "| " << Problem_name << std::endl;
	std::cout << "---------------------------------------------------------------------------" << std::endl;
	std::cout << "| Reduced Problem  | " << std::endl;	
	std::cout << "| Solution         | " << Problem.Solution.reduced_vector.transpose() << std::endl;
	std::cout << "| Objective value  | " << Problem.Objective.reduced_value << std::endl;
	std::cout << "---------------------------------------------------------------------------" << std::endl;
	std::cout << "| Original Problem | " << std::endl;
	std::cout << "| Solution         | " << Problem.Solution.orig_vector.transpose() << std::endl;
	std::cout << "| Objective value  | " << Problem.Objective.orig_value << std::endl;
	std::cout << "---------------------------------------------------------------------------" << std::endl;
	std::cout << "\n" << std::endl;
}

// Wrapping up all the process
void LP_process(LP_object &Problem, std::string Problem_name, bool result_output, bool find_sol, bool stepwise_obj, bool constraint_update, bool boundary_update, bool objective_update){
	// Generate sparse matrix for reduced inequality constraints, if needed
	if(constraint_update){
		LP_constraint_eq_redundant_deletion(Problem);	
		LP_variables_permutation(Problem, stepwise_obj);
		LP_constraint_redundant_matrix_solver(Problem);
		LP_constraint_ie_reduced_generation(Problem);
		LP_boundary_ie_reduced_generation(Problem);
		LP_objective_reduced_generation(Problem);
		LP_constraint_ie_reduced_normalization(Problem);
		LP_boundary_ie_reduced_normalization(Problem);
		LP_constraint_ie_reduced_cov_matrix_generation(Problem);
		LP_objective_ie_reduced_cov_matrix_generation(Problem);	
	}
	else{
		// Update reduced boudaries, if needed
		if(boundary_update){
			LP_boundary_ie_reduced_generation(Problem);
			LP_boundary_ie_reduced_normalization(Problem);
		}
		
		// Update reduced inequality increments along the gradient, if needed
		if(objective_update){
			LP_objective_reduced_generation(Problem);
			LP_objective_ie_reduced_cov_matrix_generation(Problem);
		}
	}

	// Solve the LP
	if(find_sol){
		LP_solution_permutation(Problem);
		LP_feasible_solution_reduced_generation(Problem);
		LP_optimization(Problem, stepwise_obj);
	}
	
	// Print results
	if(result_output){
		LP_result_print(Problem, Problem_name);
	}
}