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
void LP_variables_permutation(LP_object &Problem){
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
	Problem.Solution.orig_vector = Problem.Constraint.permutation_matrix.transpose() * Problem.Solution.orig_vector;
}

// Linear system for reducing the redundant variables
void LP_constraint_redundant_matrix_solver(LP_object &Problem){
	Problem.Solver.lu.compute(Problem.Constraint.eq_orig_matrix.rightCols(Problem.Constraints_eq_num));
	//Problem.Solver.lu.compute(Problem.Constraint.eq_orig_matrix.rightCols(Problem.Constraints_eq_num));
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
	//std::cout << std::fixed << std::setprecision(6) << Problem.Constraint.ie_reduced_matrix << "\n" << std::endl;
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

// Covariance of the reduced inequality constraints and the reduced gradient
void LP_constraint_ie_reduced_cov_matrix_generation(LP_object &Problem){
	Problem.Constraint.ie_reduced_cov_matrix = Problem.Constraint.ie_reduced_matrix * Problem.Constraint.ie_reduced_matrix.transpose();
	Problem.Objective.ie_reduced_cov_vector = Problem.Constraint.ie_reduced_matrix * Problem.Objective.reduced_vector;
}

// Main function for the optimization algorithm
void LP_optimization(LP_object &Problem){
	// Set precision for 0 detection
	double tol = pow(10, -12);
	double eps = pow(10, -8);
	
	// Declare variables for the main loop
	int active_constraint_num;
	double min_increment;
	double current_increment;
	double Previous_Obj = -std::numeric_limits<double>::infinity();
	std::vector <Eigen::Vector2i> Active_constraint_now;
	Active_constraint_now.reserve(Problem.Variables_num + Problem.Constraints_ie_num);
	std::vector <Eigen::Vector2i> Active_constraint_next;
	Active_constraint_next.reserve(Problem.Variables_num + Problem.Constraints_ie_num);
	std::vector <Trip> Subspan_trip;
	Eigen::VectorXd Projected_grad;
	Eigen::VectorXd Projected_increment;
	Eigen::MatrixXd Boundary_gap(Problem.Variables_num + Problem.Constraints_ie_num, 2);
	Eigen::SparseMatrix <double> Subspan_matrix;
	Eigen::SparseMatrix <double> Subcov_matrix;
	
	int loop_count = 0;
	//while(loop_count < 5){
	while(1){
		loop_count += 1;
		//std::cout << "---------------------------------------------------------------------------" << std::endl;
		//std::cout << "New loop" << std::endl;
		//std::cout << "---------------------------------------------------------------------------" << std::endl;
		// Clear list of current active constraints
		Active_constraint_now.clear();
		Subspan_matrix = Eigen::SparseMatrix <double> (Problem.Variables_num + Problem.Constraints_ie_num, Problem.Variables_num + Problem.Constraints_ie_num);
		Subspan_matrix.reserve(Eigen::VectorXi::Constant(Problem.Variables_num + Problem.Constraints_ie_num, 2));
		
		// Check which constraints are active currently
		Boundary_gap.col(0) = Problem.Constraint.ie_reduced_matrix * Problem.Solution.reduced_vector - Problem.Boundary.ie_reduced_matrix.col(0);
		Boundary_gap.col(1) = Problem.Boundary.ie_reduced_matrix.col(1) - Problem.Constraint.ie_reduced_matrix * Problem.Solution.reduced_vector;	
		//std::cout << std::fixed << std::setprecision(6) << Boundary_gap << "\n" << std::endl;
		for(int constraint_iter = 0; constraint_iter < Boundary_gap.rows(); ++ constraint_iter){
			if(Boundary_gap(constraint_iter, 0) < tol){
				Boundary_gap(constraint_iter, 0) = tol;
				Active_constraint_now.push_back(Eigen::Vector2i(constraint_iter, 0));
				//std::cout << Active_constraint_now[Active_constraint_now.size() - 1].transpose() << std::endl; 
			}
			else if(Boundary_gap(constraint_iter, 1) < tol){
				Boundary_gap(constraint_iter, 1) = tol;
				Active_constraint_now.push_back(Eigen::Vector2i(constraint_iter, 1));
				//std::cout << Active_constraint_now[Active_constraint_now.size() - 1].transpose() << std::endl;
			}
		}
		//std::cout << std::endl;
		
		// Check if the active constraints form a degenerate extreme point
		if(Active_constraint_now.size() > Problem.Variables_num - Problem.Constraints_eq_num){
			for(int constraint_iter = 0; constraint_iter < Active_constraint_now.size(); ++ constraint_iter){
				if(Active_constraint_now[constraint_iter](1) == 0){
					Problem.Boundary.ie_reduced_matrix(Active_constraint_now[constraint_iter](0), 0) -= eps;
				}
				else{
					Problem.Boundary.ie_reduced_matrix(Active_constraint_now[constraint_iter](0), 1) += eps;
				}
			}
			//std::cout << "\nBoundary Relaxed\n" << std::endl;
			continue;
		}
		
		// Update the active constraints along projected gradient, if the current solution is on the boundary
		if(Active_constraint_now.size() > 0){
			active_constraint_num = 0;
			for(int constraint_iter = 0; constraint_iter < Active_constraint_now.size(); ++ constraint_iter){
				//std::cout << "\n" << constraint_iter << std::endl;;
				Subspan_matrix.insert(active_constraint_num, Active_constraint_now[constraint_iter](0)) = 1;
				Subcov_matrix = Subspan_matrix.topRows(active_constraint_num + 1) * Problem.Constraint.ie_reduced_cov_matrix * Subspan_matrix.topRows(active_constraint_num + 1).transpose();
				//std::cout << Subcov_matrix << "\n" << std::endl;
				
				// Check if subspan of covariance matrix is full rank
				Problem.Solver.ldlt.compute(Subcov_matrix);
				if(abs(Problem.Solver.ldlt.determinant()) > tol){
					// If subspan of covariance matrix is full rank, solve for the projected gradient on the active constraints
					Projected_grad = Problem.Objective.reduced_vector;
					Projected_grad -= (Subspan_matrix.topRows(active_constraint_num + 1) * Problem.Constraint.ie_reduced_matrix).transpose() * Problem.Solver.ldlt.solve(Subspan_matrix.topRows(active_constraint_num + 1) * Problem.Objective.ie_reduced_cov_vector);
					Projected_grad /= Projected_grad.norm();
					
					// Check the minimum allow increment along the projected gradient is greater than 0
					Projected_increment = Problem.Constraint.ie_reduced_matrix * Projected_grad;
					min_increment = std::numeric_limits<double>::infinity();
					#pragma omp parallel
					{
						#pragma omp for reduction(min: min_increment)
						for(int constraint_iter = 0; constraint_iter < Boundary_gap.rows(); ++ constraint_iter){
							if(abs(Projected_increment(constraint_iter)) > tol){
								// Lower bound
								current_increment = -Boundary_gap(constraint_iter, 0) / Projected_increment(constraint_iter);
								if(current_increment > 0 && current_increment < min_increment){
									min_increment = current_increment;
								}
								//std::cout << std::fixed << std::setprecision(6) << current_increment << " ";
								
								// Upper bound
								current_increment = Boundary_gap(constraint_iter, 1) / Projected_increment(constraint_iter);
								if(current_increment > 0 && current_increment < min_increment){
									min_increment = current_increment;
								}
								//std::cout << std::fixed << std::setprecision(6) << current_increment << std::endl;												
							}
						}
					}
					
					// Exit loop if a feasible direction for improvement of solution is found
					if(min_increment > eps){
						//std::cout << std::fixed << std::setprecision(16) << "\n" << min_increment << std::endl;
						break;
					}
					active_constraint_num += 1;
				}
				else{
					// If subspan of covariance matrix is not full rank, remove the current entry for the subspan matrix and move on
					Subspan_matrix.coeffRef(active_constraint_num, Active_constraint_now[constraint_iter](0)) = 0;
				}
			}
			
			// Check if there are feasible directions for improvement
			if(min_increment != std::numeric_limits<double>::infinity()){
				Problem.Solution.reduced_vector += min_increment * Projected_grad;
				//std::cout << std::fixed << std::setprecision(16) << "\n" << min_increment << std::endl;
				//std::cout << std::fixed << std::setprecision(6) << Projected_grad.transpose() << "\n" << std::endl;
			}
			else{
				break;
			}			
		}
		else{
			// If the point is in the interior, use the default gradient as the direct of improvement
			Projected_grad = Problem.Objective.reduced_vector;
			Projected_increment = Problem.Constraint.ie_reduced_matrix * Projected_grad;
			min_increment = std::numeric_limits<double>::infinity();
			#pragma omp parallel
			{
				#pragma omp for reduction(min: min_increment)
				for(int constraint_iter = 0; constraint_iter < Boundary_gap.rows(); ++ constraint_iter){
					if(abs(Projected_increment(constraint_iter)) > tol){
						// Lower bound
						current_increment = -Boundary_gap(constraint_iter, 0) / Projected_increment(constraint_iter);
						if(current_increment > 0 && current_increment < min_increment){
							min_increment = current_increment;
						}
						
						// Upper bound
						current_increment = Boundary_gap(constraint_iter, 1) / Projected_increment(constraint_iter);
						if(current_increment > 0 && current_increment < min_increment){
							min_increment = current_increment;
						}												
					}
				}
			}
			Problem.Solution.reduced_vector += min_increment * Projected_grad;
		}
		
		// Check if objective value actually improved significantly
		if(Problem.Solution.reduced_vector.dot(Problem.Objective.reduced_vector) - Previous_Obj > eps){
			// If improved, update the previous solution
			Previous_Obj = Problem.Solution.reduced_vector.dot(Problem.Objective.reduced_vector);
		}
		else{
			break;
		}
		
		//std::cout << std::fixed << std::setprecision(6) << Problem.Solution.reduced_vector.transpose() << std::endl;
		//std::cout << std::fixed << std::setprecision(6) << Problem.Solution.reduced_vector.dot(Problem.Objective.reduced_vector) << "\n" << std::endl;
	}
	//std::cout << std::fixed << std::setprecision(6) << Problem.Objective.reduced_vector.transpose() << std::endl;
	//std::cout << std::fixed << std::setprecision(6) << Problem.Solution.reduced_vector.transpose() << std::endl;
	//std::cout << std::fixed << std::setprecision(6) << Problem.Solution.reduced_vector.dot(Problem.Objective.reduced_vector) << "\n" << std::endl;
	
	// Calculate the objective value of the optimal solution
	Problem.Objective.reduced_value = Previous_Obj;

	// Reconstruct the original solution and objective value from the reduced problem, if necessary
	if(Problem.Constraints_eq_num != 0){
		// Reconstruction of the solution
		// x_r = -A_r^(-1) * (A_e * x_e - b)
		Eigen::VectorXd redundant_variables = -Problem.Solver.lu.solve(Problem.Constraint.eq_orig_matrix.leftCols(Problem.Variables_num - Problem.Constraints_eq_num) * Problem.Solution.reduced_vector - Problem.Boundary.eq_vector);
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
	std::cout << std::fixed << std::setprecision(3);
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
void LP_process(LP_object &Problem, std::string Problem_name, bool result_output, bool find_sol, bool constraint_update, bool boundary_update, bool objective_update){
	// Generate sparse matrix for reduced inequality constraints, if needed
	if(constraint_update){
		LP_constraint_eq_redundant_deletion(Problem);	
		LP_variables_permutation(Problem);
		LP_constraint_redundant_matrix_solver(Problem);
		LP_constraint_ie_reduced_generation(Problem);
		LP_boundary_ie_reduced_generation(Problem);
		LP_objective_reduced_generation(Problem);
		LP_feasible_solution_reduced_generation(Problem);
		LP_constraint_ie_reduced_normalization(Problem);
		LP_boundary_ie_reduced_normalization(Problem);
		LP_constraint_ie_reduced_cov_matrix_generation(Problem);	
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
		}
	}

	// Solve the LP
	if(find_sol){
		LP_optimization(Problem);
	}
	
	// Print results
	if(result_output){
		LP_result_print(Problem, Problem_name);
	}
}