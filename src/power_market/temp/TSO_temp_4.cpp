void Flow_Based_Market_Optimization(int tick, market_inform &Market, LP_object &Problem){	
	Eigen::MatrixXd bidded_supply = Market.submitted_supply;
	Eigen::MatrixXd bidded_demand = Market.submitted_demand;
	
	// Initial market clearing within each nodes
	Eigen::VectorXi price_ID(Market.num_zone);
	Market_clearing_nodal(tick, Market, price_ID, bidded_supply, bidded_demand);
	//std::cout << price_ID.transpose() << "\n\n";
	
	// Initialization of process variables for the main optimization loop
	Eigen::MatrixXd bidded_supply_export = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_demand_export = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_supply_import = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_demand_import = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_supply_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);
	Eigen::MatrixXd bidded_demand_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);
	Eigen::MatrixXd bidded_total_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);
	Eigen::MatrixXd utility_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);

	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
		// Supply curves
		bidded_supply_export.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter)) = Market.submitted_supply.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter));
		bidded_supply_export(price_ID(zone_iter), zone_iter) = bidded_supply(price_ID(zone_iter), zone_iter);
		bidded_supply_import.col(zone_iter).array().head(price_ID(zone_iter)) = Market.submitted_supply.col(zone_iter).array().head(price_ID(zone_iter));
		bidded_supply_import(price_ID(zone_iter), zone_iter) = Market.submitted_supply(price_ID(zone_iter), zone_iter) - bidded_supply(price_ID(zone_iter), zone_iter);
		bidded_supply_aggregated(0, zone_iter) = -bidded_supply_import.col(zone_iter).sum();
		
		// Demand curves
		bidded_demand_export.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter)) = Market.submitted_demand.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter));
		bidded_demand_export(price_ID(zone_iter), zone_iter) = Market.submitted_demand(price_ID(zone_iter), zone_iter) - bidded_demand(price_ID(zone_iter), zone_iter);
		bidded_demand_import.col(zone_iter).array().head(price_ID(zone_iter)) = Market.submitted_demand.col(zone_iter).array().head(price_ID(zone_iter));
		bidded_demand_import(price_ID(zone_iter), zone_iter) = bidded_demand(price_ID(zone_iter), zone_iter);
		bidded_demand_aggregated(0, zone_iter) = -bidded_demand_import.col(zone_iter).sum();		
	}
	
	// Aggregated import-export curves for supply and demand
	for(int price_iter = 1; price_iter < Market.price_intervals + 3; ++ price_iter){
		bidded_supply_aggregated.row(price_iter) = bidded_supply_aggregated.row(price_iter - 1) + bidded_supply_import.row(price_iter - 1) + bidded_supply_export.row(price_iter - 1);
		bidded_demand_aggregated.row(price_iter) = bidded_demand_aggregated.row(price_iter - 1) + bidded_demand_import.row(price_iter - 1) + bidded_demand_export.row(price_iter - 1);
	}
	bidded_total_aggregated = bidded_supply_aggregated + bidded_demand_aggregated;
	//std::cout << price_ID(0) << "\n\n";
	//std::cout << bidded_total_aggregated.col(0).transpose() << "\n\n";
	//std::cout << bidded_total_aggregated.row(0) << "\n\n";
	//std::cout << bidded_total_aggregated.bottomRows(1) << "\n\n";
	//std::cout << bidded_total_aggregated.col(385).transpose() << "\n\n";
	//std::cout << bidded_supply_aggregated.col(385).transpose() << "\n\n";
	//std::cout << bidded_supply_import.col(385).transpose() << "\n\n";
	//std::cout << bidded_supply_export.col(385).transpose() << "\n\n";
	//std::cout << Market.submitted_supply.col(385).transpose() << "\n\n";
	//std::cout << bidded_total_aggregated.col(390).transpose() << "\n\n";
	//std::cout << bidded_total_aggregated.col(558).transpose() << "\n\n";
	//std::cout << bidded_total_aggregated.col(616).transpose() << "\n\n";
	//std::cout << bidded_total_aggregated.col(630).transpose() << "\n\n";
	
	// Update inequality constraints for import / export
	Problem.Boundary.ie_orig_matrix.col(0).segment(Market.network.num_edges, Market.network.num_vertice) = bidded_total_aggregated.row(0).transpose();
	Problem.Boundary.ie_orig_matrix.col(1).segment(Market.network.num_edges, Market.network.num_vertice) = bidded_total_aggregated.row(Market.price_intervals + 2).transpose();
	//std::cout << Problem.Boundary.ie_orig_matrix << "\n\n";

	// Aggregated utility function for each bidding zone
	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){	
		// Demand
		utility_aggregated(price_ID(zone_iter), zone_iter) = (bidded_demand_import(price_ID(zone_iter), zone_iter) + bidded_supply_import(price_ID(zone_iter), zone_iter)) * Market.bidded_price(price_ID(zone_iter));
		if(price_ID(zone_iter) > 0){
			for(int price_iter = price_ID(zone_iter) - 1; price_iter >= 0; -- price_iter){
				utility_aggregated(price_iter, zone_iter) = utility_aggregated(price_iter + 1, zone_iter) + (bidded_demand_import(price_iter, zone_iter) + bidded_supply_import(price_iter, zone_iter)) * Market.bidded_price(price_iter);
			}
		}
		
		// Supply
		utility_aggregated(price_ID(zone_iter) + 1, zone_iter) = -(bidded_demand_export(price_ID(zone_iter), zone_iter) + bidded_supply_export(price_ID(zone_iter), zone_iter)) * Market.bidded_price(price_ID(zone_iter));
		if(price_ID(zone_iter) < Market.price_intervals){
			for(int price_iter = price_ID(zone_iter) + 2; price_iter < Market.price_intervals + 3; ++ price_iter){
				utility_aggregated(price_iter, zone_iter) = utility_aggregated(price_iter - 1, zone_iter) - (bidded_demand_export(price_iter - 1, zone_iter) + bidded_supply_export(price_iter - 1, zone_iter)) * Market.bidded_price(price_iter - 1);
			}			
		}
	}
	
	// Find nodes with no power source / sink
	Eigen::VectorXi Zero_node = Eigen::VectorXi::Zero(Market.num_zone);
	std::vector <int> Non_zero;
	Non_zero.reserve(Market.num_zone);
	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
		if(abs(bidded_total_aggregated(0, zone_iter) - bidded_total_aggregated(Market.price_intervals + 2, zone_iter)) <= 0.){
			Zero_node(zone_iter) = 1;
		}
		else{
			Non_zero.push_back(zone_iter);
		}
	}
	
	// Initialize metric tensor solver for moving around the degree of freedoms
	std::vector <Trip> Metric_trip;
	Metric_trip.reserve(3 * Non_zero.size() - 5);
	for(int zone_iter = 0; zone_iter < Non_zero.size() - 1; ++ zone_iter){
		Metric_trip.push_back(Trip(zone_iter, zone_iter, 2.));
		if(zone_iter != 0){
			Metric_trip.push_back(Trip(zone_iter, zone_iter - 1, -1.));
		}
		if(zone_iter != Market.num_zone - 2){
			Metric_trip.push_back(Trip(zone_iter, zone_iter + 1, -1.));
		}
	}
	Eigen::SparseMatrix <double> Metric_matrix (Market.num_zone - 1, Market.num_zone - 1);
	Metric_matrix.setFromTriplets(Metric_trip.begin(), Metric_trip.end());
	Market.dof_metric.compute(Metric_matrix);	
	
	// Declare variables for the main loop
	bool break_flag = 0;
	double tol = pow(10., -12.);
	double eps = pow(10., -10.);	
	//double mu = 1. - eps;
	double mu = 0.;
	double dS = 10.;
	double deficit_penalty = 0.;
	double error;
	double obj = 0.;
	double obj_temp;
	Eigen::VectorXi flow_dir(Market.num_zone - 1);
	Eigen::VectorXi price_ID_temp = price_ID;
	// Initial quantity should be randomize to avoid degeneracy
	Eigen::VectorXd quan = Eigen::VectorXd::Zero(Market.num_zone);
//	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
//		quan(zone_iter) = (1 - Zero_node(zone_iter)) * (std::rand() % 1000 - 500) / 10. * eps;
//	}
//	double quan_mean = quan.sum() / (Market.num_zone - Zero_node.sum());
//	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
//		quan(zone_iter) -= (1 - Zero_node(zone_iter)) * quan_mean;
//	}
	Eigen::VectorXd quan_temp = quan;
	Eigen::VectorXd utility = Eigen::VectorXd::Zero(Market.num_zone);
	Eigen::VectorXd utility_temp = Eigen::VectorXd::Zero(Market.num_zone);
	Eigen::VectorXd increment_dot(Non_zero.size() - 1);
	Eigen::VectorXd increment_coeff(Non_zero.size() - 1);
	Eigen::VectorXd grad(Market.num_zone);
	Eigen::VectorXd voltage_temp = Eigen::VectorXd::Zero(Market.network.num_vertice);
	Eigen::VectorXd flow_temp(Market.network.num_edges);
	Eigen::MatrixXd Boundary_gap(Problem.Variables_num + Problem.Constraints_ie_num, 2);
	
	// Initialization of price, utility, and objective
	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
		// Update price after small increase / decrease
		if(price_ID_temp(zone_iter) > -1 && price_ID_temp(zone_iter) < Market.price_intervals + 2){
			if(quan_temp(zone_iter) < bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)){				
				while(quan_temp(zone_iter) < bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)){
					price_ID_temp(zone_iter) -= 1;
					if(price_ID_temp(zone_iter) == -1){
						break;
					}					
				}
			}
			else if(quan_temp(zone_iter) > bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter)){
				while(quan_temp(zone_iter) > bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter)){
					price_ID_temp(zone_iter) += 1;
					if(price_ID_temp(zone_iter) == Market.price_intervals + 2){
						//std::cout << zone_iter << ": " << bidded_total_aggregated(0, zone_iter) << " " << bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter) << " " << quan_temp(zone_iter) << " " << quan_temp(zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter) << "\n";
						break;
					}					
				}
			}				
		}
		else if(quan_temp(zone_iter) < bidded_total_aggregated(0, zone_iter)){
			price_ID_temp(zone_iter) = -1;
		}
		else if(quan_temp(zone_iter) > bidded_total_aggregated(Market.price_intervals + 2, zone_iter)){
			price_ID_temp(zone_iter) = Market.price_intervals + 2;
		}
		
		// Update utility function after small increase / decrease
		//std::cout << "Update utility function after small increase / decrease\n";
		if(price_ID_temp(zone_iter) == -1){
			utility_temp(zone_iter) = utility_aggregated(0, zone_iter);
			utility_temp(zone_iter) -= (quan_temp(zone_iter) - bidded_total_aggregated(0, zone_iter)) * ((deficit_penalty + 1) * Market.price_range_inflex(0) - deficit_penalty * Market.price_range_inflex(1));
			//std::cout << zone_iter << ": " << price_ID_temp(zone_iter) << " " << utility_temp(zone_iter) << "\n";
		}
		else if(price_ID_temp(zone_iter) == Market.price_intervals + 2){
			utility_temp(zone_iter) = utility_aggregated(Market.price_intervals + 2, zone_iter);
			utility_temp(zone_iter) -= (quan_temp(zone_iter) - bidded_total_aggregated(Market.price_intervals + 2, zone_iter)) * ((deficit_penalty + 1) * Market.price_range_inflex(1) - deficit_penalty * Market.price_range_inflex(0));
			//std::cout << zone_iter << ": " << price_ID_temp(zone_iter) << " " << utility_temp(zone_iter) << "\n";
		}
		else{
			if(bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter) != 0){
				utility_temp(zone_iter) = (bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - quan_temp(zone_iter)) * utility_aggregated(price_ID_temp(zone_iter), zone_iter)
					+ (quan_temp(zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter)) * utility_aggregated(price_ID_temp(zone_iter) + 1, zone_iter);				
				utility_temp(zone_iter) /= bidded_total_aggregated(price_ID_temp(zone_iter) + 1, zone_iter) - bidded_total_aggregated(price_ID_temp(zone_iter), zone_iter);
			}
			else{
				utility_temp(zone_iter) = utility_aggregated(price_ID_temp(zone_iter), zone_iter);
			}				
		}
	}

	// Update source / sink, voltage, and power flow
	//std::cout << "Update source / sink, voltage, and power flow\n";
	voltage_temp.tail(Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(quan_temp.tail(Market.network.num_vertice - 1));
	flow_temp = (Problem.Constraint.eq_orig_matrix).topRightCorner(Market.network.num_edges, Market.network.num_vertice) * voltage_temp;		

	// Update errors
	//std::cout << "Update errors\n";
	error = 0.;
	// Power flow errors
	for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
		error += pow(std::min(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 0), 0.) + std::max(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 1), 0.), 2.);
	}
	// Voltage errors
	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		error += pow(std::min(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 0), 0.) + std::max(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 1), 0.), 2.);
	}
	
	// Update objective, price, and utility
	obj_temp = (1. - mu) * utility_temp.sum() - mu * error;
	price_ID = price_ID_temp;
	utility = utility_temp;
	obj = obj_temp;
	//std::cout << utility_temp.transpose() << "\n\n";
	std::cout << obj << "\n\n";			

	int loop_count = 0;
	//while(loop_count < 10000){
	while(!break_flag){
		std::cout << "---------------------------------------------------------------------------\n";
		std::cout << loop_count << "\n";
		std::cout << "---------------------------------------------------------------------------\n";
		loop_count += 1;
		break_flag = 1;
		//std::cout << "Search all directions\n";
		for(int zone_iter = 0; zone_iter < Non_zero.size() - 1; ++ zone_iter){
			// Determine flow direction
			//std::cout << "Determine flow direction\n";
			if(bidded_total_aggregated(0, Non_zero[zone_iter]) + bidded_total_aggregated(Market.price_intervals + 2, Non_zero[zone_iter]) != 0. && bidded_total_aggregated(0, Non_zero[zone_iter] + 1) + bidded_total_aggregated(Market.price_intervals + 2, Non_zero[zone_iter] + 1) != 0.){
				if(price_ID(Non_zero[zone_iter]) >= price_ID(Non_zero[zone_iter] + 1)){
					flow_dir(Non_zero[zone_iter]) = -1;
					//std::cout << "Negative flow\n";
					//std::cout << price_ID(Non_zero[zone_iter]) << " " << price_ID(Non_zero[zone_iter] + 1) << "\n";
				}
				else if(price_ID(Non_zero[zone_iter]) < price_ID(Non_zero[zone_iter] + 1)){
					flow_dir(Non_zero[zone_iter]) = 1;
					//std::cout << "Positive flow\n";
					//std::cout << price_ID(Non_zero[zone_iter]) << " " << price_ID(Non_zero[zone_iter] + 1) << "\n";
				}
//				else{
//					increment_dot(zone_iter) = 0.;
//					//std::cout << "No flow\n\n";
//					//flow_dir(Non_zero[zone_iter]) = 0;
//					continue;
//				}				
			}
			else{
				increment_dot(zone_iter) = 0.;
				//std::cout << "No flow\n\n";
				//flow_dir(Non_zero[zone_iter]) = 0;
				continue;				
			}
			//flow_dir(Non_zero[zone_iter]) = 1;
					
			// Update quantity after small increase / decrease
			//std::cout << "Update quantity after small increase / decrease\n";
			quan_temp(Non_zero[zone_iter]) = quan(Non_zero[zone_iter]) + flow_dir(Non_zero[zone_iter]) * dS;
			quan_temp(Non_zero[zone_iter] + 1) = quan(Non_zero[zone_iter] + 1) - flow_dir(Non_zero[zone_iter]) * dS;
			//std::cout << quan_temp.transpose() << "\n";
			
			// Update price after small increase / decrease
			//std::cout << "Update price after small increase / decrease\n";
			if(price_ID_temp(Non_zero[zone_iter]) > -1 && price_ID_temp(Non_zero[zone_iter]) < Market.price_intervals + 2){
				if(quan_temp(Non_zero[zone_iter]) < bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])){				
					while(quan_temp(Non_zero[zone_iter]) < bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])){
						price_ID_temp(Non_zero[zone_iter]) -= 1;
						if(price_ID_temp(Non_zero[zone_iter]) == -1){
							break;
						}					
					}
				}
				else if(quan_temp(Non_zero[zone_iter]) > bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter])){
					while(quan_temp(Non_zero[zone_iter]) > bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter])){
						price_ID_temp(Non_zero[zone_iter]) += 1;
						if(price_ID_temp(Non_zero[zone_iter]) == Market.price_intervals + 2){
							break;
						}					
					}
				}				
			}
			else if(quan_temp(Non_zero[zone_iter]) < bidded_total_aggregated(0, Non_zero[zone_iter])){
				price_ID_temp(Non_zero[zone_iter]) = -1;
			}
			else if(quan_temp(Non_zero[zone_iter]) > bidded_total_aggregated(Market.price_intervals + 2, Non_zero[zone_iter])){
				price_ID_temp(Non_zero[zone_iter]) = Market.price_intervals + 2;
			}
			if(price_ID_temp(Non_zero[zone_iter] + 1) > -1 && price_ID_temp(Non_zero[zone_iter] + 1) < Market.price_intervals + 2){
				if(quan_temp(Non_zero[zone_iter] + 1) < bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1), Non_zero[zone_iter] + 1)){
					while(quan_temp(Non_zero[zone_iter] + 1) < bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1), Non_zero[zone_iter] + 1)){
						price_ID_temp(Non_zero[zone_iter] + 1) -= 1;
						if(price_ID_temp(Non_zero[zone_iter] + 1) == -1){
							break;
						}					
					}
				}
				else if(quan_temp(Non_zero[zone_iter] + 1) > bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1) + 1, Non_zero[zone_iter] + 1)){
					while(quan_temp(Non_zero[zone_iter] + 1) > bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1) + 1, Non_zero[zone_iter] + 1)){
						price_ID_temp(Non_zero[zone_iter] + 1) += 1;
						if(price_ID_temp(Non_zero[zone_iter] + 1) == Market.price_intervals + 2){
							break;
						}					
					}
				}		
			}
			else if(quan_temp(Non_zero[zone_iter] + 1) < bidded_total_aggregated(0, Non_zero[zone_iter] + 1)){
				price_ID_temp(Non_zero[zone_iter] + 1) = -1;
			}
			else if(quan_temp(Non_zero[zone_iter] + 1) > bidded_total_aggregated(Market.price_intervals + 2, Non_zero[zone_iter] + 1)){
				price_ID_temp(Non_zero[zone_iter] + 1) = Market.price_intervals + 2;
			}			
			
			// Update utility function after small increase / decrease
			//std::cout << "Update utility function after small increase / decrease\n";
			if(price_ID_temp(Non_zero[zone_iter]) == -1){
				utility_temp(Non_zero[zone_iter]) = utility_aggregated(0, Non_zero[zone_iter]);
				utility_temp(Non_zero[zone_iter]) -= (quan_temp(Non_zero[zone_iter]) - bidded_total_aggregated(0, Non_zero[zone_iter])) * ((deficit_penalty + 1) * Market.price_range_inflex(0) - deficit_penalty * Market.price_range_inflex(1));
			}
			else if(price_ID_temp(Non_zero[zone_iter]) == Market.price_intervals + 2){
				utility_temp(Non_zero[zone_iter]) = utility_aggregated(Market.price_intervals + 2, Non_zero[zone_iter]);
				utility_temp(Non_zero[zone_iter]) -= (quan_temp(Non_zero[zone_iter]) - bidded_total_aggregated(Market.price_intervals + 2, Non_zero[zone_iter])) * ((deficit_penalty + 1) * Market.price_range_inflex(1) - deficit_penalty * Market.price_range_inflex(0));
			}
			else{
				utility_temp(Non_zero[zone_iter]) = (bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]) - quan_temp(Non_zero[zone_iter])) * utility_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])
					+ (quan_temp(Non_zero[zone_iter]) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])) * utility_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]);
				if(bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]) != 0){
					utility_temp(Non_zero[zone_iter]) /= bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]);
				}				
			}
			if(price_ID_temp(Non_zero[zone_iter] + 1) == -1){
				utility_temp(Non_zero[zone_iter] + 1) = utility_aggregated(0, Non_zero[zone_iter] + 1);
				utility_temp(Non_zero[zone_iter] + 1) -= (quan_temp(Non_zero[zone_iter] + 1) - bidded_total_aggregated(0, Non_zero[zone_iter] + 1)) * ((deficit_penalty + 1) * Market.price_range_inflex(0) - deficit_penalty * Market.price_range_inflex(1));				
			}
			else if(price_ID_temp(Non_zero[zone_iter] + 1) == Market.price_intervals + 2){
				utility_temp(Non_zero[zone_iter] + 1) = utility_aggregated(Market.price_intervals + 2, Non_zero[zone_iter] + 1);
				utility_temp(Non_zero[zone_iter] + 1) -= (quan_temp(Non_zero[zone_iter] + 1) - bidded_total_aggregated(Market.price_intervals + 2, Non_zero[zone_iter] + 1)) * ((deficit_penalty + 1) * Market.price_range_inflex(1) - deficit_penalty * Market.price_range_inflex(0));				
			}
			else{
				utility_temp(Non_zero[zone_iter] + 1) = (bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1) + 1, Non_zero[zone_iter] + 1) - quan_temp(Non_zero[zone_iter] + 1)) * utility_aggregated(price_ID_temp(Non_zero[zone_iter] + 1), Non_zero[zone_iter] + 1)
					+ (quan_temp(Non_zero[zone_iter] + 1) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1), Non_zero[zone_iter] + 1)) * utility_aggregated(price_ID_temp(Non_zero[zone_iter] + 1) + 1, Non_zero[zone_iter] + 1);
				if(bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1) + 1, Non_zero[zone_iter] + 1) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1), Non_zero[zone_iter] + 1) != 0){
					utility_temp(Non_zero[zone_iter] + 1) /= bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1) + 1, Non_zero[zone_iter] + 1) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1), Non_zero[zone_iter] + 1);	
				}
			}
			
			// Update source / sink, voltage, and power flow
			//std::cout << "Update source / sink, voltage, and power flow\n";
			voltage_temp.tail(Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(quan_temp.tail(Market.network.num_vertice - 1));
			flow_temp = (Problem.Constraint.eq_orig_matrix).topRightCorner(Market.network.num_edges, Market.network.num_vertice) * voltage_temp;		

			// Update errors
			//std::cout << "Update errors\n";
			error = 0.;
			// Power flow errors
			for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
				error += pow(std::min(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 0), 0.) + std::max(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 1), 0.), 2.);
			}
			// Voltage errors
			for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
				error += pow(std::min(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 0), 0.) + std::max(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 1), 0.), 2.);
			}		
			
			// Update objective increment for this degree of freedom
			//std::cout << "Update objective increment for this degree of freedom\n\n";
			obj_temp = (1. - mu) * utility_temp.sum() - mu * error;
			increment_dot(zone_iter) = flow_dir(Non_zero[zone_iter]) * (obj_temp - obj) * (obj_temp >= obj);
			//std::cout << obj_temp - obj << " ";
			//increment_dot(zone_iter) = flow_dir(Non_zero[zone_iter]) * (obj_temp - obj);// * (obj_temp >= obj);
			//std::cout << (1. - mu) * utility_temp.sum() << "\n";
			//std::cout << mu * error << "\n\n";
			
			// Return to original values
			quan_temp(Non_zero[zone_iter]) = quan(Non_zero[zone_iter]);
			quan_temp(Non_zero[zone_iter] + 1) = quan(Non_zero[zone_iter] + 1);			
			price_ID_temp(Non_zero[zone_iter]) = price_ID(Non_zero[zone_iter]);
			price_ID_temp(Non_zero[zone_iter] + 1) = price_ID(Non_zero[zone_iter] + 1);
			utility_temp(Non_zero[zone_iter]) = utility(Non_zero[zone_iter]);
			utility_temp(Non_zero[zone_iter] + 1) = utility(Non_zero[zone_iter] + 1);																
		}
		//std::cout << "\n\n";
		
		// Update gradient direction
		//std::cout << "Update gradient direction\n";
		increment_coeff = Market.dof_metric.solve(increment_dot);
		grad = Eigen::VectorXd::Zero(Market.num_zone);
		for(int zone_iter = 0; zone_iter < Non_zero.size() - 1; ++ zone_iter){
			grad(Non_zero[zone_iter]) += increment_coeff(Non_zero[zone_iter]);
			grad(Non_zero[zone_iter] + 1) -= increment_coeff(Non_zero[zone_iter]);
		}
		std::cout << increment_dot.transpose() << "\n\n";
		std::cout << increment_coeff.transpose() << "\n\n";
		//std::cout << grad.norm() << "\n";
		//grad -= Zero_constraint.transpose() * Problem.Solver.qr.solve(Zero_constraint * grad);
		//std::cout << grad.norm() << "\n\n";
		if(grad.norm() > tol){
			grad /= grad.norm();
		}
		//std::cout << increment_coeff.transpose() << "\n";
		//std::cout << grad.transpose() << "\n\n";
		
		// Update quantity on gradient direction
		//std::cout << "Update along gradient direction\n";
		quan_temp = quan + grad * dS;
		
		// Update price after small increase / decrease
		for(int zone_iter = 0; zone_iter < Non_zero.size() - 1; ++ zone_iter){
			//std::cout << "Update price after small increase / decrease\n";
			if(price_ID_temp(Non_zero[zone_iter]) > -1 && price_ID_temp(Non_zero[zone_iter]) < Market.price_intervals + 2){
				if(quan_temp(Non_zero[zone_iter]) < bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])){
					while(quan_temp(Non_zero[zone_iter]) < bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])){
						price_ID_temp(Non_zero[zone_iter]) -= 1;
						if(price_ID_temp(Non_zero[zone_iter]) == -1){
							break;
						}					
					}
				}
				else if(quan_temp(Non_zero[zone_iter]) > bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter])){
					while(quan_temp(Non_zero[zone_iter]) > bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter])){
						price_ID_temp(Non_zero[zone_iter]) += 1;
						if(price_ID_temp(Non_zero[zone_iter]) == Market.price_intervals + 2){
							break;
						}					
					}
				}				
			}
			else if(quan_temp(Non_zero[zone_iter]) < bidded_total_aggregated(0, Non_zero[zone_iter])){
				price_ID_temp(Non_zero[zone_iter]) = -1;
			}
			else if(quan_temp(Non_zero[zone_iter]) > bidded_total_aggregated(Market.price_intervals + 2, Non_zero[zone_iter])){
				price_ID_temp(Non_zero[zone_iter]) = Market.price_intervals + 2;
			}				
			if(price_ID_temp(Non_zero[zone_iter] + 1) > -1 && price_ID_temp(Non_zero[zone_iter] + 1) < Market.price_intervals + 2){
				if(quan_temp(Non_zero[zone_iter] + 1) < bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1), Non_zero[zone_iter] + 1)){
					while(quan_temp(Non_zero[zone_iter] + 1) < bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1), Non_zero[zone_iter] + 1)){
						price_ID_temp(Non_zero[zone_iter] + 1) -= 1;
						if(price_ID_temp(Non_zero[zone_iter] + 1) == -1){
							break;
						}					
					}
				}
				else if(quan_temp(Non_zero[zone_iter] + 1) > bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1) + 1, Non_zero[zone_iter] + 1)){
					while(quan_temp(Non_zero[zone_iter] + 1) > bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1) + 1, Non_zero[zone_iter] + 1)){
						price_ID_temp(Non_zero[zone_iter] + 1) += 1;
						if(price_ID_temp(Non_zero[zone_iter] + 1) == Market.price_intervals + 2){
							break;
						}					
					}
				}				
			}
			else if(quan_temp(Non_zero[zone_iter] + 1) < bidded_total_aggregated(0, Non_zero[zone_iter] + 1)){
				price_ID_temp(Non_zero[zone_iter] + 1) = -1;
			}
			else if(quan_temp(Non_zero[zone_iter] + 1) > bidded_total_aggregated(Market.price_intervals + 2, Non_zero[zone_iter] + 1)){
				price_ID_temp(Non_zero[zone_iter] + 1) = Market.price_intervals + 2;
			}
			
			// Update utility function after small increase / decrease
			//std::cout << "Update utility function after small increase / decrease\n";
			if(price_ID_temp(Non_zero[zone_iter]) == -1){
				utility_temp(Non_zero[zone_iter]) = utility_aggregated(0, Non_zero[zone_iter]);
				utility_temp(Non_zero[zone_iter]) -= (quan_temp(Non_zero[zone_iter]) - bidded_total_aggregated(0, Non_zero[zone_iter])) * ((deficit_penalty + 1) * Market.price_range_inflex(0) - deficit_penalty * Market.price_range_inflex(1));
			}
			else if(price_ID_temp(Non_zero[zone_iter]) == Market.price_intervals + 2){
				utility_temp(Non_zero[zone_iter]) = utility_aggregated(Market.price_intervals + 2, Non_zero[zone_iter]);
				utility_temp(Non_zero[zone_iter]) -= (quan_temp(Non_zero[zone_iter]) - bidded_total_aggregated(Market.price_intervals + 2, Non_zero[zone_iter])) * ((deficit_penalty + 1) * Market.price_range_inflex(1) - deficit_penalty * Market.price_range_inflex(0));
			}
			else{
				utility_temp(Non_zero[zone_iter]) = (bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]) - quan_temp(Non_zero[zone_iter])) * utility_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])
					+ (quan_temp(Non_zero[zone_iter]) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])) * utility_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]);
				if(bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]) != 0){
					utility_temp(Non_zero[zone_iter]) /= bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]);					
				}				
			}
			if(price_ID_temp(Non_zero[zone_iter] + 1) == -1){
				utility_temp(Non_zero[zone_iter] + 1) = utility_aggregated(0, Non_zero[zone_iter] + 1);
				utility_temp(Non_zero[zone_iter] + 1) -= (quan_temp(Non_zero[zone_iter] + 1) - bidded_total_aggregated(0, Non_zero[zone_iter] + 1)) * ((deficit_penalty + 1) * Market.price_range_inflex(0) - deficit_penalty * Market.price_range_inflex(1));				
			}
			else if(price_ID_temp(Non_zero[zone_iter] + 1) == Market.price_intervals + 2){
				utility_temp(Non_zero[zone_iter] + 1) = utility_aggregated(Market.price_intervals + 2, Non_zero[zone_iter] + 1);
				utility_temp(Non_zero[zone_iter] + 1) -= (quan_temp(Non_zero[zone_iter] + 1) - bidded_total_aggregated(Market.price_intervals + 2, Non_zero[zone_iter] + 1)) * ((deficit_penalty + 1) * Market.price_range_inflex(1) - deficit_penalty * Market.price_range_inflex(0));				
			}
			else{
				utility_temp(Non_zero[zone_iter] + 1) = (bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1) + 1, Non_zero[zone_iter] + 1) - quan_temp(Non_zero[zone_iter] + 1)) * utility_aggregated(price_ID_temp(Non_zero[zone_iter] + 1), Non_zero[zone_iter] + 1)
					+ (quan_temp(Non_zero[zone_iter] + 1) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1), Non_zero[zone_iter] + 1)) * utility_aggregated(price_ID_temp(Non_zero[zone_iter] + 1) + 1, Non_zero[zone_iter] + 1);
				if(bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1) + 1, Non_zero[zone_iter] + 1) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1), Non_zero[zone_iter] + 1) != 0.){
					utility_temp(Non_zero[zone_iter] + 1) /= bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1) + 1, Non_zero[zone_iter] + 1) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter] + 1), Non_zero[zone_iter] + 1);
				}
			}		
		}
		
		// Update source / sink, voltage, and power flow
		//std::cout << "Update source / sink, voltage, and power flow\n";
		voltage_temp.tail(Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(quan_temp.tail(Market.network.num_vertice - 1));
		flow_temp = (Problem.Constraint.eq_orig_matrix).topRightCorner(Market.network.num_edges, Market.network.num_vertice) * voltage_temp;		

		// Update errors
		//std::cout << "Update errors\n";
		error = 0.;
		// Power flow errors
		for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
			error += pow(std::min(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 0), 0.) + std::max(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 1), 0.), 2.);
		}
		// Voltage errors
		for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
			error += pow(std::min(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 0), 0.) + std::max(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 1), 0.), 2.);
		}
		
		// Update objective
		//std::cout << "Update objective\n";
		obj_temp = (1. - mu) * utility_temp.sum() - mu * error;
		
		// Check if objective function is actually improved
		//std::cout << "Check if objective function is actually improved\n";
		if(obj_temp > obj){			
			quan = quan_temp;
			price_ID = price_ID_temp;
			utility = utility_temp;
			obj = obj_temp;
			std::cout << "Objective function updated\n";
			std::cout << obj << "\n\n";
			//std::cout << price_ID.transpose() << "\n";
			//std::cout << quan.transpose() << "\n";
			//std::cout << obj_temp - obj << "\n\n";			
		}
		else{
			dS /= 2.;
			std::cout << "Objective function not updated\n";
			//std::cout << quan_temp.transpose() << "\n";
			//std::cout << quan.transpose() << "\n";
			//std::cout << utility_temp.transpose() << "\n";
			//std::cout << utility_temp.minCoeff() << "\n";
			std::cout << obj_temp << "\n\n";
			//std::cout << dS << "\n\n";
			
			// Return to original values
			quan_temp = quan;		
			price_ID_temp = price_ID;
			utility_temp = utility;
		}
		
		if(dS > .0001){
			break_flag = 0;
		}
	}
	
	// Update final source / sink, voltage, and power flow
	//std::cout << "Update final source / sink, voltage, and power flow\n";
	Problem.Solution.orig_vector.segment(Market.network.num_edges, Market.network.num_vertice) = quan;
	Problem.Solution.orig_vector.tail(Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(Problem.Solution.orig_vector.segment(Market.network.num_edges + 1, Market.network.num_vertice - 1));
	Problem.Solution.orig_vector.head(Market.network.num_edges) = (Problem.Constraint.eq_orig_matrix).topRightCorner(Market.network.num_edges, Market.network.num_vertice) * Problem.Solution.orig_vector.tail(Market.network.num_vertice);		
	//std::cout << price_ID.transpose() << "\n\n";
	std::cout << price_ID.maxCoeff() << "\n";
	std::cout << price_ID.minCoeff() << "\n\n";
	//std::cout << Problem.Solution.orig_vector.segment(Market.network.num_edges, Market.network.num_vertice).transpose() << "\n\n";
	//std::cout << Problem.Solution.orig_vector.transpose() << "\n\n";
	std::cout << quan.minCoeff() << " " << quan.maxCoeff() << "\n";
	std::cout << (Problem.Solution.orig_vector.tail(Market.network.num_vertice - 1)).minCoeff() << " " << (Problem.Solution.orig_vector.tail(Market.network.num_vertice - 1)).maxCoeff() << "\n";
	std::cout << (Problem.Solution.orig_vector.head(Market.network.num_edges)).minCoeff() << " " << (Problem.Solution.orig_vector.head(Market.network.num_edges)).maxCoeff() << "\n\n";		
}

void Flow_Based_Market_Optimization_Test_4(int tick, market_inform &Market, LP_object &Problem){
	Eigen::MatrixXd bidded_supply = Market.submitted_supply;
	Eigen::MatrixXd bidded_demand = Market.submitted_demand;
	
	// Initial market clearing within each nodes
	Eigen::VectorXi price_ID(Market.num_zone);
	Market_clearing_nodal(tick, Market, price_ID, bidded_supply, bidded_demand);
	//std::cout << price_ID.transpose() << "\n\n";
	
	// Initialization of process variables for the main optimization loop
	Eigen::MatrixXd bidded_supply_export = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_demand_export = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_supply_import = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_demand_import = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_supply_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);
	Eigen::MatrixXd bidded_demand_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);
	Eigen::MatrixXd bidded_total_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);
	Eigen::MatrixXd utility_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);

	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
		// Supply curves
		bidded_supply_export.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter)) = Market.submitted_supply.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter));
		bidded_supply_export(price_ID(zone_iter), zone_iter) = bidded_supply(price_ID(zone_iter), zone_iter);
		bidded_supply_import.col(zone_iter).array().head(price_ID(zone_iter)) = Market.submitted_supply.col(zone_iter).array().head(price_ID(zone_iter));
		bidded_supply_import(price_ID(zone_iter), zone_iter) = Market.submitted_supply(price_ID(zone_iter), zone_iter) - bidded_supply(price_ID(zone_iter), zone_iter);
		bidded_supply_aggregated(0, zone_iter) = -bidded_supply_import.col(zone_iter).sum();
		
		// Demand curves
		bidded_demand_export.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter)) = Market.submitted_demand.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter));
		bidded_demand_export(price_ID(zone_iter), zone_iter) = Market.submitted_demand(price_ID(zone_iter), zone_iter) - bidded_demand(price_ID(zone_iter), zone_iter);
		bidded_demand_import.col(zone_iter).array().head(price_ID(zone_iter)) = Market.submitted_demand.col(zone_iter).array().head(price_ID(zone_iter));
		bidded_demand_import(price_ID(zone_iter), zone_iter) = bidded_demand(price_ID(zone_iter), zone_iter);
		bidded_demand_aggregated(0, zone_iter) = -bidded_demand_import.col(zone_iter).sum();		
	}
	
	// Aggregated import-export curves for supply and demand
	for(int price_iter = 1; price_iter < Market.price_intervals + 3; ++ price_iter){
		bidded_supply_aggregated.row(price_iter) = bidded_supply_aggregated.row(price_iter - 1) + bidded_supply_import.row(price_iter - 1) + bidded_supply_export.row(price_iter - 1);
		bidded_demand_aggregated.row(price_iter) = bidded_demand_aggregated.row(price_iter - 1) + bidded_demand_import.row(price_iter - 1) + bidded_demand_export.row(price_iter - 1);
	}
	bidded_total_aggregated = bidded_supply_aggregated + bidded_demand_aggregated;
	
	// Update inequality constraints for import / export
	Problem.Boundary.ie_orig_matrix.col(0).segment(Market.network.num_edges, Market.network.num_vertice) = bidded_total_aggregated.row(0).transpose();
	Problem.Boundary.ie_orig_matrix.col(1).segment(Market.network.num_edges, Market.network.num_vertice) = bidded_total_aggregated.row(Market.price_intervals + 2).transpose();

	// Aggregated utility function for each bidding zone
	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){	
		// Demand
		utility_aggregated(price_ID(zone_iter), zone_iter) = (bidded_demand_import(price_ID(zone_iter), zone_iter) + bidded_supply_import(price_ID(zone_iter), zone_iter)) * Market.bidded_price(price_ID(zone_iter));
		if(price_ID(zone_iter) > 0){
			for(int price_iter = price_ID(zone_iter) - 1; price_iter >= 0; -- price_iter){
				utility_aggregated(price_iter, zone_iter) = utility_aggregated(price_iter + 1, zone_iter) + (bidded_demand_import(price_iter, zone_iter) + bidded_supply_import(price_iter, zone_iter)) * Market.bidded_price(price_iter);
			}
		}
		
		// Supply
		utility_aggregated(price_ID(zone_iter) + 1, zone_iter) = -(bidded_demand_export(price_ID(zone_iter), zone_iter) + bidded_supply_export(price_ID(zone_iter), zone_iter)) * Market.bidded_price(price_ID(zone_iter));
		if(price_ID(zone_iter) < Market.price_intervals){
			for(int price_iter = price_ID(zone_iter) + 2; price_iter < Market.price_intervals + 3; ++ price_iter){
				utility_aggregated(price_iter, zone_iter) = utility_aggregated(price_iter - 1, zone_iter) - (bidded_demand_export(price_iter - 1, zone_iter) + bidded_supply_export(price_iter - 1, zone_iter)) * Market.bidded_price(price_iter - 1);
			}			
		}
	}

	// Find nodes with no power source / sink
	std::vector <int> Non_zero;
	Non_zero.reserve(Market.num_zone);
	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
		if(abs(bidded_total_aggregated(0, zone_iter) - bidded_total_aggregated(Market.price_intervals + 2, zone_iter)) > 0.){
			Non_zero.push_back(zone_iter);
		}
	}
	
	// Declare variables for the initial solver
	double ratio;
	double source_error_weight = 1.;
	double error;
	double eps = pow(10., -12.);
	double mu = 1. - eps;
//	double dmu = .00001;
//	double mu = .99;
	double obj; 	
	Eigen::VectorXd quan(Market.num_zone);
	Eigen::VectorXd utility(Market.num_zone);
	Eigen::VectorXd voltage = Eigen::VectorXd::Zero(Market.num_zone);
	Eigen::VectorXd flow = Eigen::VectorXd::Zero(Market.network.num_edges);
	
	// Find optimal solution without voltage and power flow constraints
	Eigen::VectorXd bidded_total_aggregated_all = bidded_total_aggregated * Eigen::VectorXd::Constant(Market.num_zone, 1.);
	for(int price_iter = 0; price_iter < Market.price_intervals + 2; ++ price_iter){
		if(bidded_total_aggregated_all(price_iter) <= 0 && bidded_total_aggregated_all(price_iter + 1) >= 0){
			price_ID = Eigen::VectorXi::Constant(Market.num_zone, price_iter);
			if(bidded_total_aggregated_all(price_iter + 1) > bidded_total_aggregated_all(price_iter)){
				ratio = bidded_total_aggregated_all(price_iter + 1) / (bidded_total_aggregated_all(price_iter + 1) - bidded_total_aggregated_all(price_iter));				
			}
			else{
				ratio = .5;
			}
			quan = ratio * bidded_total_aggregated.row(price_iter) + (1. - ratio) * bidded_total_aggregated.row(price_iter + 1);
			
			break;
		}
	}
	
	// Find initial values for utility, voltage, and power flow
	utility = ratio * utility_aggregated.row(price_ID(0)) + (1. - ratio) * utility_aggregated.row(price_ID(0) + 1);
	voltage.tail(Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(quan.tail(Market.network.num_vertice - 1));
	flow = (Problem.Constraint.eq_orig_matrix).topRightCorner(Market.network.num_edges, Market.network.num_vertice) * voltage;
	
	// Find initial values for error and objective
	error = 0.;
	// Power flow errors
	for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
		error += pow(std::min(flow(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 0), 0.) + std::max(flow(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 1), 0.), 2.);
	}
	// Power errors
	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		error += source_error_weight * pow(std::min(quan(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + node_iter, 0), 0.) + std::max(quan(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + node_iter, 1), 0.), 2.);
	}
	// Voltage errors
	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		error += pow(std::min(voltage(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 0), 0.) + std::max(voltage(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 1), 0.), 2.);
	}
	// Exit function if no constraint is violated
	if(error == 0.){
		// Update confirmed bids and prices
		return;
	}
	obj = (1. - mu) * utility.sum() - mu * error;
	
	// Declare variables for the loop
	int updated_num;
	int node_ref_ID;
	int loop_count;
	double dS_max = .1 * (bidded_total_aggregated.bottomRows(1)).lpNorm <1> () / bidded_total_aggregated.cols();
	double dS = dS_max;
	double dQ;
	//double alpha;
	double obj_temp;
	double error_temp;
	Eigen::VectorXi price_ID_temp = price_ID;
	Eigen::VectorXd quan_temp = quan;
	Eigen::VectorXd voltage_temp = voltage;
	Eigen::VectorXd flow_temp = flow;
	Eigen::VectorXd utility_temp = utility;
	
	std::cout << "\n" << loop_count << ": " << error << "; " << dS << "\n";	
	loop_count = 0;
	//while(dS > eps * dS_max && loop_count < 1000){
	while(loop_count < 1000){
		loop_count += 1;
		updated_num = 0;
		for(int zone_iter = 0; zone_iter < Non_zero.size(); ++ zone_iter){
			node_ref_ID = std::rand() % Non_zero.size();
			while(zone_iter == node_ref_ID){
				node_ref_ID = std::rand() % Non_zero.size();
			}
			
			for(int dir_iter = 0; dir_iter < 2; ++ dir_iter){
				// Find possible increment within a price range
				if(dir_iter == 0){
					if(quan_temp(Non_zero[zone_iter]) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]) >= bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID]) - quan_temp(Non_zero[node_ref_ID])){
						dQ = quan_temp(Non_zero[zone_iter]) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]);
					}
					else{
						dQ = bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID]) - quan_temp(Non_zero[node_ref_ID]);
					}
				}
				else{
					if(bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]) - quan_temp(Non_zero[zone_iter]) >= quan_temp(Non_zero[node_ref_ID]) - bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID])){
						dQ = bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]) - quan_temp(Non_zero[zone_iter]);
					}
					else{
						dQ = quan_temp(Non_zero[node_ref_ID]) - bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID]);
					}
				}
				dQ *= .5;
				dQ += eps;
				if(dQ < 0.){
					dQ = -eps;
				}
							
				// Update quantity after small increase / decrease
				quan_temp(Non_zero[zone_iter]) += 2 * (dir_iter - .5) * dS; 
				quan_temp(Non_zero[node_ref_ID]) -= 2 * (dir_iter - .5) * dS;
				//quan_temp(Non_zero[zone_iter]) += 2 * (dir_iter - .5) * std::max(dS, dQ);
				//quan_temp(Non_zero[node_ref_ID]) -= 2 * (dir_iter - .5) * std::max(dS, dQ);
				
//				// Check if direction is possbile
//				if(quan_temp(Non_zero[zone_iter]) < bidded_total_aggregated(0, Non_zero[zone_iter]) - 100. * eps || quan_temp(Non_zero[zone_iter]) > bidded_total_aggregated(Market.price_intervals + 2, Non_zero[zone_iter]) + 100. * eps){
//					continue;
//				}
//				if(quan_temp(Non_zero[node_ref_ID]) < bidded_total_aggregated(0, Non_zero[node_ref_ID]) - 100. * eps || quan_temp(Non_zero[node_ref_ID]) > bidded_total_aggregated(Market.price_intervals + 2, Non_zero[node_ref_ID]) + 100. * eps){
//					continue;
//				}					
				
				// Update price after small increase / decrease
				if(quan_temp(Non_zero[zone_iter]) < bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])){
					while(quan_temp(Non_zero[zone_iter]) < bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])){
						if(price_ID_temp(Non_zero[zone_iter]) == 0){
							break;
						}
						price_ID_temp(Non_zero[zone_iter]) -= 1;
					}
				}
				else if(quan_temp(Non_zero[zone_iter]) > bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter])){
					while(quan_temp(Non_zero[zone_iter]) > bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter])){
						if(price_ID_temp(Non_zero[zone_iter]) == Market.price_intervals + 1){
							break;
						}
						price_ID_temp(Non_zero[zone_iter]) += 1;
					}
				}
				if(quan_temp(Non_zero[node_ref_ID]) < bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID])){
					while(quan_temp(Non_zero[node_ref_ID]) < bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID])){
						if(price_ID_temp(Non_zero[node_ref_ID]) == 0){
							break;
						}
						price_ID_temp(Non_zero[node_ref_ID]) -= 1;
					}
				}
				else if(quan_temp(Non_zero[node_ref_ID]) > bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID])){
					while(quan_temp(Non_zero[node_ref_ID]) > bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID])){
						if(price_ID_temp(Non_zero[node_ref_ID]) == Market.price_intervals + 1){
							break;
						}
						price_ID_temp(Non_zero[node_ref_ID]) += 1;
					}
				}
				
				// Update utility after small increase / decrease
				if(quan_temp(Non_zero[zone_iter]) < bidded_total_aggregated(0, Non_zero[zone_iter])){
					utility_temp(Non_zero[zone_iter]) = utility_aggregated(0, Non_zero[zone_iter]);
					utility_temp(Non_zero[zone_iter]) -= (quan_temp(Non_zero[zone_iter]) - bidded_total_aggregated(0, Non_zero[zone_iter])) * Market.price_range_inflex(0);
				}
				else if(quan_temp(Non_zero[zone_iter]) > bidded_total_aggregated(Market.price_intervals + 2, Non_zero[zone_iter]))	{
					utility_temp(Non_zero[zone_iter]) = utility_aggregated(Market.price_intervals + 2, Non_zero[zone_iter]);
					utility_temp(Non_zero[zone_iter]) -= (quan_temp(Non_zero[zone_iter]) - bidded_total_aggregated(Market.price_intervals + 2, Non_zero[zone_iter])) * Market.price_range_inflex(1);
				}
				else{
					if(bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]) != 0){
						utility_temp(Non_zero[zone_iter]) = (bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]) - quan_temp(Non_zero[zone_iter])) * utility_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])
							+ (quan_temp(Non_zero[zone_iter]) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])) * utility_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]);
						utility_temp(Non_zero[zone_iter]) /= bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]);	
					}
					else{
						utility_temp(Non_zero[zone_iter]) = utility_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]);
					}
				}
				if(quan_temp(Non_zero[node_ref_ID]) < bidded_total_aggregated(0, Non_zero[node_ref_ID])){
					utility_temp(Non_zero[node_ref_ID]) = utility_aggregated(0, Non_zero[node_ref_ID]);
					utility_temp(Non_zero[node_ref_ID]) -= (quan_temp(Non_zero[node_ref_ID]) - bidded_total_aggregated(0, Non_zero[node_ref_ID])) * Market.price_range_inflex(0);
				}
				else if(quan_temp(Non_zero[node_ref_ID]) > bidded_total_aggregated(Market.price_intervals + 2, Non_zero[node_ref_ID]))	{
					utility_temp(Non_zero[node_ref_ID]) = utility_aggregated(Market.price_intervals + 2, Non_zero[node_ref_ID]);
					utility_temp(Non_zero[node_ref_ID]) -= (quan_temp(Non_zero[node_ref_ID]) - bidded_total_aggregated(Market.price_intervals + 2, Non_zero[node_ref_ID])) * Market.price_range_inflex(1);
				}
				else{
					if(bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID]) - bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID]) != 0){
						utility_temp(Non_zero[node_ref_ID]) = (bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID]) - quan_temp(Non_zero[node_ref_ID])) * utility_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID])
							+ (quan_temp(Non_zero[node_ref_ID]) - bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID])) * utility_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID]);
						utility_temp(Non_zero[node_ref_ID]) /= bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID]) - bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID]);	
					}
					else{
						utility_temp(Non_zero[node_ref_ID]) = utility_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID]);
					}
				}
				
				// Update voltage and power flow after small increase / decrease
				voltage_temp.tail(Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(quan_temp.tail(Market.network.num_vertice - 1));
				flow_temp = (Problem.Constraint.eq_orig_matrix).topRightCorner(Market.network.num_edges, Market.network.num_vertice) * voltage_temp;		
		
				// Update errors
				error_temp = 0.;
				// Power flow errors
				for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
					error_temp += pow(std::min(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 0), 0.) + std::max(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 1), 0.), 2.);
				}
				// Power errors
				for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
					error_temp += source_error_weight * pow(std::min(quan_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + node_iter, 0), 0.) + std::max(quan_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + node_iter, 1), 0.), 2.);
				}
				// Voltage errors
				for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
					error_temp += pow(std::min(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 0), 0.) + std::max(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 1), 0.), 2.);
				}
				
				// Calculate objective 
				obj_temp = (1. - mu) * utility_temp.sum() - mu * error_temp;	
				
				// Check if objective should be updated
				if(obj_temp > obj && error > eps){
				//if(obj_temp > obj){
					// Update if objective has been improved
					quan(Non_zero[zone_iter]) = quan_temp(Non_zero[zone_iter]);
					quan(Non_zero[node_ref_ID]) = quan_temp(Non_zero[node_ref_ID]);
					voltage = voltage_temp;
					flow = flow_temp;
					price_ID(Non_zero[zone_iter]) = price_ID_temp(Non_zero[zone_iter]);
					price_ID(Non_zero[node_ref_ID]) = price_ID_temp(Non_zero[node_ref_ID]);
					utility(Non_zero[zone_iter]) = utility_temp(Non_zero[zone_iter]);
					utility(Non_zero[node_ref_ID]) = utility_temp(Non_zero[node_ref_ID]);
					obj = obj_temp;
					error = error_temp;	
					updated_num += 1;
					break;
				}
				else{
					// Return to origin
					quan_temp(Non_zero[zone_iter]) = quan(Non_zero[zone_iter]);
					quan_temp(Non_zero[node_ref_ID]) = quan(Non_zero[node_ref_ID]);
					price_ID_temp(Non_zero[zone_iter]) = price_ID(Non_zero[zone_iter]);
					price_ID_temp(Non_zero[node_ref_ID]) = price_ID(Non_zero[node_ref_ID]);
					utility_temp(Non_zero[zone_iter]) = utility(Non_zero[zone_iter]);
					utility_temp(Non_zero[node_ref_ID]) = utility(Non_zero[node_ref_ID]);					
				}													
			}		
		}
		
		// Check if increment should be smaller
		//if(updated_num <= Non_zero.size() / 100){
		if(updated_num == 0){
			dS /= 2.;
			//std::cout << dS << "\n";
		}
//		else if(updated_num >= Non_zero.size() / 2){
//			dS *= 2.;
//		}
		
		if(loop_count % 100 == 0){
			std::cout << loop_count << ": " << error << "; " << dS << "\n";	
		}
	}

	Problem.Solution.orig_vector.segment(Market.network.num_edges, Market.network.num_vertice) = quan;
	Problem.Solution.orig_vector.tail(Market.network.num_vertice) = voltage;
	Problem.Solution.orig_vector.head(Market.network.num_edges) = flow;	
	std::cout << "\n";
	std::cout << quan.minCoeff() << " " << quan.maxCoeff() << " " << .5 * quan.array().abs().sum() << "\n";
	std::cout << voltage.minCoeff() << " " << voltage.maxCoeff() << "\n";
	std::cout << flow.minCoeff() << " " << flow.maxCoeff() << "\n";	
	std::cout << "\n" << quan.transpose() << "\n\n";
	std::cout << "\n" << voltage.transpose() << "\n\n";
}

void Flow_Based_Market_Optimization_Test_5(int tick, market_inform &Market, LP_object &Problem){
	Eigen::MatrixXd bidded_supply = Market.submitted_supply;
	Eigen::MatrixXd bidded_demand = Market.submitted_demand;
	
	// Initial market clearing within each nodes
	Eigen::VectorXi price_ID(Market.num_zone);
	Market_clearing_nodal(tick, Market, price_ID, bidded_supply, bidded_demand);
	//std::cout << price_ID.transpose() << "\n\n";
	
	// Initialization of process variables for the main optimization loop
	Eigen::MatrixXd bidded_supply_export = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_demand_export = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_supply_import = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_demand_import = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_supply_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);
	Eigen::MatrixXd bidded_demand_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);
	Eigen::MatrixXd bidded_total_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);
	Eigen::MatrixXd utility_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);

	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
		// Supply curves
		bidded_supply_export.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter)) = Market.submitted_supply.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter));
		bidded_supply_export(price_ID(zone_iter), zone_iter) = bidded_supply(price_ID(zone_iter), zone_iter);
		bidded_supply_import.col(zone_iter).array().head(price_ID(zone_iter)) = Market.submitted_supply.col(zone_iter).array().head(price_ID(zone_iter));
		bidded_supply_import(price_ID(zone_iter), zone_iter) = Market.submitted_supply(price_ID(zone_iter), zone_iter) - bidded_supply(price_ID(zone_iter), zone_iter);
		bidded_supply_aggregated(0, zone_iter) = -bidded_supply_import.col(zone_iter).sum();
		
		// Demand curves
		bidded_demand_export.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter)) = Market.submitted_demand.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter));
		bidded_demand_export(price_ID(zone_iter), zone_iter) = Market.submitted_demand(price_ID(zone_iter), zone_iter) - bidded_demand(price_ID(zone_iter), zone_iter);
		bidded_demand_import.col(zone_iter).array().head(price_ID(zone_iter)) = Market.submitted_demand.col(zone_iter).array().head(price_ID(zone_iter));
		bidded_demand_import(price_ID(zone_iter), zone_iter) = bidded_demand(price_ID(zone_iter), zone_iter);
		bidded_demand_aggregated(0, zone_iter) = -bidded_demand_import.col(zone_iter).sum();		
	}
	
	// Aggregated import-export curves for supply and demand
	for(int price_iter = 1; price_iter < Market.price_intervals + 3; ++ price_iter){
		bidded_supply_aggregated.row(price_iter) = bidded_supply_aggregated.row(price_iter - 1) + bidded_supply_import.row(price_iter - 1) + bidded_supply_export.row(price_iter - 1);
		bidded_demand_aggregated.row(price_iter) = bidded_demand_aggregated.row(price_iter - 1) + bidded_demand_import.row(price_iter - 1) + bidded_demand_export.row(price_iter - 1);
	}
	bidded_total_aggregated = bidded_supply_aggregated + bidded_demand_aggregated;
	
	// Update inequality constraints for import / export
	Problem.Boundary.ie_orig_matrix.col(0).segment(Market.network.num_edges, Market.network.num_vertice) = bidded_total_aggregated.row(0).transpose();
	Problem.Boundary.ie_orig_matrix.col(1).segment(Market.network.num_edges, Market.network.num_vertice) = bidded_total_aggregated.row(Market.price_intervals + 2).transpose();

	// Aggregated utility function for each bidding zone
	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){	
		// Demand
		utility_aggregated(price_ID(zone_iter), zone_iter) = (bidded_demand_import(price_ID(zone_iter), zone_iter) + bidded_supply_import(price_ID(zone_iter), zone_iter)) * Market.bidded_price(price_ID(zone_iter));
		if(price_ID(zone_iter) > 0){
			for(int price_iter = price_ID(zone_iter) - 1; price_iter >= 0; -- price_iter){
				utility_aggregated(price_iter, zone_iter) = utility_aggregated(price_iter + 1, zone_iter) + (bidded_demand_import(price_iter, zone_iter) + bidded_supply_import(price_iter, zone_iter)) * Market.bidded_price(price_iter);
			}
		}
		
		// Supply
		utility_aggregated(price_ID(zone_iter) + 1, zone_iter) = -(bidded_demand_export(price_ID(zone_iter), zone_iter) + bidded_supply_export(price_ID(zone_iter), zone_iter)) * Market.bidded_price(price_ID(zone_iter));
		if(price_ID(zone_iter) < Market.price_intervals){
			for(int price_iter = price_ID(zone_iter) + 2; price_iter < Market.price_intervals + 3; ++ price_iter){
				utility_aggregated(price_iter, zone_iter) = utility_aggregated(price_iter - 1, zone_iter) - (bidded_demand_export(price_iter - 1, zone_iter) + bidded_supply_export(price_iter - 1, zone_iter)) * Market.bidded_price(price_iter - 1);
			}			
		}
	}

	// Find nodes with no power source / sink
	std::vector <int> Non_zero;
	Non_zero.reserve(Market.num_zone);
	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
		if(abs(bidded_total_aggregated(0, zone_iter) - bidded_total_aggregated(Market.price_intervals + 2, zone_iter)) > 0.){
			Non_zero.push_back(zone_iter);
		}
	}
	
	// Declare variables for the initial solver
	double eps = pow(10., -12.);
	double ratio;
	double error;
	double utility_sum;
	Eigen::Vector2d quan_boundary;
	Eigen::VectorXd quan(Market.num_zone);
	Eigen::VectorXd utility(Market.num_zone);
	Eigen::VectorXd voltage = Eigen::VectorXd::Zero(Market.num_zone);
	Eigen::VectorXd flow = Eigen::VectorXd::Zero(Market.network.num_edges);
	Eigen::MatrixXd quan_interval(Market.num_zone, 2);
	Eigen::MatrixXd utility_interval(Market.num_zone, 2);
	
	// Find optimal solution without voltage and power flow constraints
	Eigen::VectorXd bidded_total_aggregated_all = bidded_total_aggregated * Eigen::VectorXd::Constant(Market.num_zone, 1.);
	for(int price_iter = 0; price_iter < Market.price_intervals + 2; ++ price_iter){
		if(bidded_total_aggregated_all(price_iter) <= 0 && bidded_total_aggregated_all(price_iter + 1) >= 0){
			price_ID = Eigen::VectorXi::Constant(Market.num_zone, price_iter);
			if(bidded_total_aggregated_all(price_iter + 1) > bidded_total_aggregated_all(price_iter)){
				ratio = bidded_total_aggregated_all(price_iter + 1) / (bidded_total_aggregated_all(price_iter + 1) - bidded_total_aggregated_all(price_iter));				
			}
			else{
				ratio = .5;
			}
			quan = ratio * bidded_total_aggregated.row(price_iter) + (1. - ratio) * bidded_total_aggregated.row(price_iter + 1);
			quan_interval.col(0) = bidded_total_aggregated.row(price_iter);
			quan_interval.col(1) = bidded_total_aggregated.row(price_iter + 1);
			utility_interval.col(0) = utility_aggregated.row(price_iter);
			utility_interval.col(1) = utility_aggregated.row(price_iter + 1);
			quan_boundary << bidded_total_aggregated_all(price_iter), bidded_total_aggregated_all(price_iter + 1);
			
			break;
		}
	}
	
	// Find initial values for utility, voltage, and power flow
	utility = ratio * utility_aggregated.row(price_ID(0)) + (1. - ratio) * utility_aggregated.row(price_ID(0) + 1);
	utility_sum = utility.sum();
	voltage.tail(Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(quan.tail(Market.network.num_vertice - 1));
	flow = (Problem.Constraint.eq_orig_matrix).topRightCorner(Market.network.num_edges, Market.network.num_vertice) * voltage;
	
	// Find initial values for error and objective
	error = 0.;
	// Power flow errors
	for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
		error += pow(std::min(flow(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 0), 0.) + std::max(flow(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 1), 0.), 2.);
	}
	// Voltage errors
	for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
		error += pow(std::min(voltage(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 0), 0.) + std::max(voltage(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 1), 0.), 2.);
	}
	// Exit function if no constraint is violated
	if(error == 0.){
		// Update confirmed bids and prices
		return;
	}
	
	// Declare variables for the loop
	int update_count = 0;
	double ratio_temp;
	double error_temp;
	double utility_sum_temp = utility_sum;
	Eigen::VectorXi price_ID_temp = price_ID;
	Eigen::Vector2d quan_boundary_temp = quan_boundary;
	Eigen::VectorXd quan_temp = quan;
	Eigen::VectorXd utility_temp = utility;
	Eigen::VectorXd voltage_temp = voltage;
	Eigen::VectorXd flow_temp = flow;
	Eigen::MatrixXd quan_interval_temp = quan_interval;
	Eigen::MatrixXd utility_interval_temp = utility_interval;
	
	int loop_count = 0;
	std::cout << loop_count << ": " << utility_sum << "; " << error << "\n";
	while(loop_count < 500){
		loop_count += 1;
		for(int zone_iter = 0; zone_iter < Non_zero.size(); ++ zone_iter){
			for(int dir_iter = 0; dir_iter < 2; ++ dir_iter){
				quan_boundary_temp(0) -= bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]);
				quan_boundary_temp(1) -= bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]);
				price_ID_temp(Non_zero[zone_iter]) += 2 * dir_iter - 1;
				
				// Check if price_ID_temp is out of boundary
				if(price_ID_temp(Non_zero[zone_iter]) < 0 || price_ID_temp(Non_zero[zone_iter]) > Market.price_intervals + 1){
					// Price_temp out of boundary, return to origin value
					price_ID_temp(Non_zero[zone_iter]) = price_ID(Non_zero[zone_iter]);
					quan_boundary_temp = quan_boundary;
					continue;					
				}
				// Update price_ID_temp to non-zero quantity interval
				if(price_ID_temp(Non_zero[zone_iter]) + 2 * dir_iter - 1 >= 0 && price_ID_temp(Non_zero[zone_iter]) + 2 * dir_iter - 1 <= Market.price_intervals + 1){
					while(bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]) == bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 2 * dir_iter - 1, Non_zero[zone_iter])){
						price_ID_temp(Non_zero[zone_iter]) += 2 * dir_iter - 1;
						if(price_ID_temp(Non_zero[zone_iter]) + 2 * dir_iter - 1 < 0 || price_ID_temp(Non_zero[zone_iter]) + 2 * dir_iter - 1 > Market.price_intervals + 1){
							break;
						}
					}				
				}
				
				// Update quantity and utility boundary and check if 0 is within the boundaries
				quan_interval_temp.row(Non_zero[zone_iter]) << bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]), bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]);
				utility_interval_temp.row(Non_zero[zone_iter]) << utility_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]), utility_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]);
				quan_boundary_temp(0) += bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]);
				quan_boundary_temp(1) += bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]);
				if(quan_boundary_temp(0) > 0. || quan_boundary_temp(1) < 0.){
					// 0 not in quantity boundary, return to origin value
					price_ID_temp(Non_zero[zone_iter]) = price_ID(Non_zero[zone_iter]);
					quan_interval_temp.row(Non_zero[zone_iter]) = quan_interval.row(Non_zero[zone_iter]);
					utility_interval_temp.row(Non_zero[zone_iter]) = utility_interval.row(Non_zero[zone_iter]);
					quan_boundary_temp = quan_boundary;
					continue;				
				}
				
				// Find the correct ratio for quantities
				if(quan_boundary_temp(1) > quan_boundary_temp(0)){
					ratio_temp = quan_boundary_temp(1) / (quan_boundary_temp(1) - quan_boundary_temp(0));
				}
				else{
					ratio_temp = .5;
				}
				
				// Update the quantity and utility at each node
				quan_temp = ratio_temp * quan_interval_temp.col(0) + (1. - ratio_temp) * quan_interval_temp.col(1);
				utility_temp = ratio_temp * utility_interval_temp.col(0) + (1. - ratio_temp) * utility_interval_temp.col(1);
				utility_sum_temp = utility_temp.sum();
				
				// Update the voltage and power flow
				voltage_temp.tail(Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(quan_temp.tail(Market.network.num_vertice - 1));
				flow_temp = (Problem.Constraint.eq_orig_matrix).topRightCorner(Market.network.num_edges, Market.network.num_vertice) * voltage_temp;
				
				// Update errors
				error_temp = 0.;
				// Power flow errors
				for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
					error_temp += pow(std::min(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 0), 0.) + std::max(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 1), 0.), 2.);
				}
				// Voltage errors
				for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
					error_temp += pow(std::min(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 0), 0.) + std::max(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 1), 0.), 2.);
				}
				//std::cout << zone_iter << " " << dir_iter << ": " << error_temp - error << "; " << quan_temp.sum() << "\n";
				
				// Check if direction is plausible
				if(error_temp < error && error_temp > 0. && utility_sum_temp < utility_sum){
				//if(error_temp < error && error_temp > 0.){
					// Update solution if plausible
					price_ID(Non_zero[zone_iter]) = price_ID_temp(Non_zero[zone_iter]);
					quan_interval.row(Non_zero[zone_iter]) = quan_interval_temp.row(Non_zero[zone_iter]);
					utility_interval.row(Non_zero[zone_iter]) = utility_interval_temp.row(Non_zero[zone_iter]);
					ratio = ratio_temp;
					error = error_temp;
					quan_boundary = quan_boundary_temp;
					quan = quan_temp;				
					utility_sum = utility_sum_temp;
					utility = utility_temp;			
					voltage = voltage_temp;
					flow = flow_temp;
					break;								
				}
				else{
					// Return to original value
					price_ID_temp(Non_zero[zone_iter]) = price_ID(Non_zero[zone_iter]);
					quan_interval_temp.row(Non_zero[zone_iter]) = quan_interval.row(Non_zero[zone_iter]);
					utility_interval_temp.row(Non_zero[zone_iter]) = utility_interval.row(Non_zero[zone_iter]);
					quan_boundary_temp = quan_boundary;
					quan_temp = quan;
					utility_temp = utility;
				}	
			}
		}
		
		if(loop_count % 10 == 0){
			std::cout << loop_count << ": " << utility_sum << "; " << error << "\n";
		}		
	}
	
	std::cout << "\n";
	std::cout << std::setprecision(6) << quan.transpose() << "\n\n";
	std::cout << price_ID.transpose() << "\n\n";	
}

void Flow_Based_Market_Optimization(int tick, market_inform &Market, LP_object &Problem){
	Eigen::MatrixXd bidded_supply = Market.submitted_supply;
	Eigen::MatrixXd bidded_demand = Market.submitted_demand;
	
	// Initial market clearing within each nodes
	Eigen::VectorXi price_ID(Market.num_zone);
	Market_clearing_nodal(tick, Market, price_ID, bidded_supply, bidded_demand);
	//std::cout << price_ID.transpose() << "\n\n";
	
	// Initialization of process variables for the main optimization loop
	Eigen::MatrixXd bidded_supply_export = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_demand_export = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_supply_import = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_demand_import = Eigen::MatrixXd::Zero(Market.price_intervals + 2, Market.num_zone);
	Eigen::MatrixXd bidded_supply_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);
	Eigen::MatrixXd bidded_demand_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);
	Eigen::MatrixXd bidded_total_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);
	Eigen::MatrixXd utility_aggregated = Eigen::MatrixXd::Zero(Market.price_intervals + 3, Market.num_zone);

	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
		// Supply curves
		bidded_supply_export.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter)) = Market.submitted_supply.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter));
		bidded_supply_export(price_ID(zone_iter), zone_iter) = bidded_supply(price_ID(zone_iter), zone_iter);
		bidded_supply_import.col(zone_iter).array().head(price_ID(zone_iter)) = Market.submitted_supply.col(zone_iter).array().head(price_ID(zone_iter));
		bidded_supply_import(price_ID(zone_iter), zone_iter) = Market.submitted_supply(price_ID(zone_iter), zone_iter) - bidded_supply(price_ID(zone_iter), zone_iter);
		bidded_supply_aggregated(0, zone_iter) = -bidded_supply_import.col(zone_iter).sum();
		
		// Demand curves
		bidded_demand_export.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter)) = Market.submitted_demand.col(zone_iter).array().tail(Market.price_intervals + 1 - price_ID(zone_iter));
		bidded_demand_export(price_ID(zone_iter), zone_iter) = Market.submitted_demand(price_ID(zone_iter), zone_iter) - bidded_demand(price_ID(zone_iter), zone_iter);
		bidded_demand_import.col(zone_iter).array().head(price_ID(zone_iter)) = Market.submitted_demand.col(zone_iter).array().head(price_ID(zone_iter));
		bidded_demand_import(price_ID(zone_iter), zone_iter) = bidded_demand(price_ID(zone_iter), zone_iter);
		bidded_demand_aggregated(0, zone_iter) = -bidded_demand_import.col(zone_iter).sum();		
	}
	
	// Aggregated import-export curves for supply and demand
	for(int price_iter = 1; price_iter < Market.price_intervals + 3; ++ price_iter){
		bidded_supply_aggregated.row(price_iter) = bidded_supply_aggregated.row(price_iter - 1) + bidded_supply_import.row(price_iter - 1) + bidded_supply_export.row(price_iter - 1);
		bidded_demand_aggregated.row(price_iter) = bidded_demand_aggregated.row(price_iter - 1) + bidded_demand_import.row(price_iter - 1) + bidded_demand_export.row(price_iter - 1);
	}
	bidded_total_aggregated = bidded_supply_aggregated + bidded_demand_aggregated;
	
	// Update inequality constraints for import / export
	Problem.Boundary.ie_orig_matrix.col(0).segment(Market.network.num_edges, Market.network.num_vertice) = bidded_total_aggregated.row(0).transpose();
	Problem.Boundary.ie_orig_matrix.col(1).segment(Market.network.num_edges, Market.network.num_vertice) = bidded_total_aggregated.row(Market.price_intervals + 2).transpose();

	// Aggregated utility function for each bidding zone
	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){	
		// Demand
		utility_aggregated(price_ID(zone_iter), zone_iter) = (bidded_demand_import(price_ID(zone_iter), zone_iter) + bidded_supply_import(price_ID(zone_iter), zone_iter)) * Market.bidded_price(price_ID(zone_iter));
		if(price_ID(zone_iter) > 0){
			for(int price_iter = price_ID(zone_iter) - 1; price_iter >= 0; -- price_iter){
				utility_aggregated(price_iter, zone_iter) = utility_aggregated(price_iter + 1, zone_iter) + (bidded_demand_import(price_iter, zone_iter) + bidded_supply_import(price_iter, zone_iter)) * Market.bidded_price(price_iter);
			}
		}
		
		// Supply
		utility_aggregated(price_ID(zone_iter) + 1, zone_iter) = -(bidded_demand_export(price_ID(zone_iter), zone_iter) + bidded_supply_export(price_ID(zone_iter), zone_iter)) * Market.bidded_price(price_ID(zone_iter));
		if(price_ID(zone_iter) < Market.price_intervals){
			for(int price_iter = price_ID(zone_iter) + 2; price_iter < Market.price_intervals + 3; ++ price_iter){
				utility_aggregated(price_iter, zone_iter) = utility_aggregated(price_iter - 1, zone_iter) - (bidded_demand_export(price_iter - 1, zone_iter) + bidded_supply_export(price_iter - 1, zone_iter)) * Market.bidded_price(price_iter - 1);
			}			
		}
	}

	// Find nodes with no power source / sink
	Eigen::VectorXi Zero_node = Eigen::VectorXi::Zero(Market.num_zone);
	std::vector <int> Non_zero;
	Non_zero.reserve(Market.num_zone);
	for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
		if(abs(bidded_total_aggregated(0, zone_iter) - bidded_total_aggregated(Market.price_intervals + 2, zone_iter)) <= 0.){
			Zero_node(zone_iter) = 1;
		}
		else{
			Non_zero.push_back(zone_iter);
		}
	}
	
	// Declare variables for the main loop
	int updated_num;
	int node_ref_ID;
	int zone_iter;
	//double tol = pow(10., -12.);
	double eps = pow(10., -10.);
	double mu = 1 - eps;
	double dS_max = .2 * (bidded_total_aggregated.bottomRows(1)).lpNorm <1> () / bidded_total_aggregated.cols();
	//double dS_max = 10.;
	double dS = dS_max;
	double error;
	double obj = 0.;
	double obj_temp;
	Eigen::VectorXi price_ID_temp = price_ID;
	Eigen::VectorXd voltage = Eigen::VectorXd::Zero(Market.num_zone);
	Eigen::VectorXd voltage_temp = voltage;
	Eigen::VectorXd quan = Eigen::VectorXd::Zero(Market.num_zone);	
	Eigen::VectorXd quan_temp = quan;	
	Eigen::VectorXd flow = Eigen::VectorXd::Zero(Market.network.num_edges);
	Eigen::VectorXd flow_temp = flow;
	Eigen::VectorXd utility = Eigen::VectorXd::Zero(Market.num_zone);
	Eigen::VectorXd utility_temp = utility;
	//Eigen::VectorXd grad = Eigen::VectorXd::Zero(Market.num_zone);
	
	int loop_count = 0;
	//while(loop_count < 200){
	while(dS > .01 * dS_max){
		loop_count += 1;
		updated_num = 0;
		//for(int zone_iter = 0; zone_iter < Market.num_zone; ++ zone_iter){
		for(int rand_draw = 0; rand_draw < 100000; ++ rand_draw){
			zone_iter = std::rand() % Non_zero.size();
			node_ref_ID = std::rand() % Non_zero.size();
			if(zone_iter == node_ref_ID){
				continue;
			}			
			
			for(int dir_iter = 0; dir_iter < 2; ++ dir_iter){
				//Check if direction is plausible
				if(price_ID_temp(Non_zero[zone_iter]) < price_ID_temp(Non_zero[node_ref_ID]) && dir_iter == 0){
					continue;
				}
				else if(price_ID_temp(Non_zero[zone_iter]) > price_ID_temp(Non_zero[node_ref_ID]) && dir_iter == 1){
					continue;
				}
				
				// Update quantity after small increase / decrease
				quan_temp(Non_zero[zone_iter]) += 2 * (dir_iter - .5) * dS;
				quan_temp(Non_zero[node_ref_ID]) -= 2 * (dir_iter - .5) * dS;
				
				// Update price after small increase / decrease
				if(quan_temp(Non_zero[zone_iter]) < bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])){
					while(quan_temp(Non_zero[zone_iter]) < bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])){
						if(price_ID_temp(Non_zero[zone_iter]) == 0){
							break;
						}
						price_ID_temp(Non_zero[zone_iter]) -= 1;
					}
				}
				else if(quan_temp(Non_zero[zone_iter]) > bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter])){
					while(quan_temp(Non_zero[zone_iter]) > bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter])){
						if(price_ID_temp(Non_zero[zone_iter]) == Market.price_intervals + 1){
							break;
						}
						price_ID_temp(Non_zero[zone_iter]) += 1;
					}
				}
				if(quan_temp(Non_zero[node_ref_ID]) < bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID])){
					while(quan_temp(Non_zero[node_ref_ID]) < bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID])){
						if(price_ID_temp(Non_zero[node_ref_ID]) == 0){
							break;
						}
						price_ID_temp(Non_zero[node_ref_ID]) -= 1;
					}
				}
				else if(quan_temp(Non_zero[node_ref_ID]) > bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID])){
					while(quan_temp(Non_zero[node_ref_ID]) > bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID])){
						if(price_ID_temp(Non_zero[node_ref_ID]) == Market.price_intervals + 1){
							break;
						}
						price_ID_temp(Non_zero[node_ref_ID]) += 1;
					}
				}				
				
				// Update utility after small increase / decrease
				if(quan_temp(Non_zero[zone_iter]) < bidded_total_aggregated(0, Non_zero[zone_iter])){
					utility_temp(Non_zero[zone_iter]) = utility_aggregated(0, Non_zero[zone_iter]);
					utility_temp(Non_zero[zone_iter]) -= (quan_temp(Non_zero[zone_iter]) - bidded_total_aggregated(0, Non_zero[zone_iter])) * Market.price_range_inflex(0);
				}
				else if(quan_temp(Non_zero[zone_iter]) > bidded_total_aggregated(Market.price_intervals + 2, Non_zero[zone_iter]))	{
					utility_temp(Non_zero[zone_iter]) = utility_aggregated(Market.price_intervals + 2, Non_zero[zone_iter]);
					utility_temp(Non_zero[zone_iter]) -= (quan_temp(Non_zero[zone_iter]) - bidded_total_aggregated(Market.price_intervals + 2, Non_zero[zone_iter])) * Market.price_range_inflex(1);
				}
				else{
					if(bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]) != 0){
						utility_temp(Non_zero[zone_iter]) = (bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]) - quan_temp(Non_zero[zone_iter])) * utility_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])
							+ (quan_temp(Non_zero[zone_iter]) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter])) * utility_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]);
						utility_temp(Non_zero[zone_iter]) /= bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]) + 1, Non_zero[zone_iter]) - bidded_total_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]);	
					}
					else{
						utility_temp(Non_zero[zone_iter]) = utility_aggregated(price_ID_temp(Non_zero[zone_iter]), Non_zero[zone_iter]);
					}
				}
				if(quan_temp(Non_zero[node_ref_ID]) < bidded_total_aggregated(0, Non_zero[node_ref_ID])){
					utility_temp(Non_zero[node_ref_ID]) = utility_aggregated(0, Non_zero[node_ref_ID]);
					utility_temp(Non_zero[node_ref_ID]) -= (quan_temp(Non_zero[node_ref_ID]) - bidded_total_aggregated(0, Non_zero[node_ref_ID])) * Market.price_range_inflex(0);
				}
				else if(quan_temp(Non_zero[node_ref_ID]) > bidded_total_aggregated(Market.price_intervals + 2, Non_zero[node_ref_ID]))	{
					utility_temp(Non_zero[node_ref_ID]) = utility_aggregated(Market.price_intervals + 2, Non_zero[node_ref_ID]);
					utility_temp(Non_zero[node_ref_ID]) -= (quan_temp(Non_zero[node_ref_ID]) - bidded_total_aggregated(Market.price_intervals + 2, Non_zero[node_ref_ID])) * Market.price_range_inflex(1);
				}
				else{
					if(bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID]) - bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID]) != 0){
						utility_temp(Non_zero[node_ref_ID]) = (bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID]) - quan_temp(Non_zero[node_ref_ID])) * utility_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID])
							+ (quan_temp(Non_zero[node_ref_ID]) - bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID])) * utility_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID]);
						utility_temp(Non_zero[node_ref_ID]) /= bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]) + 1, Non_zero[node_ref_ID]) - bidded_total_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID]);	
					}
					else{
						utility_temp(Non_zero[node_ref_ID]) = utility_aggregated(price_ID_temp(Non_zero[node_ref_ID]), Non_zero[node_ref_ID]);
					}
				}
				
				// Update voltage and power flow after small increase / decrease
				voltage_temp.tail(Market.network.num_vertice - 1) = Problem.Solver.ldlt.solve(quan_temp.tail(Market.network.num_vertice - 1));
				flow_temp = (Problem.Constraint.eq_orig_matrix).topRightCorner(Market.network.num_edges, Market.network.num_vertice) * voltage_temp;		
		
				// Update errors
				error = 0.;
				// Power flow errors
				for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
					error += pow(std::min(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 0), 0.) + std::max(flow_temp(edge_iter) - Problem.Boundary.ie_orig_matrix(edge_iter, 1), 0.), 2.);
				}
				// Power errors
				for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
					error += pow(std::min(quan_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + node_iter, 0), 0.) + std::max(quan_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + node_iter, 1), 0.), 2.);
				}
				// Voltage errors
				for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
					error += pow(std::min(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 0), 0.) + std::max(voltage_temp(node_iter) - Problem.Boundary.ie_orig_matrix(Market.network.num_edges + Market.network.num_vertice + node_iter, 1), 0.), 2.);
				}
				
				// Calculate objective 
				obj_temp = (1. - mu) * utility_temp.sum() - mu * error;
				//std::cout << Non_zero[zone_iter] << " " << Non_zero[node_ref_ID] << " " << obj_temp << "\n";
				
				// Check if objective should be updated
				if(obj_temp > obj){
					// Update if objective has been improved
					quan(Non_zero[zone_iter]) = quan_temp(Non_zero[zone_iter]);
					quan(Non_zero[node_ref_ID]) = quan_temp(Non_zero[node_ref_ID]);
					voltage = voltage_temp;
					flow = flow_temp;
					price_ID(Non_zero[zone_iter]) = price_ID_temp(Non_zero[zone_iter]);
					price_ID(Non_zero[node_ref_ID]) = price_ID_temp(Non_zero[node_ref_ID]);
					utility(Non_zero[zone_iter]) = utility_temp(Non_zero[zone_iter]);
					utility(Non_zero[node_ref_ID]) = utility_temp(Non_zero[node_ref_ID]);
					obj = obj_temp;	
					updated_num += 1;
					break;
				}
				else{
					// Return to origin
					quan_temp(Non_zero[zone_iter]) = quan(Non_zero[zone_iter]);
					quan_temp(Non_zero[node_ref_ID]) = quan(Non_zero[node_ref_ID]);
					price_ID_temp(Non_zero[zone_iter]) = price_ID(Non_zero[zone_iter]);
					price_ID_temp(Non_zero[node_ref_ID]) = price_ID(Non_zero[node_ref_ID]);
					utility_temp(Non_zero[zone_iter]) = utility(Non_zero[zone_iter]);
					utility_temp(Non_zero[node_ref_ID]) = utility(Non_zero[node_ref_ID]);					
				}																					
			}	
		}
		
		// Check if increment should be smaller
		if(updated_num == 0){
			dS /= 2.;
			//std::cout << dS << "\n";
		}
//		else if(dS < .1){
//			dS *= 2;
//		}
		
		std::cout << loop_count << ": " << quan.minCoeff() << " " << quan.maxCoeff() << " " << .5 * quan.array().abs().sum() << " " << obj << "\n";
	}

	Problem.Solution.orig_vector.segment(Market.network.num_edges, Market.network.num_vertice) = quan;
	Problem.Solution.orig_vector.tail(Market.network.num_vertice) = voltage;
	Problem.Solution.orig_vector.head(Market.network.num_edges) = flow;	
	std::cout << quan.minCoeff() << " " << quan.maxCoeff() << " " << .5 * quan.array().abs().sum() << "\n";
	std::cout << voltage.minCoeff() << " " << voltage.maxCoeff() << "\n";
	std::cout << flow.minCoeff() << " " << flow.maxCoeff() << "\n";	
	//std::cout << "\n" << grad.transpose() << "\n\n";
	//std::cout << "\n" << obj << "\n\n";
	//std::cout << "\n" << quan.transpose() << "\n\n";
	//std::cout << "\n" << voltage.transpose() << "\n\n";
	//std::cout << "\n" << flow.transpose() << "\n\n";
}