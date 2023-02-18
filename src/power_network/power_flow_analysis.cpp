// Source file for power flow analysis
#include "power_flow_analysis.h"

// HELM method
// P-Q buses
// Set {V}(s) = \sum {a}_n * s^n
// {\hat V}(s) = \sum {b}_n * s^n
// {1. / V}(s) = \sum {c}_n * s^n
// {1. / \hat V}(s) = \sum {d}_n * s^n
// Equations:
// [Y] {V}(s) = s * {Conj(S)} . {1. / \hat V}(s)
// [Conj(Y)] {\hat V}(s) = s * {S} . {1. / V}(s)
//
// P-U buses
// Set {a}, {b}, {c}, {d} the same as P-Q bus
// {Q}(s) = \sum {e}_n * s^n
// Equations:
// [Y] {V}(s) = (s * {P} - j {Q}(s)) . {1. / \hat V}(s)
// [Conj(Y)] {\hat V}(s) = (s * {P} + j {Q}(s)) . {1. / V}(s)
// {V}(s) * {\hat V}(s) = {U_0^2} + s * ({U^2} - {U_0^2})
//
// Conservation law for currents
// {S} . {1. / V}(s)  = 0.
// {Conj(S)} . {1. / \hat V}(s) = 0.

namespace local{
	Eigen::VectorXcd Pade(Eigen::VectorXcd c){
		int n = c.size() / 2;
		int m = c.size() - n;

		// Find the Denominator Terms (Normalized)
		Eigen::MatrixXcd M_D = Eigen::MatrixXcd::Zero(n, n);
		for(int row_ID = 0; row_ID < n; row_ID ++){
			M_D.row(row_ID) = c.segment(m + row_ID - n, n).reverse();
		}
		Eigen::VectorXcd b(n + 1);
		b << 1, M_D.colPivHouseholderQr().solve(-c.segment(m, n));
		b /= b.norm();		// Normalization to avoid singularity and ill-conditioning

		// Find the Nominator Terms
		Eigen::VectorXcd a = Eigen::VectorXcd::Zero(m);
		for(int item_ID = 0; item_ID < m; item_ID ++){
			for(int col_ID = 0; col_ID <= item_ID; col_ID ++){
			  a(item_ID) += c(item_ID - col_ID) * b(col_ID);
			}
		}

		Eigen::VectorXcd result(m + n + 1);
		result << a, b;

		return result;
	}
}

void power_network::HELM_Transmission_Set(network_inform &Power_network_inform, power_market::market_whole_inform &Power_market_inform){
	int node_num = Power_network_inform.nodes.bidding_zone.size();
	int edge_num = Power_network_inform.edges.distance.size();

	// -------------------------------------------------------------------------------
	// Set the nodal admittance matrix
	// -------------------------------------------------------------------------------
	std::vector<Eigen::TripletXcd> Y_n_trip;
	Y_n_trip.reserve(node_num + 2 * edge_num);
	Eigen::VectorXcd Y_n_Diag = Eigen::VectorXcd::Zero(node_num);
	Power_market_inform.TSO_Market.power_flow.edge_admittance = Eigen::VectorXcd(edge_num);

	for(int edge_iter = 0; edge_iter < edge_num; ++ edge_iter){
		int from_ID = Power_network_inform.edges.from(edge_iter);
		int to_ID = Power_network_inform.edges.to(edge_iter);
		int voltage = Power_network_inform.edges.voltage_base(edge_iter);

		// Series admittance
		std::complex <double> y_series(1., 0.);
		y_series /= Power_network_inform.edges.distance(edge_iter);
		y_series /= Power_network_inform.tech_parameters.z_trans_series;
		y_series *= Power_network_inform.tech_parameters.impedenace_base_levels[voltage];
		Power_market_inform.TSO_Market.power_flow.edge_admittance(edge_iter) = y_series;

		// Shunt admittance
		std::complex <double> y_shunt(1., 0.);
		y_shunt *= Power_network_inform.edges.distance(edge_iter);
		y_shunt *= Power_network_inform.tech_parameters.y_trans_shunt;
		y_shunt *= Power_network_inform.tech_parameters.impedenace_base_levels[voltage];

		// Triplet for series impedence
		Y_n_trip.push_back(Eigen::TripletXcd(from_ID, to_ID, -y_series));
		Y_n_trip.push_back(Eigen::TripletXcd(to_ID, from_ID, -y_series));

		// Update diagonal terms
		Y_n_Diag(from_ID) += y_series + .5 * y_shunt;
		Y_n_Diag(to_ID) += y_series + .5 * y_shunt;
	}
	//std::cout << Power_market_inform.TSO_Market.power_flow.edge_admittance << "\n\n";

	// Triplet for diagonal terms
	for(int node_iter = 0; node_iter < node_num; ++ node_iter){
		Y_n_trip.push_back(Eigen::TripletXcd(node_iter, node_iter, Y_n_Diag(node_iter)));
	}

	// Store the nodal admittance matrix
	Power_market_inform.TSO_Market.power_flow.nodal_admittance = Eigen::SparseMatrix <std::complex <double>> (node_num, node_num);
	Power_market_inform.TSO_Market.power_flow.nodal_admittance.setFromTriplets(Y_n_trip.begin(), Y_n_trip.end());

	// -------------------------------------------------------------------------------
	// Determine type of buses
	// -------------------------------------------------------------------------------
	// Consider only HV power suppliers
	Power_market_inform.TSO_Market.power_flow.PQ_bus.reserve(node_num);
	Power_market_inform.TSO_Market.power_flow.PU_bus.reserve(node_num);
	Power_market_inform.TSO_Market.power_flow.ref_bus.reserve(node_num);
	Eigen::VectorXi bus_type = Eigen::VectorXi::Zero(node_num);
	int node_ref_num;
	for(int edge_iter = 0; edge_iter < Power_network_inform.cbt.entry_bz.size(); ++ edge_iter){
		if(Power_network_inform.cbt.entry_node_num(edge_iter) == 0){
			continue;
		}
		node_ref_num = (int) Power_network_inform.cbt.entry_nodes(edge_iter, 0);
		break;
	}
	bus_type(node_ref_num) = 2;

	int hydro_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		if(bus_type(node_ID) != 0){
			continue;
		}

		bus_type(node_ID) = 1;
	}

	int wind_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		if(bus_type(node_ID) != 0){
			continue;
		}

		bus_type(node_ID) = 1;
	}

	int pump_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV.size();
	for(int agent_iter = 0; agent_iter < pump_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		if(bus_type(node_ID) != 0){
			continue;
		}

		bus_type(node_ID) = 1;
	}

	for(int bus_iter = 0; bus_iter < node_num; ++ bus_iter){
		if(bus_type(bus_iter) == 0){
			Power_market_inform.TSO_Market.power_flow.PQ_bus.push_back(bus_iter);
		}
		else if(bus_type(bus_iter) == 1){
			Power_market_inform.TSO_Market.power_flow.PU_bus.push_back(bus_iter);
		}
		else{
			Power_market_inform.TSO_Market.power_flow.ref_bus.push_back(bus_iter);
		}
	}

	// -------------------------------------------------------------------------------
	// Set the solver for iterative process
	// -------------------------------------------------------------------------------
	// Entries from original nodal admittance matrix, reordered and delete reference bus
	std::vector<Eigen::TripletXcd> Permut_Y_n_trip;
	Permut_Y_n_trip.reserve(node_num);
	for(int bus_iter = 0; bus_iter < Power_market_inform.TSO_Market.power_flow.PU_bus.size(); ++ bus_iter){
		int row_ID = bus_iter;
		int bus_ID = Power_market_inform.TSO_Market.power_flow.PU_bus[bus_iter];
		Permut_Y_n_trip.push_back(Eigen::TripletXcd(row_ID, bus_ID, 1.));
	}
	for(int bus_iter = 0; bus_iter < Power_market_inform.TSO_Market.power_flow.PQ_bus.size(); ++ bus_iter){
		int row_ID = Power_market_inform.TSO_Market.power_flow.PU_bus.size() + bus_iter;
		int bus_ID = Power_market_inform.TSO_Market.power_flow.PQ_bus[bus_iter];
		Permut_Y_n_trip.push_back(Eigen::TripletXcd(row_ID, bus_ID, 1.));
	}
	Eigen::SparseMatrix <std::complex <double>> Permut_Y_n(Permut_Y_n_trip.size(), node_num);
	Permut_Y_n.setFromTriplets(Permut_Y_n_trip.begin(), Permut_Y_n_trip.end());
	Eigen::SparseMatrix <std::complex <double>> Y_n_small = Permut_Y_n * Power_market_inform.TSO_Market.power_flow.nodal_admittance * Permut_Y_n.transpose();

	// Matrix for solver
	std::vector<Eigen::TripletXcd> Mat_trip;
	Mat_trip.reserve(2 * Y_n_small.nonZeros() + 4 * Power_market_inform.TSO_Market.power_flow.PU_bus.size() + 4 * (node_num - node_ref_num));

	// Entries from reduced nodal admittance matrix
	for(int col_iter = 0; col_iter < Y_n_small.outerSize(); ++ col_iter){
		for(Eigen::SparseMatrix<std::complex <double>>::InnerIterator inner_iter(Y_n_small, col_iter); inner_iter; ++ inner_iter){
			std::complex <double> y_conj = inner_iter.value();
			y_conj = std::conj(y_conj);

			Mat_trip.push_back(Eigen::TripletXcd(inner_iter.row(), inner_iter.col(), inner_iter.value()));
			Mat_trip.push_back(Eigen::TripletXcd(2 * Y_n_small.rows() + inner_iter.row(), 2 * Y_n_small.rows() + inner_iter.col(), y_conj));
		}
	}

	// Reciporal relation for V and 1/V
	for(int bus_iter = 0; bus_iter < Y_n_small.rows(); ++ bus_iter){
		int row_ID = Y_n_small.rows() + bus_iter;
		int col_ID_1 = bus_iter;
		int col_ID_2 = Y_n_small.rows() + bus_iter;
		Mat_trip.push_back(Eigen::TripletXcd(row_ID, col_ID_1, 1.));
		Mat_trip.push_back(Eigen::TripletXcd(row_ID, col_ID_2, 1.));

		row_ID += 2 * Y_n_small.rows();
		col_ID_1 += 2 * Y_n_small.rows();
		col_ID_2 += 2 * Y_n_small.rows();
		Mat_trip.push_back(Eigen::TripletXcd(row_ID, col_ID_1, 1.));
		Mat_trip.push_back(Eigen::TripletXcd(row_ID, col_ID_2, 1.));
	}

	// Additional terms for PU buses
	for(int node_iter = 0; node_iter < Power_market_inform.TSO_Market.power_flow.PU_bus.size(); ++ node_iter){
		// Entries from reactive power at PU buses
		int col_ID = 4 * Y_n_small.rows() + node_iter;
		std::complex <double> root_i(0., 1.);

		Mat_trip.push_back(Eigen::TripletXcd(node_iter, col_ID, root_i));
		Mat_trip.push_back(Eigen::TripletXcd(node_iter + 2 * Y_n_small.rows(), col_ID, -root_i));

		// Entries for voltage magnitude constraint
		// Assume trivial solution is V = 1. everywhere
		int row_ID = 4 * Y_n_small.rows() + node_iter;
		Mat_trip.push_back(Eigen::TripletXcd(row_ID, node_iter, 1.));
		Mat_trip.push_back(Eigen::TripletXcd(row_ID, node_iter + 2 * Y_n_small.rows(), 1.));
	}
	Eigen::SparseMatrix <std::complex <double>> Mat(4 * Y_n_small.rows() + Power_market_inform.TSO_Market.power_flow.PU_bus.size(), 4 * Y_n_small.rows() + Power_market_inform.TSO_Market.power_flow.PU_bus.size());
	Mat.setFromTriplets(Mat_trip.begin(), Mat_trip.end());
//	std::cout << Mat.row(0) << "\n\n";
//	std::cout << Mat.row(220) << "\n\n";
//	std::cout << Mat.row(440) << "\n\n";
//	std::cout << Mat.row(660) << "\n\n";
//	std::cout << Mat.row(880) << "\n\n";
	Power_market_inform.TSO_Market.power_flow.solver.compute(Mat);
	//std::cout << Mat.rows() << "\t" << Power_market_inform.TSO_Market.power_flow.solver.rank() << "\n\n";

	// -------------------------------------------------------------------------------
	// Initialize
	// -------------------------------------------------------------------------------
	int Time = configuration::parameters::Time();

	Power_market_inform.TSO_Market.power_flow.P_node = Eigen::MatrixXd::Zero(Time, node_num);
	Power_market_inform.TSO_Market.power_flow.Q_node = Eigen::MatrixXd::Zero(Time, node_num);
	Power_market_inform.TSO_Market.power_flow.voltage_abs = Eigen::MatrixXd::Zero(Time, node_num);
	Power_market_inform.TSO_Market.power_flow.voltage_arg = Eigen::MatrixXd::Zero(Time, node_num);
	Power_market_inform.TSO_Market.power_flow.current_abs = Eigen::MatrixXd::Zero(Time, edge_num);
	Power_market_inform.TSO_Market.power_flow.current_arg = Eigen::MatrixXd::Zero(Time, edge_num);
}

void power_network::HELM_Transmission_Update(int tick, network_inform &Power_network_inform, power_market::market_whole_inform &Power_market_inform){
	int node_num = Power_network_inform.nodes.bidding_zone.size();
	int point_num = Power_network_inform.points.bidding_zone.size();

	// Update from cross-border flows
	int edge_num = Power_market_inform.agent_profiles.cross_border.size();
	for(int edge_iter = 0; edge_iter < edge_num; ++ edge_iter){
		int node_num = Power_market_inform.agent_profiles.cross_border[edge_iter].node_num;
		int entry_bz_ID = Power_market_inform.agent_profiles.cross_border[edge_iter].entry_bz_ID;

		if(node_num == 0){
			continue;
		}
		for(int node_iter = 0; node_iter < node_num; ++ node_iter){
			int node_ID = Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].node_ID;
			double real_power;
			real_power = Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results.actual_supply;
			real_power -= Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results.actual_demand;
			Power_market_inform.TSO_Market.power_flow.P_node(tick, node_ID) += real_power;
		}
	}

	// Update from end-users
	{
		//int sample_num = agent::end_user::parameters::sample_num();
		int sample_num = Power_market_inform.agent_profiles.end_user_type.cols();
		double apparent_power_base = pow(1. - pow(agent::end_user::parameters::power_factor(), 2.), .5);
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			int node_ID = Power_network_inform.points.node(point_iter);

			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				double real_power;
				real_power = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_supply;
				Power_market_inform.TSO_Market.power_flow.P_node(tick, node_ID) += real_power;
				real_power = -Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_demand;
				Power_market_inform.TSO_Market.power_flow.P_node(tick, node_ID) += real_power;
				Power_market_inform.TSO_Market.power_flow.Q_node(tick, node_ID) += real_power * apparent_power_base;
			}
		}
	}

	// Update from industrial demand
	{
		int industrial_HV_num = Power_market_inform.agent_profiles.industrial.HV.size();
		double apparent_power_base = pow(1. - pow(agent::end_user::parameters::power_factor(), 2.), .5);
		for(int agent_iter = 0; agent_iter < industrial_HV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.industrial.HV[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);

			double real_power = -Power_market_inform.agent_profiles.industrial.HV[agent_iter].results.actual_demand;
			Power_market_inform.TSO_Market.power_flow.P_node(tick, node_ID) += real_power;
			Power_market_inform.TSO_Market.power_flow.Q_node(tick, node_ID) += real_power * apparent_power_base;
		}
	}

	// Update from power suppliers
	int hydro_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		double real_power = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		Power_market_inform.TSO_Market.power_flow.P_node(tick, node_ID) += real_power;
	}

	int hydro_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		double real_power = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		Power_market_inform.TSO_Market.power_flow.P_node(tick, node_ID) += real_power;
	}

	int wind_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		double real_power = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		Power_market_inform.TSO_Market.power_flow.P_node(tick, node_ID) += real_power;
	}

	int wind_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		double real_power = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		Power_market_inform.TSO_Market.power_flow.P_node(tick, node_ID) += real_power;
	}

	int pump_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV.size();
	for(int agent_iter = 0; agent_iter < pump_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		double real_power = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		real_power -= Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].results.actual_demand;
		Power_market_inform.TSO_Market.power_flow.P_node(tick, node_ID) += real_power;
	}

	int pump_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
	for(int agent_iter = 0; agent_iter < pump_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		double real_power = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		real_power -= Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].results.actual_demand;
		Power_market_inform.TSO_Market.power_flow.P_node(tick, node_ID) += real_power;
	}
}

void power_network::HELM_Transmission_Solve(int tick, network_inform& Power_network_inform, power_market::market_whole_inform& Power_market_inform){
	int node_num = Power_network_inform.nodes.bidding_zone.size();
	int PU_bus_num = Power_market_inform.TSO_Market.power_flow.PU_bus.size();
	int PQ_bus_num = Power_market_inform.TSO_Market.power_flow.PQ_bus.size();
	int ref_bus_num = Power_market_inform.TSO_Market.power_flow.ref_bus.size();
	int node_small_num = node_num - ref_bus_num;
	std::complex <double> root_i(0., 1.);

//	double test_factor = 1.;
//	Power_market_inform.TSO_Market.power_flow.P_node.row(tick) *= test_factor;
//	Power_market_inform.TSO_Market.power_flow.Q_node.row(tick) *= test_factor;
	//std::cout << Power_market_inform.TSO_Market.power_flow.P_node.row(tick) << "\n\n";

	// -------------------------------------------------------------------------------
	// Initialization of power series coefficients
	// -------------------------------------------------------------------------------
	int power_terms = 200;
	Eigen::MatrixXcd V_up_reg = Eigen::MatrixXcd::Zero(node_small_num, power_terms);
	Eigen::MatrixXcd V_up_hat = Eigen::MatrixXcd::Zero(node_small_num, power_terms);
	Eigen::MatrixXcd V_down_reg = Eigen::MatrixXcd::Zero(node_small_num, power_terms);
	Eigen::MatrixXcd V_down_hat = Eigen::MatrixXcd::Zero(node_small_num, power_terms);
	Eigen::MatrixXcd Q_node = Eigen::MatrixXcd::Zero(PU_bus_num, power_terms);

	// -------------------------------------------------------------------------------
	// Trivial solution (V = 1 and S = 0 everywhere)
	// -------------------------------------------------------------------------------
	V_up_reg.col(0) = Eigen::VectorXcd::Ones(node_small_num);
	V_up_hat.col(0) = Eigen::VectorXcd::Ones(node_small_num);
	V_down_reg.col(0) = Eigen::VectorXcd::Ones(node_small_num);
	V_down_hat.col(0) = Eigen::VectorXcd::Ones(node_small_num);

	// -------------------------------------------------------------------------------
	// Iteratively solve the linear equations
	// -------------------------------------------------------------------------------
	for(int term_iter = 1; term_iter < power_terms; ++ term_iter){
		//std::cout << term_iter << "\n\n";
		Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(4 * node_small_num + PU_bus_num);

		// Rhs for PU Bus
		for(int bus_iter = 0; bus_iter < PU_bus_num; ++ bus_iter){
			int row_ID = bus_iter;
			int bus_ID = Power_market_inform.TSO_Market.power_flow.PU_bus[bus_iter];

			double P_node = Power_market_inform.TSO_Market.power_flow.P_node(tick, bus_ID);
			rhs(row_ID) = P_node * V_down_hat(bus_iter , term_iter - 1);
			for(int term_iter_2 = 1; term_iter_2 < term_iter; ++ term_iter_2){
				rhs(row_ID) += -root_i * Q_node(bus_iter , term_iter_2) * V_down_hat(bus_iter , term_iter - term_iter_2);
			}

			row_ID += 2 * node_small_num;
			rhs(row_ID) = P_node * V_down_reg(bus_iter , term_iter - 1);
			for(int term_iter_2 = 1; term_iter_2 < term_iter; ++ term_iter_2){
				rhs(row_ID) += root_i * Q_node(bus_iter , term_iter_2) * V_down_reg(bus_iter, term_iter - term_iter_2);
			}
		}

		// Rhs for PQ Bus
		for(int bus_iter = 0; bus_iter < PQ_bus_num; ++ bus_iter){
			int row_ID = PU_bus_num + bus_iter;
			int bus_ID = Power_market_inform.TSO_Market.power_flow.PQ_bus[bus_iter];

			std::complex <double> S_node = Power_market_inform.TSO_Market.power_flow.P_node(tick, bus_ID);
			S_node += -root_i * Power_market_inform.TSO_Market.power_flow.Q_node(tick, bus_ID);
			rhs(row_ID) = S_node * V_down_hat(bus_iter, term_iter - 1);

			row_ID += 2 * node_small_num;
			S_node = Power_market_inform.TSO_Market.power_flow.P_node(tick, bus_ID);
			S_node += root_i * Power_market_inform.TSO_Market.power_flow.Q_node(tick, bus_ID);
			rhs(row_ID) = S_node * V_down_reg(bus_iter, term_iter - 1);
		}

		// Rhs for reciporal relation for V and 1/V
		for(int term_iter_2 = 1; term_iter_2 < term_iter; ++ term_iter_2){
			rhs.segment(node_small_num, node_small_num) += -V_up_reg.col(term_iter_2).cwiseProduct(V_down_reg.col(term_iter - term_iter_2));
			rhs.segment(3 * node_small_num, node_small_num) += -V_up_hat.col(term_iter_2).cwiseProduct(V_down_hat.col(term_iter - term_iter_2));
		}

		// Rhs for Voltage magnitude constraint
		// Assume all reference voltage level are 1.
		// Rhs for PU Bus
		if(PU_bus_num != 0){
			for(int term_iter_2 = 1; term_iter_2 < term_iter; ++ term_iter_2){
				rhs.tail(PU_bus_num) += -V_up_reg.col(term_iter_2).head(PU_bus_num).cwiseProduct(V_up_hat.col(term_iter - term_iter_2).head(PU_bus_num));
			}
		}
		Eigen::VectorXcd result_temp = Power_market_inform.TSO_Market.power_flow.solver.solve(rhs);

		// -------------------------------------------------------------------------------
		// Update power series
		// -------------------------------------------------------------------------------
		V_up_reg.col(term_iter) = result_temp.head(node_small_num);
		V_up_hat.col(term_iter) = result_temp.segment(node_small_num, node_small_num);
		V_down_reg.col(term_iter) = result_temp.segment(2 * node_small_num, node_small_num);
		V_down_hat.col(term_iter) = result_temp.segment(3 * node_small_num, node_small_num);
		Q_node.col(term_iter) = result_temp.segment(4 * node_small_num, PU_bus_num);
	}
	// Sanity check
	Eigen::VectorXcd V_reg_dir = V_up_reg * Eigen::VectorXcd::Ones(power_terms);
	Eigen::VectorXcd Q_node_dir = Q_node * Eigen::VectorXcd::Ones(power_terms);
	std::cout << V_reg_dir.array().abs().minCoeff() << "\t" << V_reg_dir.array().abs().maxCoeff() << "\n\n";
	std::cout << Q_node_dir.array().real().minCoeff() << "\t" << Q_node_dir.array().real().maxCoeff() << "\n\n";

//	// -------------------------------------------------------------------------------
//	// Pade-approximant
//	// -------------------------------------------------------------------------------
//	// Keep this code commented until spurious poles are fixed!!
//	Eigen::MatrixXcd V_reg_pade = Eigen::MatrixXcd::Zero(node_small_num, power_terms + 1);
//	for(int node_iter = 0; node_iter < node_small_num; ++ node_iter){
//		V_reg_pade.row(node_iter) = local::Pade(V_up_reg.row(node_iter));
//	}
//
//	// Resulting voltage
//	Eigen::VectorXcd V_result_num = Eigen::VectorXcd::Zero(node_small_num);
//	Eigen::VectorXcd V_result_den = Eigen::VectorXcd::Zero(node_small_num);
//	for(int terms_iter = 0; terms_iter < power_terms / 2; ++ terms_iter){
//		V_result_num += V_reg_pade.col(terms_iter);
//		V_result_den += V_reg_pade.col(terms_iter + power_terms / 2);
//	}
//	V_result_den += V_reg_pade.col(power_terms);
//	Eigen::VectorXcd V_result  = V_result_num.array() / V_result_den.array();
//
//	// -------------------------------------------------------------------------------
//	// Roots of the Pade Approximant
//	// -------------------------------------------------------------------------------
//	// Keep this code commented until spurious poles are fixed!!
//	std::vector<Eigen::TripletXcd> Roots_Pade_trip;
//	Eigen::MatrixXcd Roots_Pade_Nodes(node_small_num, power_terms / 2);
//	Eigen::SparseMatrix <std::complex <double>> Roots_Pade(power_terms / 2, power_terms / 2);
//	Eigen::ComplexSchur<Eigen::MatrixXcd> Roots_Pade_Eigen_Solver;
//
//	for(int node_iter = 0; node_iter < node_small_num; ++ node_iter){
//		Roots_Pade_trip.clear();
//		Roots_Pade_trip.reserve(power_terms - 1);
//		for(int terms_iter = 0; terms_iter < power_terms / 2 - 1; ++ terms_iter){
//			Roots_Pade_trip.push_back(Eigen::TripletXcd(terms_iter, terms_iter + 1, 1));
//			Roots_Pade_trip.push_back(Eigen::TripletXcd(power_terms / 2 - 1, terms_iter, V_reg_pade(node_iter, terms_iter + power_terms / 2) / V_reg_pade(node_iter, power_terms)));
//		}
//		Roots_Pade_trip.push_back(Eigen::TripletXcd(power_terms / 2 - 1, power_terms / 2 - 1, V_reg_pade(node_iter, power_terms - 1) / V_reg_pade(node_iter, power_terms)));
//		Roots_Pade.setFromTriplets(Roots_Pade_trip.begin(), Roots_Pade_trip.end());
//		Roots_Pade_Eigen_Solver.compute(Roots_Pade);
//		Roots_Pade_Nodes.row(node_iter) = Roots_Pade_Eigen_Solver.matrixT().diagonal();
//	}
//	// Sanity check
//	std::cout << Roots_Pade_Nodes.array().abs().minCoeff() << "\t" << Roots_Pade_Nodes.array().abs().maxCoeff() << "\n\n";

	// -------------------------------------------------------------------------------
	// Store the results
	// -------------------------------------------------------------------------------
	// P-U Buses
	for(int node_iter = 0; node_iter < PU_bus_num; ++ node_iter){
		int node_ID = Power_market_inform.TSO_Market.power_flow.PU_bus[node_iter];
		Power_market_inform.TSO_Market.power_flow.voltage_abs(tick, node_ID) = abs(V_reg_dir(node_iter));
		Power_market_inform.TSO_Market.power_flow.voltage_arg(tick, node_ID) = arg(V_reg_dir(node_iter));
		Power_market_inform.TSO_Market.power_flow.Q_node(tick, node_ID) = Q_node_dir(node_iter).real();
	}

	// P-Q Buses
	for(int node_iter = 0; node_iter < PQ_bus_num; ++ node_iter){
		int node_ID = Power_market_inform.TSO_Market.power_flow.PQ_bus[node_iter];
		Power_market_inform.TSO_Market.power_flow.voltage_abs(tick, node_ID) = abs(V_reg_dir(PU_bus_num + node_iter));
		Power_market_inform.TSO_Market.power_flow.voltage_arg(tick, node_ID) = arg(V_reg_dir(PU_bus_num + node_iter));
	}

	// Reference Buses
	for(int node_iter = 0; node_iter < ref_bus_num; ++ node_iter){
		int node_ID = Power_market_inform.TSO_Market.power_flow.ref_bus[node_iter];
		Power_market_inform.TSO_Market.power_flow.voltage_abs(tick, node_ID) = 1.;
	}

	// Current on edges
	int edge_num = Power_network_inform.edges.distance.size();
	for(int edge_iter = 0; edge_iter < edge_num; ++ edge_iter){
		int from_ID = Power_network_inform.edges.from(edge_iter);
		int to_ID = Power_network_inform.edges.to(edge_iter);

		double V_abs_from = Power_market_inform.TSO_Market.power_flow.voltage_abs(tick, from_ID);
		double V_arg_from = Power_market_inform.TSO_Market.power_flow.voltage_arg(tick, from_ID);
		std::complex <double> V_from = V_abs_from;
		V_from *= std::complex <double> (cos(V_arg_from), sin(V_arg_from));

		double V_abs_to = Power_market_inform.TSO_Market.power_flow.voltage_abs(tick, to_ID);
		double V_arg_to = Power_market_inform.TSO_Market.power_flow.voltage_arg(tick, to_ID);
		std::complex <double> V_to = V_abs_to;
		V_to *= std::complex <double> (cos(V_arg_to), sin(V_arg_to));

		std::complex <double> current = (V_from - V_to) * Power_market_inform.TSO_Market.power_flow.edge_admittance(edge_iter);
		Power_market_inform.TSO_Market.power_flow.current_abs(tick, edge_iter) = abs(current);
		Power_market_inform.TSO_Market.power_flow.current_arg(tick, edge_iter) = arg(current);
	}
}

void power_network::HELM_Set(network_inform &Power_network_inform, power_market::market_whole_inform &Power_market_inform){
	int node_num = Power_network_inform.nodes.bidding_zone.size();
	int point_num = Power_network_inform.points.bidding_zone.size();
	int bus_num = node_num + point_num;
	int edge_trans_num = Power_network_inform.edges.distance.size();
	int edge_distr_num = 0;
	int DSO_num = Power_network_inform.DSO_cluster.size();
	for(int DSO_iter = 0; DSO_iter < DSO_num; ++ DSO_iter){
		edge_distr_num += Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() * (Power_network_inform.DSO_cluster[DSO_iter].points_ID.size() - 1) / 2;
	}
	double pi = boost::math::constants::pi<double>();

	// -------------------------------------------------------------------------------
	// Set the nodal admittance matrix
	// -------------------------------------------------------------------------------
	std::vector<Eigen::TripletXcd> Y_n_trip;
	Y_n_trip.reserve(bus_num + 2 * (edge_trans_num + edge_distr_num + point_num));
	Eigen::VectorXcd Y_n_Diag = Eigen::VectorXcd::Zero(bus_num);
	Power_network_inform.power_flow.edge_admittance = Eigen::VectorXcd(edge_trans_num);

	// Transmission level
	for(int edge_iter = 0; edge_iter < edge_trans_num; ++ edge_iter){
		int from_ID = Power_network_inform.edges.from(edge_iter);
		int to_ID = Power_network_inform.edges.to(edge_iter);
		int voltage = Power_network_inform.edges.voltage_base(edge_iter);

		// Series admittance
		std::complex <double> y_series(1., 0.);
		y_series /= Power_network_inform.edges.distance(edge_iter);
		y_series /= Power_network_inform.tech_parameters.z_trans_series;
		y_series *= Power_network_inform.tech_parameters.impedenace_base_levels[voltage];
		Power_network_inform.power_flow.edge_admittance(edge_iter) = y_series;

		// Shunt admittance
		std::complex <double> y_shunt(1., 0.);
		y_shunt *= Power_network_inform.edges.distance(edge_iter);
		y_shunt *= Power_network_inform.tech_parameters.y_trans_shunt;
		y_shunt *= Power_network_inform.tech_parameters.impedenace_base_levels[voltage];

		// Triplet for series impedence
		Y_n_trip.push_back(Eigen::TripletXcd(from_ID, to_ID, -y_series));
		Y_n_trip.push_back(Eigen::TripletXcd(to_ID, from_ID, -y_series));

		// Update diagonal terms
		Y_n_Diag(from_ID) += y_series + .5 * y_shunt;
		Y_n_Diag(to_ID) += y_series + .5 * y_shunt;
	}

	// Distribution level
	double z_base_low = Power_network_inform.tech_parameters.z_base_distr();
	double z_base_high = Power_network_inform.tech_parameters.z_base_conn();
	for(int DSO_iter = 0; DSO_iter < DSO_num; ++ DSO_iter){
		int DSO_point_num = Power_network_inform.DSO_cluster[DSO_iter].points_ID.size();
		int DSO_node_num = Power_network_inform.DSO_cluster[DSO_iter].nodes_ID.size();

		double partition_func = 0.;
		Eigen::MatrixXd num_line = Eigen::MatrixXd::Ones(DSO_point_num, DSO_point_num);
		Eigen::MatrixXd distance = Eigen::MatrixXd::Zero(DSO_point_num, DSO_point_num);

		// Connection between points
		for(int row_iter = 0; row_iter < DSO_point_num - 1; ++ row_iter){
			for(int col_iter = row_iter + 1; col_iter < DSO_point_num ; ++ col_iter){
				int point_ID_1 = Power_network_inform.DSO_cluster[DSO_iter].points_ID[row_iter];
				int point_ID_2 = Power_network_inform.DSO_cluster[DSO_iter].points_ID[col_iter];

				distance(row_iter, col_iter) = Power_network_inform.points.distance(point_ID_1, point_ID_2);
				num_line(row_iter, col_iter) = 1. / pow(distance(row_iter, col_iter) * 1E-3, 1. + Power_network_inform.tech_parameters.fraction_dim_distr);
				partition_func += num_line(row_iter, col_iter);
			}
		}
		num_line /= partition_func;
		num_line *= Power_network_inform.tech_parameters.line_density_distr * DSO_point_num;

		for(int row_iter = 0; row_iter < DSO_point_num - 1; ++ row_iter){
			for(int col_iter = row_iter + 1; col_iter < DSO_point_num ; ++ col_iter){
				int point_ID_1 = Power_network_inform.DSO_cluster[DSO_iter].points_ID[row_iter];
				int point_ID_2 = Power_network_inform.DSO_cluster[DSO_iter].points_ID[col_iter];

				// Series admittance
				std::complex <double> y_series(1., 0.);
				y_series /= distance(row_iter, col_iter);
				y_series /= Power_network_inform.tech_parameters.z_distr_series;
				y_series *= z_base_low;
				y_series *= num_line(row_iter, col_iter);

				// Shunt admittance
				std::complex <double> y_shunt(1., 0.);
				y_shunt *= distance(row_iter, col_iter);
				y_shunt *= Power_network_inform.tech_parameters.y_distr_shunt;
				y_shunt *= z_base_low;
				y_shunt *= num_line(row_iter, col_iter);

				// Triplet for series impedence
				Y_n_trip.push_back(Eigen::TripletXcd(node_num + point_ID_1, node_num + point_ID_2, -y_series));
				Y_n_trip.push_back(Eigen::TripletXcd(node_num + point_ID_2, node_num + point_ID_1, -y_series));

				// Update diagonal terms
				Y_n_Diag(node_num + point_ID_1) += y_series + .5 * y_shunt;
				Y_n_Diag(node_num + point_ID_2) += y_series + .5 * y_shunt;
			}
		}

		// Connection between nodes and points
		for(int point_iter = 0; point_iter < DSO_point_num; ++ point_iter){
			int point_ID = Power_network_inform.DSO_cluster[DSO_iter].points_ID[point_iter];
			int node_ID = Power_network_inform.points.node(point_ID);

			Eigen::Vector2d point_coor = Eigen::Vector2d(Power_network_inform.points.lon(point_ID), Power_network_inform.points.lat(point_ID));
			Eigen::Vector2d node_coor = Eigen::Vector2d(Power_network_inform.nodes.lon(node_ID), Power_network_inform.nodes.lat(node_ID));
			point_coor *= pi / 180.;
			node_coor *= pi / 180.;
			double distance_temp = spatial_field::geodist(point_coor, node_coor);

			// Series admittance
			std::complex <double> y_series(1., 0.);
			y_series /= distance_temp;
			y_series /= Power_network_inform.tech_parameters.z_conn_series;
			y_series *= z_base_high;
			y_series *= Power_network_inform.tech_parameters.line_density_conn;

			// Shunt admittance
			std::complex <double> y_shunt(1., 0.);
			y_shunt *= distance_temp;
			y_shunt *= Power_network_inform.tech_parameters.y_conn_shunt;
			y_shunt *= z_base_high;
			y_shunt *= Power_network_inform.tech_parameters.line_density_conn;

			// Triplet for series impedence
			Y_n_trip.push_back(Eigen::TripletXcd(node_ID, node_num + point_ID, -y_series));
			Y_n_trip.push_back(Eigen::TripletXcd(node_num + point_ID, node_ID, -y_series));

			// Update diagonal terms
			Y_n_Diag(node_num + point_ID) += y_series + .5 * y_shunt;
			Y_n_Diag(node_ID) += y_series + .5 * y_shunt;
		}
	}

	// Triplet for diagonal terms
	for(int node_iter = 0; node_iter < bus_num; ++ node_iter){
		Y_n_trip.push_back(Eigen::TripletXcd(node_iter, node_iter, Y_n_Diag(node_iter)));
	}

	// Store the nodal admittance matrix
	Power_network_inform.power_flow.nodal_admittance = Eigen::SparseMatrix <std::complex <double>> (bus_num, bus_num);
	Power_network_inform.power_flow.nodal_admittance.setFromTriplets(Y_n_trip.begin(), Y_n_trip.end());

	// -------------------------------------------------------------------------------
	// Determine type of buses
	// -------------------------------------------------------------------------------
	// Assume only HV power suppliers can give reactive power
	Power_network_inform.power_flow.PQ_bus.reserve(bus_num);
	Power_network_inform.power_flow.PU_bus.reserve(bus_num);
	Power_network_inform.power_flow.ref_bus.reserve(bus_num);
	Eigen::VectorXi bus_type = Eigen::VectorXi::Zero(bus_num);
	int node_ref_num;
	for(int edge_iter = 0; edge_iter < Power_network_inform.cbt.entry_bz.size(); ++ edge_iter){
		if(Power_network_inform.cbt.entry_node_num(edge_iter) == 0){
			continue;
		}
		node_ref_num = (int) Power_network_inform.cbt.entry_nodes(edge_iter, 0);
		break;
	}
	bus_type(node_ref_num) = 2;

	int hydro_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		if(bus_type(node_ID) != 0){
			continue;
		}

		bus_type(node_ID) = 1;
	}

	int wind_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		if(bus_type(node_ID) != 0){
			continue;
		}

		bus_type(node_ID) = 1;
	}

	int pump_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV.size();
	for(int agent_iter = 0; agent_iter < pump_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		if(bus_type(node_ID) != 0){
			continue;
		}

		bus_type(node_ID) = 1;
	}

	for(int bus_iter = 0; bus_iter < bus_num; ++ bus_iter){
		if(bus_type(bus_iter) == 0){
			Power_network_inform.power_flow.PQ_bus.push_back(bus_iter);
		}
		else if(bus_type(bus_iter) == 1){
			Power_network_inform.power_flow.PU_bus.push_back(bus_iter);
		}
		else{
			Power_network_inform.power_flow.ref_bus.push_back(bus_iter);
		}
	}

	// -------------------------------------------------------------------------------
	// Set the solver for iterative process
	// -------------------------------------------------------------------------------
	// Entries from original nodal admittance matrix, reordered and delete reference bus
	std::vector<Eigen::TripletXcd> Permut_Y_n_trip;
	Permut_Y_n_trip.reserve(bus_num);
	for(int bus_iter = 0; bus_iter < Power_network_inform.power_flow.PU_bus.size(); ++ bus_iter){
		int row_ID = bus_iter;
		int bus_ID = Power_network_inform.power_flow.PU_bus[bus_iter];
		Permut_Y_n_trip.push_back(Eigen::TripletXcd(row_ID, bus_ID, 1.));
	}
	for(int bus_iter = 0; bus_iter < Power_network_inform.power_flow.PQ_bus.size(); ++ bus_iter){
		int row_ID = Power_network_inform.power_flow.PU_bus.size() + bus_iter;
		int bus_ID = Power_network_inform.power_flow.PQ_bus[bus_iter];
		Permut_Y_n_trip.push_back(Eigen::TripletXcd(row_ID, bus_ID, 1.));
	}
	Eigen::SparseMatrix <std::complex <double>> Permut_Y_n(Permut_Y_n_trip.size(), bus_num);
	Permut_Y_n.setFromTriplets(Permut_Y_n_trip.begin(), Permut_Y_n_trip.end());
	Eigen::SparseMatrix <std::complex <double>> Y_n_small = Permut_Y_n * Power_network_inform.power_flow.nodal_admittance * Permut_Y_n.transpose();

	// Matrix for solver
	std::vector<Eigen::TripletXcd> Mat_trip;
	Mat_trip.reserve(2 * Y_n_small.nonZeros() + 4 * Power_network_inform.power_flow.PU_bus.size() + 4 * (bus_num));

	// Entries from reduced nodal admittance matrix
	for(int col_iter = 0; col_iter < Y_n_small.outerSize(); ++ col_iter){
		for(Eigen::SparseMatrix<std::complex <double>>::InnerIterator inner_iter(Y_n_small, col_iter); inner_iter; ++ inner_iter){
			std::complex <double> y_conj = inner_iter.value();
			y_conj = std::conj(y_conj);

			Mat_trip.push_back(Eigen::TripletXcd(inner_iter.row(), inner_iter.col(), inner_iter.value()));
			Mat_trip.push_back(Eigen::TripletXcd(2 * Y_n_small.rows() + inner_iter.row(), 2 * Y_n_small.rows() + inner_iter.col(), y_conj));
		}
	}

	// Reciporal relation for V and 1/V
	for(int bus_iter = 0; bus_iter < Y_n_small.rows(); ++ bus_iter){
		int row_ID = Y_n_small.rows() + bus_iter;
		int col_ID_1 = bus_iter;
		int col_ID_2 = Y_n_small.rows() + bus_iter;
		Mat_trip.push_back(Eigen::TripletXcd(row_ID, col_ID_1, 1.));
		Mat_trip.push_back(Eigen::TripletXcd(row_ID, col_ID_2, 1.));

		row_ID += 2 * Y_n_small.rows();
		col_ID_1 += 2 * Y_n_small.rows();
		col_ID_2 += 2 * Y_n_small.rows();
		Mat_trip.push_back(Eigen::TripletXcd(row_ID, col_ID_1, 1.));
		Mat_trip.push_back(Eigen::TripletXcd(row_ID, col_ID_2, 1.));
	}

	// Additional terms for PU buses
	for(int node_iter = 0; node_iter < Power_network_inform.power_flow.PU_bus.size(); ++ node_iter){
		// Entries from reactive power at PU buses
		int col_ID = 4 * Y_n_small.rows() + node_iter;
		std::complex <double> root_i(0., 1.);

		Mat_trip.push_back(Eigen::TripletXcd(node_iter, col_ID, root_i));
		Mat_trip.push_back(Eigen::TripletXcd(node_iter + 2 * Y_n_small.rows(), col_ID, -root_i));

		// Entries for voltage magnitude constraint
		// Assume trivial solution is V = 1. everywhere
		int row_ID = 4 * Y_n_small.rows() + node_iter;
		Mat_trip.push_back(Eigen::TripletXcd(row_ID, node_iter, 1.));
		Mat_trip.push_back(Eigen::TripletXcd(row_ID, node_iter + 2 * Y_n_small.rows(), 1.));
	}
	Eigen::SparseMatrix <std::complex <double>> Mat(4 * Y_n_small.rows() + Power_network_inform.power_flow.PU_bus.size(), 4 * Y_n_small.rows() + Power_network_inform.power_flow.PU_bus.size());
	Mat.setFromTriplets(Mat_trip.begin(), Mat_trip.end());
	Power_network_inform.power_flow.solver.compute(Mat);

	// -------------------------------------------------------------------------------
	// Initialize
	// -------------------------------------------------------------------------------
	int Time = configuration::parameters::Time();

	Power_network_inform.power_flow.P_node = Eigen::MatrixXd::Zero(Time, bus_num);
	Power_network_inform.power_flow.Q_node = Eigen::MatrixXd::Zero(Time, bus_num);
	Power_network_inform.power_flow.voltage_abs = Eigen::MatrixXd::Zero(Time, bus_num);
	Power_network_inform.power_flow.voltage_arg = Eigen::MatrixXd::Zero(Time, bus_num);
	Power_network_inform.power_flow.current_abs = Eigen::MatrixXd::Zero(Time, edge_trans_num);
	Power_network_inform.power_flow.current_arg = Eigen::MatrixXd::Zero(Time, edge_trans_num);

	Power_market_inform.TSO_Market.power_flow.P_node = Eigen::MatrixXd::Zero(Time, node_num);
	Power_market_inform.TSO_Market.power_flow.Q_node = Eigen::MatrixXd::Zero(Time, node_num);
	Power_market_inform.TSO_Market.power_flow.voltage_abs = Eigen::MatrixXd::Zero(Time, node_num);
	Power_market_inform.TSO_Market.power_flow.voltage_arg = Eigen::MatrixXd::Zero(Time, node_num);
	Power_market_inform.TSO_Market.power_flow.current_abs = Eigen::MatrixXd::Zero(Time, edge_trans_num);
	Power_market_inform.TSO_Market.power_flow.current_arg = Eigen::MatrixXd::Zero(Time, edge_trans_num);
}

void power_network::HELM_Node_Update(int tick, network_inform &Power_network_inform, power_market::market_whole_inform &Power_market_inform){
	int node_num = Power_network_inform.nodes.bidding_zone.size();
	int point_num = Power_network_inform.points.bidding_zone.size();

	// Update from cross-border flows
	int edge_num = Power_market_inform.agent_profiles.cross_border.size();
	for(int edge_iter = 0; edge_iter < edge_num; ++ edge_iter){
		int node_num = Power_market_inform.agent_profiles.cross_border[edge_iter].node_num;
		int entry_bz_ID = Power_market_inform.agent_profiles.cross_border[edge_iter].entry_bz_ID;

		if(node_num == 0){
			continue;
		}
		for(int node_iter = 0; node_iter < node_num; ++ node_iter){
			int node_ID = Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].node_ID;
			double real_power;
			real_power = Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results.actual_supply;
			real_power -= Power_market_inform.agent_profiles.cross_border[edge_iter].profiles[node_iter].results.actual_demand;
			Power_network_inform.power_flow.P_node(tick, node_ID) += real_power;
		}
	}

	// Update from end-users
	{
		//int sample_num = agent::end_user::parameters::sample_num();
		int sample_num = Power_market_inform.agent_profiles.end_user_type.cols();
		double apparent_power_base = pow(1. - pow(agent::end_user::parameters::power_factor(), 2.), .5);
		for(int point_iter = 0; point_iter < point_num; ++ point_iter){
			for(int sample_iter = 0; sample_iter < sample_num; ++ sample_iter){
				double real_power;
				real_power = Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_supply;
				Power_network_inform.power_flow.P_node(tick, node_num + point_iter) += real_power;
				real_power = -Power_market_inform.agent_profiles.end_users[point_iter][sample_iter].operation.results.actual_demand;
				Power_network_inform.power_flow.P_node(tick, node_num + point_iter) += real_power;
				Power_network_inform.power_flow.Q_node(tick, node_num + point_iter) += real_power * apparent_power_base;
			}
		}
	}

	// Update from industrial demand
	{
		int industrial_HV_num = Power_market_inform.agent_profiles.industrial.HV.size();
		double apparent_power_base = pow(1. - pow(agent::end_user::parameters::power_factor(), 2.), .5);
		for(int agent_iter = 0; agent_iter < industrial_HV_num; ++ agent_iter){
			int point_ID = Power_market_inform.agent_profiles.industrial.HV[agent_iter].point_ID;
			int node_ID = Power_network_inform.points.node(point_ID);

			double real_power = -Power_market_inform.agent_profiles.industrial.HV[agent_iter].results.actual_demand;
			Power_network_inform.power_flow.P_node(tick, node_ID) += real_power;
			Power_network_inform.power_flow.Q_node(tick, node_ID) += real_power * apparent_power_base;
		}
	}

	// Update from power suppliers
	int hydro_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		double real_power = Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		Power_network_inform.power_flow.P_node(tick, node_ID) += real_power;
	}

	int hydro_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant.size();
	for(int agent_iter = 0; agent_iter < hydro_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].point_ID;

		double real_power = Power_market_inform.agent_profiles.power_supplier.hydro.LV_plant[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		Power_network_inform.power_flow.P_node(tick, node_num + point_ID) += real_power;
	}

	int wind_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		double real_power = Power_market_inform.agent_profiles.power_supplier.wind.HV_plant[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		Power_network_inform.power_flow.P_node(tick, node_ID) += real_power;
	}

	int wind_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant.size();
	for(int agent_iter = 0; agent_iter < wind_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].point_ID;

		double real_power = Power_market_inform.agent_profiles.power_supplier.wind.LV_plant[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		Power_network_inform.power_flow.P_node(tick, node_num + point_ID) += real_power;
	}

	int pump_HV_plant_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV.size();
	for(int agent_iter = 0; agent_iter < pump_HV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		double real_power = Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		real_power -= Power_market_inform.agent_profiles.power_supplier.pump_storage.HV[agent_iter].results.actual_demand;
		Power_network_inform.power_flow.P_node(tick, node_ID) += real_power;
	}

	int pump_LV_plant_num = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV.size();
	for(int agent_iter = 0; agent_iter < pump_LV_plant_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].point_ID;

		double real_power = Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		real_power -= Power_market_inform.agent_profiles.power_supplier.pump_storage.LV[agent_iter].results.actual_demand;
		Power_network_inform.power_flow.P_node(tick, node_num + point_ID) += real_power;
	}

	int slack_HV_num = Power_market_inform.agent_profiles.power_supplier.slack.HV_plant.size();
	for(int agent_iter = 0; agent_iter < slack_HV_num; ++ agent_iter){
		int point_ID = Power_market_inform.agent_profiles.power_supplier.slack.HV_plant[agent_iter].point_ID;
		int node_ID = Power_network_inform.points.node(point_ID);

		double real_power = Power_market_inform.agent_profiles.power_supplier.slack.HV_plant[agent_iter].results.actual_supply;
		real_power /= 1. - power_network::parameters::loss_factor();
		real_power -= Power_market_inform.agent_profiles.power_supplier.slack.HV_plant[agent_iter].results.actual_demand;
		Power_network_inform.power_flow.P_node(tick, node_ID) += real_power;
	}
}

void power_network::HELM_Solve(int tick, network_inform &Power_network_inform, power_market::market_whole_inform& Power_market_inform){
	int node_num = Power_network_inform.nodes.bidding_zone.size();
	int point_num = Power_network_inform.points.bidding_zone.size();
	int PU_bus_num = Power_network_inform.power_flow.PU_bus.size();
	int PQ_bus_num = Power_network_inform.power_flow.PQ_bus.size();
	int ref_bus_num = Power_network_inform.power_flow.ref_bus.size();
	int bus_num = node_num + point_num;
	int bus_small_num = bus_num - ref_bus_num;
	std::complex <double> root_i(0., 1.);

	// -------------------------------------------------------------------------------
	// Initialization of power series coefficients
	// -------------------------------------------------------------------------------
	int power_terms = 20;
	Eigen::MatrixXcd V_up_reg = Eigen::MatrixXcd::Zero(bus_small_num, power_terms);
	Eigen::MatrixXcd V_up_hat = Eigen::MatrixXcd::Zero(bus_small_num, power_terms);
	Eigen::MatrixXcd V_down_reg = Eigen::MatrixXcd::Zero(bus_small_num, power_terms);
	Eigen::MatrixXcd V_down_hat = Eigen::MatrixXcd::Zero(bus_small_num, power_terms);
	Eigen::MatrixXcd Q_node = Eigen::MatrixXcd::Zero(PU_bus_num, power_terms);

	// -------------------------------------------------------------------------------
	// Trivial solution (V = 1 and S = 0 everywhere)
	// -------------------------------------------------------------------------------
	V_up_reg.col(0) = Eigen::VectorXcd::Ones(bus_small_num);
	V_up_hat.col(0) = Eigen::VectorXcd::Ones(bus_small_num);
	V_down_reg.col(0) = Eigen::VectorXcd::Ones(bus_small_num);
	V_down_hat.col(0) = Eigen::VectorXcd::Ones(bus_small_num);

	// -------------------------------------------------------------------------------
	// Iteratively solve the linear equations
	// -------------------------------------------------------------------------------
	for(int term_iter = 1; term_iter < power_terms; ++ term_iter){
		//std::cout << term_iter << "\n\n";
		Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(4 * bus_small_num + PU_bus_num);

		// Rhs for PU Bus
		for(int bus_iter = 0; bus_iter < PU_bus_num; ++ bus_iter){
			int row_ID = bus_iter;
			int bus_ID = Power_network_inform.power_flow.PU_bus[bus_iter];

			double P_node = Power_network_inform.power_flow.P_node(tick, bus_ID);
			rhs(row_ID) = P_node * V_down_hat(bus_iter , term_iter - 1);
			for(int term_iter_2 = 1; term_iter_2 < term_iter; ++ term_iter_2){
				rhs(row_ID) += -root_i * Q_node(bus_iter , term_iter_2) * V_down_hat(bus_iter , term_iter - term_iter_2);
			}

			row_ID += 2 * bus_small_num;
			rhs(row_ID) = P_node * V_down_reg(bus_iter , term_iter - 1);
			for(int term_iter_2 = 1; term_iter_2 < term_iter; ++ term_iter_2){
				rhs(row_ID) += root_i * Q_node(bus_iter , term_iter_2) * V_down_reg(bus_iter, term_iter - term_iter_2);
			}
		}

		// Rhs for PQ Bus
		for(int bus_iter = 0; bus_iter < PQ_bus_num; ++ bus_iter){
			int row_ID = PU_bus_num + bus_iter;
			int bus_ID = Power_network_inform.power_flow.PQ_bus[bus_iter];

			std::complex <double> S_node = Power_network_inform.power_flow.P_node(tick, bus_ID);
			S_node += -root_i * Power_network_inform.power_flow.Q_node(tick, bus_ID);
			rhs(row_ID) = S_node * V_down_hat(bus_iter, term_iter - 1);

			row_ID += 2 * bus_small_num;
			S_node = Power_network_inform.power_flow.P_node(tick, bus_ID);
			S_node += root_i * Power_network_inform.power_flow.Q_node(tick, bus_ID);
			rhs(row_ID) = S_node * V_down_reg(bus_iter, term_iter - 1);
		}

		// Rhs for reciporal relation for V and 1/V
		for(int term_iter_2 = 1; term_iter_2 < term_iter; ++ term_iter_2){
			rhs.segment(bus_small_num, bus_small_num) += -V_up_reg.col(term_iter_2).cwiseProduct(V_down_reg.col(term_iter - term_iter_2));
			rhs.segment(3 * bus_small_num, bus_small_num) += -V_up_hat.col(term_iter_2).cwiseProduct(V_down_hat.col(term_iter - term_iter_2));
		}

		// Rhs for Voltage magnitude constraint
		// Assume all reference voltage level are 1.
		// Rhs for PU Bus
		if(PU_bus_num != 0){
			for(int term_iter_2 = 1; term_iter_2 < term_iter; ++ term_iter_2){
				rhs.tail(PU_bus_num) += -V_up_reg.col(term_iter_2).head(PU_bus_num).cwiseProduct(V_up_hat.col(term_iter - term_iter_2).head(PU_bus_num));
			}
		}
		Eigen::VectorXcd result_temp = Power_network_inform.power_flow.solver.solve(rhs);

		// -------------------------------------------------------------------------------
		// Update power series
		// -------------------------------------------------------------------------------
		V_up_reg.col(term_iter) = result_temp.head(bus_small_num);
		V_up_hat.col(term_iter) = result_temp.segment(bus_small_num, bus_small_num);
		V_down_reg.col(term_iter) = result_temp.segment(2 * bus_small_num, bus_small_num);
		V_down_hat.col(term_iter) = result_temp.segment(3 * bus_small_num, bus_small_num);
		Q_node.col(term_iter) = result_temp.segment(4 * bus_small_num, PU_bus_num);
	}
	// Sanity check
	Eigen::VectorXcd V_reg_dir = V_up_reg * Eigen::VectorXcd::Ones(power_terms);
	Eigen::VectorXcd Q_node_dir = Q_node * Eigen::VectorXcd::Ones(power_terms);
	std::cout << V_reg_dir.array().abs().minCoeff() << "\t" << V_reg_dir.array().abs().maxCoeff() << "\n\n";
	std::cout << Q_node_dir.array().real().minCoeff() << "\t" << Q_node_dir.array().real().maxCoeff() << "\n\n";

//	// -------------------------------------------------------------------------------
//	// Pade-approximant
//	// -------------------------------------------------------------------------------
//	// Keep this code commented until spurious poles are fixed!!
//	Eigen::MatrixXcd V_reg_pade = Eigen::MatrixXcd::Zero(bus_small_num, power_terms + 1);
//	for(int bus_iter = 0; bus_iter < bus_small_num; ++ bus_iter){
//		V_reg_pade.row(bus_iter) = local::Pade(V_up_reg.row(bus_iter));
//	}
//
//	// Resulting voltage
//	Eigen::VectorXcd V_result_num = Eigen::VectorXcd::Zero(bus_small_num);
//	Eigen::VectorXcd V_result_den = Eigen::VectorXcd::Zero(bus_small_num);
//	for(int terms_iter = 0; terms_iter < power_terms / 2; ++ terms_iter){
//		V_result_num += V_reg_pade.col(terms_iter);
//		V_result_den += V_reg_pade.col(terms_iter + power_terms / 2);
//	}
//	V_result_den += V_reg_pade.col(power_terms);
//	Eigen::VectorXcd V_result  = V_result_num.array() / V_result_den.array();
//
//	// -------------------------------------------------------------------------------
//	// Roots of the Pade Approximant
//	// -------------------------------------------------------------------------------
//	// Keep this code commented until spurious poles are fixed!!
//	std::vector<Eigen::TripletXcd> Roots_Pade_trip;
//	Eigen::MatrixXcd Roots_Pade_buss(bus_small_num, power_terms / 2);
//	Eigen::SparseMatrix <std::complex <double>> Roots_Pade(power_terms / 2, power_terms / 2);
//	Eigen::ComplexSchur<Eigen::MatrixXcd> Roots_Pade_Eigen_Solver;
//
//	for(int bus_iter = 0; bus_iter < bus_small_num; ++ bus_iter){
//		Roots_Pade_trip.clear();
//		Roots_Pade_trip.reserve(power_terms - 1);
//		for(int terms_iter = 0; terms_iter < power_terms / 2 - 1; ++ terms_iter){
//			Roots_Pade_trip.push_back(Eigen::TripletXcd(terms_iter, terms_iter + 1, 1));
//			Roots_Pade_trip.push_back(Eigen::TripletXcd(power_terms / 2 - 1, terms_iter, V_reg_pade(bus_iter, terms_iter + power_terms / 2) / V_reg_pade(bus_iter, power_terms)));
//		}
//		Roots_Pade_trip.push_back(Eigen::TripletXcd(power_terms / 2 - 1, power_terms / 2 - 1, V_reg_pade(bus_iter, power_terms - 1) / V_reg_pade(bus_iter, power_terms)));
//		Roots_Pade.setFromTriplets(Roots_Pade_trip.begin(), Roots_Pade_trip.end());
//		Roots_Pade_Eigen_Solver.compute(Roots_Pade);
//		Roots_Pade_buss.row(bus_iter) = Roots_Pade_Eigen_Solver.matrixT().diagonal();
//	}
//	// Sanity check
//	std::cout << Roots_Pade_buss.array().abs().minCoeff() << "\t" << Roots_Pade_buss.array().abs().maxCoeff() << "\n\n";

	// -------------------------------------------------------------------------------
	// Store the results, whole
	// -------------------------------------------------------------------------------
	// P-U Buses
	for(int node_iter = 0; node_iter < PU_bus_num; ++ node_iter){
		int node_ID = Power_network_inform.power_flow.PU_bus[node_iter];

		Power_network_inform.power_flow.voltage_abs(tick, node_ID) = abs(V_reg_dir(node_iter));
		Power_network_inform.power_flow.voltage_arg(tick, node_ID) = arg(V_reg_dir(node_iter));
		Power_network_inform.power_flow.Q_node(tick, node_ID) = Q_node_dir(node_iter).real();

		if(node_ID < node_num){
			Power_market_inform.TSO_Market.power_flow.voltage_abs(tick, node_ID) = Power_network_inform.power_flow.voltage_abs(tick, node_ID);
			Power_market_inform.TSO_Market.power_flow.voltage_arg(tick, node_ID) = Power_network_inform.power_flow.voltage_arg(tick, node_ID);
			Power_market_inform.TSO_Market.power_flow.P_node(tick, node_ID) += Power_network_inform.power_flow.P_node(tick, node_ID);
			Power_market_inform.TSO_Market.power_flow.Q_node(tick, node_ID) += Power_network_inform.power_flow.Q_node(tick, node_ID);
		}
		else{
			int point_ID = node_ID - node_num;
			int trans_node_ID = Power_network_inform.points.node(trans_node_ID);
			Power_market_inform.TSO_Market.power_flow.P_node(tick, trans_node_ID) += Power_network_inform.power_flow.P_node(tick, node_ID);
			Power_market_inform.TSO_Market.power_flow.Q_node(tick, trans_node_ID) += Power_network_inform.power_flow.Q_node(tick, node_ID);
		}
	}

	// P-Q Buses
	for(int node_iter = 0; node_iter < PQ_bus_num; ++ node_iter){
		int node_ID = Power_network_inform.power_flow.PQ_bus[node_iter];
		Power_network_inform.power_flow.voltage_abs(tick, node_ID) = abs(V_reg_dir(PU_bus_num + node_iter));
		Power_network_inform.power_flow.voltage_arg(tick, node_ID) = arg(V_reg_dir(PU_bus_num + node_iter));

		if(node_ID < node_num){
			Power_market_inform.TSO_Market.power_flow.voltage_abs(tick, node_ID) = Power_network_inform.power_flow.voltage_abs(tick, node_ID);
			Power_market_inform.TSO_Market.power_flow.voltage_arg(tick, node_ID) = Power_network_inform.power_flow.voltage_arg(tick, node_ID);
			Power_market_inform.TSO_Market.power_flow.P_node(tick, node_ID) += Power_network_inform.power_flow.P_node(tick, node_ID);
			Power_market_inform.TSO_Market.power_flow.Q_node(tick, node_ID) += Power_network_inform.power_flow.Q_node(tick, node_ID);
		}
		else{
			int point_ID = node_ID - node_num;
			int trans_node_ID = Power_network_inform.points.node(point_ID);
			Power_market_inform.TSO_Market.power_flow.P_node(tick, trans_node_ID) += Power_network_inform.power_flow.P_node(tick, node_ID);
			Power_market_inform.TSO_Market.power_flow.Q_node(tick, trans_node_ID) += Power_network_inform.power_flow.Q_node(tick, node_ID);
		}
	}

	// Reference Buses
	for(int node_iter = 0; node_iter < ref_bus_num; ++ node_iter){
		int node_ID = Power_network_inform.power_flow.ref_bus[node_iter];
		Power_network_inform.power_flow.voltage_abs(tick, node_ID) = 1.;

		if(node_ID < node_num){
			Power_market_inform.TSO_Market.power_flow.voltage_abs(tick, node_ID) = 1.;
		}
	}

	// Current on edges
	int edge_trans_num = Power_network_inform.edges.distance.size();
	for(int edge_iter = 0; edge_iter < edge_trans_num; ++ edge_iter){
		int from_ID = Power_network_inform.edges.from(edge_iter);
		int to_ID = Power_network_inform.edges.to(edge_iter);

		double V_abs_from = Power_network_inform.power_flow.voltage_abs(tick, from_ID);
		double V_arg_from = Power_network_inform.power_flow.voltage_arg(tick, from_ID);
		std::complex <double> V_from = V_abs_from;
		V_from *= std::complex <double> (cos(V_arg_from), sin(V_arg_from));

		double V_abs_to = Power_network_inform.power_flow.voltage_abs(tick, to_ID);
		double V_arg_to = Power_network_inform.power_flow.voltage_arg(tick, to_ID);
		std::complex <double> V_to = V_abs_to;
		V_to *= std::complex <double> (cos(V_arg_to), sin(V_arg_to));

		std::complex <double> current = (V_from - V_to) * Power_network_inform.power_flow.edge_admittance(edge_iter);
		Power_network_inform.power_flow.current_abs(tick, edge_iter) = abs(current);
		Power_network_inform.power_flow.current_arg(tick, edge_iter) = arg(current);
		Power_market_inform.TSO_Market.power_flow.current_abs(tick, edge_iter) = abs(current);
		Power_market_inform.TSO_Market.power_flow.current_arg(tick, edge_iter) = arg(current);
	}
}
