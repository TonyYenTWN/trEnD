#include "contingency_analysis.h"

namespace local{
	void Gibbs_sampling(int duration, double beta, Eigen::MatrixXi &sample, Eigen::MatrixXd &temporal_hamiltonian, Eigen::MatrixXd &temporal_prob,  Eigen::MatrixXd &spatial_hamiltonian){
        int max_int = 1E12;
	    std::random_device rd;  // a seed source for the random number engine
	    std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<> distrib(0, max_int);

//     // Update when spatial correlation between components is added
//        if(spatial_hamiltonian.size() != 0){
//        	//
//		}

		// Metropolis walk for each components
		for(int row_iter = 0; row_iter < sample.rows(); ++ row_iter){

			// Gibbs sampling for the first time step
			{
				int tick = 0;
				// Randomly set the boundary condition from theoretical stationary probability distribution
				int s_0 = (distrib(gen) < temporal_prob(row_iter, 0) / temporal_prob(row_iter, 1) * max_int);

				// Calculate energy gap and the resulting conditional likelihood of s = 1 happening
				int s_neighbor_sum = s_0 + sample(row_iter, tick + 1);
				double energy_gap = temporal_hamiltonian(row_iter, 0) - temporal_hamiltonian(row_iter, 1) * s_neighbor_sum;
				double ratio = exp(-beta * energy_gap);

				// Update state according to the conditional probability
				sample(row_iter, tick) = (distrib(gen) <= ratio / (1. + ratio) * max_int);
			}

			// Gibbs sampling for the steps
			for(int tick = 1; tick < duration - 1; ++ tick){
				// Calculate energy gap and the resulting conditional likelihood of s = 1 happening
				int s_neighbor_sum = sample(row_iter, tick - 1) + sample(row_iter, tick + 1);
				double energy_gap = temporal_hamiltonian(row_iter, 0) - temporal_hamiltonian(row_iter, 1) * s_neighbor_sum;
				double ratio = exp(-beta * energy_gap);

				// Update state according to the conditional probability
				sample(row_iter, tick) = (distrib(gen) <= ratio / (1. + ratio) * max_int);
			}

			// Gibbs sampling for the final step
			{
				int tick = duration - 1;

				// Randomly set the boundary condition from theoretical stationary probability distribution
				int s_T = (distrib(gen) < temporal_prob(row_iter, 0) / temporal_prob(row_iter, 1) * max_int);

				// Calculate energy gap and the resulting conditional likelihood of s = 1 happening
				int s_neighbor_sum = sample(row_iter, tick - 1) + s_T;
				double energy_gap = temporal_hamiltonian(row_iter, 0) - temporal_hamiltonian(row_iter, 1) * s_neighbor_sum;
				double ratio = exp(-beta * energy_gap);

				// Update state according to the conditional probability
				sample(row_iter, tick) = (distrib(gen) <= ratio / (1. + ratio) * max_int);
			}
		}
	}

	void contingency_analysis_LP_set(int sample_ID, int tick, power_network::contingency_analysis_struct &contingency_analysis, power_market::market_inform &Market){
        // -------------------------------------------------------------------------------
        // LP Solver initialization for flow-based market optimization
        // Need to set up each time a contingency occurs (constraints are different)
        // Variables are sorted as {V}, {S} (power injection), {S'} (supply), {D}_{inf}, {D}_{flex}, {BESS}, {EV}, {D}_{shift, t - tau : t + tau}, {I} (current)
        // for {BESS} and {EV}, positive  = discharge
        // for consistency, {D} are always negative
        // -------------------------------------------------------------------------------

        // -------------------------------------------------------------------------------
        // Set matrix for general constraints
        // -------------------------------------------------------------------------------
        // Construct node admittance matrix
        std::vector <Eigen::TripletXd> Y_n_trip;
        Y_n_trip.reserve(2 * Market.network.num_edges + Market.network.num_vertice);
        Eigen::SparseMatrix <double, Eigen::RowMajor> Y_n(Market.network.num_vertice, Market.network.num_vertice);
        Eigen::VectorXd Y_n_diag = Eigen::VectorXd::Zero(Market.network.num_vertice);
        Eigen::VectorXpd Connection_num = Eigen::VectorXpd::Ones(Market.network.num_vertice);
        for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
            int from_ID = Market.network.incidence[edge_iter](0);
            int to_ID = Market.network.incidence[edge_iter](1);

            double y_edge = Market.network.admittance[edge_iter];
            // consider the contingency: failure of transformer on either end or the line itself
            y_edge *= contingency_analysis.samples[sample_ID](Market.network.num_vertice + edge_iter, tick);
            y_edge *= contingency_analysis.samples[sample_ID](from_ID, tick);
            y_edge *= contingency_analysis.samples[sample_ID](to_ID, tick);

            // Equality constraints of voltage - source / sink at the nodes, off-diagonal terms
            Y_n_trip.push_back(Eigen::TripletXd(from_ID, to_ID, -y_edge));
            Y_n_trip.push_back(Eigen::TripletXd(to_ID, from_ID, -y_edge));
            Connection_num(from_ID) += 1;
            Connection_num(to_ID) += 1;

            // Equality constraints of voltage - source / sink at the nodes, diagonal terms
            Y_n_diag(from_ID) += y_edge;
            Y_n_diag(to_ID) += y_edge;
        }
        for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
            // Equality constraints of voltage - source / sink at the nodes, summed diagonal terms
            Y_n_trip.push_back(Eigen::TripletXd(node_iter, node_iter, Y_n_diag(node_iter)));
        }
        Y_n.setFromTriplets(Y_n_trip.begin(), Y_n_trip.end());

        // Generate sparse matrix for general (equality) constraints of voltage, power flow, and source / sink summation
        int constrant_num = 2 * Market.network.num_vertice + Market.network.num_edges;
        int variable_num = Market.network.num_vertice * (Market.flex_stat.unfulfilled_demand.rows() + 6) + Market.network.num_edges;
        Eigen::VectorXpd non_zero_num(constrant_num);
        non_zero_num << Connection_num + Eigen::VectorXpd::Ones(Market.network.num_vertice), Eigen::VectorXpd::Constant(Market.network.num_edges, 3), Eigen::VectorXpd::Constant(Market.network.num_vertice, Market.flex_stat.unfulfilled_demand.rows() + 5);
        alglib::integer_1d_array row_sizes_general;
        row_sizes_general.setcontent(non_zero_num.size(), non_zero_num.data());
        alglib::sparsematrix constraint_general;
        alglib::sparsecreatecrs(constrant_num, variable_num, row_sizes_general, constraint_general);
        // Rows for voltage - source / sink equalities
        for(int row_iter = 0; row_iter < Y_n.outerSize(); ++ row_iter){
            for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator inner_iter(Y_n, row_iter); inner_iter; ++ inner_iter){
                alglib::sparseset(constraint_general, inner_iter.row(), inner_iter.col(), inner_iter.value());
            }
            // Update the columns right to the node admittance block
            alglib::sparseset(constraint_general, row_iter, Market.network.num_vertice + row_iter, -1.);
        }
        // Rows for voltage - power flow equalities
        for(int edge_iter = 0; edge_iter < Market.network.num_edges; ++ edge_iter){
            int from_ID = Market.network.incidence[edge_iter](0);
            int to_ID = Market.network.incidence[edge_iter](1);

            double y_edge = Market.network.admittance[edge_iter];
            // consider the contingency: failure of transformer on either end or the line itself
            y_edge *= contingency_analysis.samples[sample_ID](Market.network.num_vertice + edge_iter, tick);
            y_edge *= contingency_analysis.samples[sample_ID](from_ID, tick);
            y_edge *= contingency_analysis.samples[sample_ID](to_ID, tick);

            if( from_ID < to_ID){
                alglib::sparseset(constraint_general, Y_n.rows() + edge_iter,  from_ID, y_edge);
                alglib::sparseset(constraint_general, Y_n.rows() + edge_iter, to_ID, -y_edge);
            }
            else{
                alglib::sparseset(constraint_general, Y_n.rows() + edge_iter,  from_ID, -y_edge);
                alglib::sparseset(constraint_general, Y_n.rows() + edge_iter, to_ID, y_edge);
            }
            alglib::sparseset(constraint_general, Y_n.rows() + edge_iter, Market.network.num_vertice * (Market.flex_stat.unfulfilled_demand.rows() + 6) + edge_iter, -1.);
        }
        // Rows for source / sink summation at each node
        // {S} - {S'} (supply) - {D}_{inf} - {D}_{flex} - {D}_{shift, t - tau : t + tau} - {BESS} - {EV} = 0
        for(int node_iter = 0; node_iter < Y_n.rows(); ++ node_iter){
            int row_ID = Y_n.rows() + Market.network.num_edges + node_iter;
            alglib::sparseset(constraint_general, row_ID, Y_n.rows() + node_iter, 1.);
            int col_start = 2 * Y_n.rows() + node_iter * (Market.flex_stat.unfulfilled_demand.rows() + 4);
            int col_end = 2 * Y_n.rows() + (node_iter + 1) * (Market.flex_stat.unfulfilled_demand.rows() + 4) - 1;
            for(int col_iter = col_start; col_iter <= col_end; ++ col_iter){
                alglib::sparseset(constraint_general, row_ID, col_iter, -1.);
            }
        }

        // -------------------------------------------------------------------------------
        // Set bounds for general and box constraints
        // -------------------------------------------------------------------------------
        Eigen::MatrixXd bound_general = Eigen::MatrixXd::Zero(constrant_num, 2);
        Eigen::MatrixXd bound_box = Eigen::MatrixXd::Zero(variable_num, 2);
        bound_box.topRows(Market.network.num_vertice) = Market.network.voltage_constraint;
        bound_box.middleRows(Market.network.num_vertice, Market.network.num_vertice).col(0) = Eigen::VectorXd::Constant(Market.network.num_vertice, -std::numeric_limits<double>::infinity());
        bound_box.middleRows(Market.network.num_vertice, Market.network.num_vertice).col(1) = Eigen::VectorXd::Constant(Market.network.num_vertice, std::numeric_limits<double>::infinity());
        for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
            int row_start = 2 * Market.network.num_vertice + node_iter * (Market.flex_stat.unfulfilled_demand.rows() + 4);

            // default constraint on supply (before considering failure power plants)
            bound_box.row(row_start) << 0., Market.flex_stat.supply_inflex(tick, node_iter) + Market.flex_stat.supply_flex(tick, node_iter);

            // constraint on inflexible demand
            // shiftable demand that is deferred to the last possible time interval is considered inflexible demand
            bound_box.row(row_start + 1) << -Market.flex_stat.demand_inflex(tick, node_iter) - Market.flex_stat.unfulfilled_demand(0, node_iter), 0.;

            // constraint on flexible demand
            bound_box.row(row_start + 2) << -Market.flex_stat.demand_flex(tick, node_iter), 0.;

            // constraint on power from BESS
            double BESS_ub = Market.flex_stat.BESS_soc.soc_current(node_iter) - Market.flex_stat.BESS_soc.soc_min(tick, node_iter);
            BESS_ub = std::min(BESS_ub, Market.flex_stat.BESS_soc.capacity_max(node_iter));
            double BESS_lb = Market.flex_stat.BESS_soc.soc_current(node_iter) - Market.flex_stat.BESS_soc.soc_max(tick, node_iter);
            BESS_lb = std::max(BESS_lb, -Market.flex_stat.BESS_soc.capacity_max(node_iter));
            bound_box.row(row_start + 3) << BESS_lb, BESS_ub;

            // constraint on power from EV
            double EV_ub = Market.flex_stat.EV_soc.soc_current(node_iter) - Market.flex_stat.EV_soc.soc_min(tick, node_iter);
            EV_ub = std::min(EV_ub, Market.flex_stat.EV_soc.capacity_max(node_iter));
            double EV_lb = Market.flex_stat.EV_soc.soc_current(node_iter) - Market.flex_stat.EV_soc.soc_max(tick, node_iter);
            EV_lb = std::max(EV_lb, -Market.flex_stat.EV_soc.capacity_max(node_iter));
            bound_box.row(row_start + 4) << EV_lb, EV_ub;

            // constraint on shiftable demand
            for(int tock = 1; tock < Market.flex_stat.unfulfilled_demand.rows(); ++ tock){
                bound_box.row(row_start + 4 + tock) << -Market.flex_stat.unfulfilled_demand(tock, node_iter), 0.;
            }
        }
        bound_box.bottomRows(Market.network.num_edges) = Market.network.power_constraint;

        // Bounds of general constraints
        alglib::real_1d_array lb_general;
        alglib::real_1d_array ub_general;
        lb_general.setcontent(bound_general.rows(), bound_general.col(0).data());
        ub_general.setcontent(bound_general.rows(), bound_general.col(1).data());

        // Bounds of box constraints
        alglib::real_1d_array lb_box;
        alglib::real_1d_array ub_box;
        lb_box.setcontent(bound_box.rows(), bound_box.col(0).data());
        ub_box.setcontent(bound_box.rows(), bound_box.col(1).data());

        // -------------------------------------------------------------------------------
        // Set scale of variables
        // -------------------------------------------------------------------------------
        Eigen::VectorXd scale_vec = Eigen::VectorXd::Ones(variable_num);
        scale_vec.head(Market.network.num_vertice) = bound_box.col(1).head(Market.network.num_vertice) - bound_box.col(0).head(Market.network.num_vertice);
        scale_vec.tail(Market.network.num_edges) = bound_box.col(1).tail(Market.network.num_edges) - bound_box.col(0).tail(Market.network.num_edges);
        alglib::real_1d_array scale;
        scale.setcontent(scale_vec.size(), scale_vec.data());

        // -------------------------------------------------------------------------------
        // Set objective coefficients of variables
        // -------------------------------------------------------------------------------
        Eigen::VectorXd obj_vec = Eigen::VectorXd::Zero(variable_num);
        for(int node_iter = 0; node_iter < Market.network.num_vertice; ++ node_iter){
            int row_start = 2 * Market.network.num_vertice + node_iter * (Market.flex_stat.unfulfilled_demand.rows() + 4);

            // Only inflexible demand is considered in the LP problem; should fulfill as much as possible
            obj_vec(row_start + 1) = 1.;
        }
        alglib::real_1d_array obj_coeff;
        obj_coeff.setcontent(obj_vec.size(), obj_vec.data());

        // -------------------------------------------------------------------------------
        // Set the LP problem object
        // -------------------------------------------------------------------------------
        alglib::minlpcreate(variable_num, contingency_analysis.problem);
        alglib::minlpsetcost(contingency_analysis.problem, obj_coeff);
        alglib::minlpsetlc2(contingency_analysis.problem, constraint_general, lb_general, ub_general, constrant_num);
        alglib::minlpsetbc(contingency_analysis.problem, lb_box, ub_box);
        alglib::minlpsetscale(contingency_analysis.problem, scale);
        alglib::minlpsetalgodss(contingency_analysis.problem, 0.);

        // -------------------------------------------------------------------------------
        // Solve the problem
        // -------------------------------------------------------------------------------
        alglib::minlpoptimize(contingency_analysis.problem);
	}

	void contingency_analysis_update(int sample_ID, int tick, power_network::contingency_analysis_struct &contingency_analysis, power_market::market_inform &Market){

	}
}

namespace power_network{
    void flex_stat_input(power_market::market_whole_inform &Power_market_inform, power_network::network_inform &Power_network_inform, configuration::process_config &process_par){

        if(!process_par.simulation_flag){
             int Time = process_par.total_time;

            // Initialization of the TSO
            power_market::TSO_Market_Set(Power_market_inform.TSO_Market, Power_network_inform, Time);

            // Read flex stat data
            std::string dir_name = "csv/case/" + process_par.folder_name + "/processed/power_market/flex_stat/";

            auto fin = dir_name + "demand_flex.csv";
            auto fin_dim = basic::get_file_dim(fin);
            Power_market_inform.TSO_Market.flex_stat.demand_flex = basic::read_file(fin_dim[0], fin_dim[1], fin);

            fin = dir_name + "demand_inflex.csv";
            fin_dim = basic::get_file_dim(fin);
            Power_market_inform.TSO_Market.flex_stat.demand_flex = basic::read_file(fin_dim[0], fin_dim[1], fin);

            fin = dir_name + "demand_shiftable.csv";
            fin_dim = basic::get_file_dim(fin);
            Power_market_inform.TSO_Market.flex_stat.demand_flex = basic::read_file(fin_dim[0], fin_dim[1], fin);

            fin = dir_name + "supply_flex.csv";
            fin_dim = basic::get_file_dim(fin);
            Power_market_inform.TSO_Market.flex_stat.demand_flex = basic::read_file(fin_dim[0], fin_dim[1], fin);

            fin = dir_name + "supply_inflex.csv";
            fin_dim = basic::get_file_dim(fin);
            Power_market_inform.TSO_Market.flex_stat.demand_flex = basic::read_file(fin_dim[0], fin_dim[1], fin);

            // Read end-user data
            // Necessary when BESS and EV data are required
            // use the new read config function to read the data
            std::string end_user_type_file = "csv/case/" + process_par.folder_name + "/input/agent/end_user_types.csv";
        }

        // Set shiftable demand data
        int load_shift_time = agent::end_user::parameters::load_shift_time();
        int foresight_time = agent::end_user::parameters::foresight_time();
        load_shift_time = std::min(load_shift_time, foresight_time / 2);
        // setting for trivial case
        Power_market_inform.TSO_Market.flex_stat.unfulfilled_demand = Eigen::MatrixXd::Zero(2 * load_shift_time + 1, Power_market_inform.TSO_Market.num_zone);

        // Set BESS and EV soc range
        // setting for trivial case
        Power_market_inform.TSO_Market.flex_stat.BESS_soc.soc_min = Eigen::MatrixXd::Zero(process_par.time_boundary[1], Power_market_inform.TSO_Market.num_zone);
        Power_market_inform.TSO_Market.flex_stat.BESS_soc.soc_max = Eigen::MatrixXd::Zero(process_par.time_boundary[1], Power_market_inform.TSO_Market.num_zone);
        Power_market_inform.TSO_Market.flex_stat.BESS_soc.soc_current = .5 * Power_market_inform.TSO_Market.flex_stat.BESS_soc.soc_min.col(0);
        Power_market_inform.TSO_Market.flex_stat.BESS_soc.soc_current += .5 * Power_market_inform.TSO_Market.flex_stat.BESS_soc.soc_min.col(0);
        Power_market_inform.TSO_Market.flex_stat.BESS_soc.capacity_max = Eigen::VectorXd::Zero(Power_market_inform.TSO_Market.num_zone);

        Power_market_inform.TSO_Market.flex_stat.EV_soc.soc_min = Eigen::MatrixXd::Zero(process_par.time_boundary[1], Power_market_inform.TSO_Market.num_zone);
        Power_market_inform.TSO_Market.flex_stat.EV_soc.soc_max = Eigen::MatrixXd::Zero(process_par.time_boundary[1], Power_market_inform.TSO_Market.num_zone);
        Power_market_inform.TSO_Market.flex_stat.EV_soc.soc_current = .5 * Power_market_inform.TSO_Market.flex_stat.EV_soc.soc_min.col(0);
        Power_market_inform.TSO_Market.flex_stat.EV_soc.soc_current += .5 * Power_market_inform.TSO_Market.flex_stat.EV_soc.soc_min.col(0);
        Power_market_inform.TSO_Market.flex_stat.EV_soc.capacity_max = Eigen::VectorXd::Zero(Power_market_inform.TSO_Market.num_zone);
    }

    void contingency_analysis_set(contingency_analysis_struct &contingency_analysis, power_market::market_whole_inform &Power_market_inform, configuration::process_config &process_par){
        // Initialization of number of components
        // index the components in the following order:
        // transformers (transmission nodes), power lines, HV hydro, HV_wind, HV_PSPP
        contingency_analysis.num_component  = Power_market_inform.TSO_Market.network.num_vertice;
        contingency_analysis.num_component += Power_market_inform.TSO_Market.network.num_edges;
        contingency_analysis.num_component += Power_market_inform.agent_profiles.power_supplier.hydro.HV_hybrid.size();
        contingency_analysis.num_component += Power_market_inform.agent_profiles.power_supplier.hydro.HV_plant.size();
        contingency_analysis.num_component += Power_market_inform.agent_profiles.power_supplier.wind.HV_hybrid.size();
        contingency_analysis.num_component += Power_market_inform.agent_profiles.power_supplier.wind.HV_plant.size();
        contingency_analysis.num_component += Power_market_inform.agent_profiles.power_supplier.pump_storage.HV.size();
        int num_component = contingency_analysis.num_component;

        // Initialization of duration
        contingency_analysis.duration = process_par.time_boundary[1];

        // Initialization of temporal probability
    	Eigen::Vector2d transition_prob;
    	transition_prob << 1E-4, .1;
    	contingency_analysis.temporal_prob_0 = Eigen::MatrixXd (num_component, 2);
        contingency_analysis.temporal_prob_0.col(0) << transition_prob[0] * Eigen::VectorXd::Ones(num_component);
        contingency_analysis.temporal_prob_0.col(1) << transition_prob[1] * Eigen::VectorXd::Ones(num_component);

        // Initialization of the LP Problem
    }

    // Contingency sampling using MCMC
    void contigency_sampling(contingency_analysis_struct &contingency_analysis, int num_sample){
        contingency_analysis.num_sample = num_sample;

        int max_int = 1E12;
	    std::random_device rd;  // a seed source for the random number engine
	    std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution <> distrib(0, max_int);

        contingency_analysis.samples_set();
        Eigen::MatrixXi sample = Eigen::MatrixXi::Zero(contingency_analysis.num_component, contingency_analysis.duration);

        // Transform temporal conditional probabilities to hamiltonian
        contingency_analysis.hamiltonian_get();

        // Burn-up step
        for(int burn_up_iter = 0; burn_up_iter < contingency_analysis.num_burn_up; ++ burn_up_iter){
        	local::Gibbs_sampling(contingency_analysis.duration, contingency_analysis.beta, sample, contingency_analysis.temporal_hamiltonian, contingency_analysis.temporal_prob_0, contingency_analysis.spatial_hamiltonian);
		}

		// Actual sampling step (assuming stationary probability has been reached)
		#pragma omp parallel
		{
			#pragma omp for
			for(int sample_iter = 0; sample_iter < contingency_analysis.num_sample; ++ sample_iter){
				local::Gibbs_sampling(contingency_analysis.duration, contingency_analysis.beta, sample, contingency_analysis.temporal_hamiltonian, contingency_analysis.temporal_prob_0, contingency_analysis.spatial_hamiltonian);
				contingency_analysis.samples[sample_iter] = sample;
			}
		}
    }

    // Optimal power flow for different contingencies
    void contingency_analysis_solve(contingency_analysis_struct &contingency_analysis, power_market::market_whole_inform &Power_market_inform, configuration::process_config &process_par){
        int sample_ID = 0;
        int tick = 0;
        local::contingency_analysis_LP_set(sample_ID, tick, contingency_analysis, Power_market_inform.TSO_Market);
    }

}

//int main(){
//	int num_component = 1;
//	int duration = 1000;
//	int num_sample = 100000;
//	int num_burn_up = 1000;
//	double beta = 1.;
//	Eigen::Vector2d transition_prob;
//	transition_prob << 1E-4, .1;
//
//	Eigen::MatrixXd spatial_hamiltonian;
//	Eigen::MatrixXd temporal_prob(num_component, 2);
//	temporal_prob.col(0) << transition_prob[0] * Eigen::VectorXd::Ones(num_component);
//	temporal_prob.col(1) << transition_prob[1] * Eigen::VectorXd::Ones(num_component);
//
//	std::vector <Eigen::MatrixXi> MCMC_sampling = local::sampling(num_sample, duration, beta, temporal_prob, spatial_hamiltonian, num_burn_up);
//
//	//std::cin.get();
//    return 0;
//}

