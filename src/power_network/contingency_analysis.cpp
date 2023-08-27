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
}

namespace power_network{
    void contingency_analysis_set(contingency_analysis_struct &contingency_analysis, power_market::market_whole_inform &Power_market_inform, configuration::process_config &process_par){
        // Initialization of number of components
        // index the components in the following order!!
        contingency_analysis.num_component = Power_market_inform.TSO_Market.network.num_edges;
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

