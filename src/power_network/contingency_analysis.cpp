#include "contingency_analysis.h"

namespace local{
	static inline Eigen::Vector2d hamiltonian_get(Eigen::Vector2d transition_prob, double beta){
        Eigen::Vector2d hamiltonian;
		hamiltonian[0] = -std::log(transition_prob[0] * transition_prob[1] / std::pow(1 - transition_prob[0], 2)) / beta;
		hamiltonian[1] = std::log((1 - transition_prob[0]) * (1 - transition_prob[1]) / transition_prob[0] / transition_prob[1]) / beta;

		return hamiltonian;
	}

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

    std::vector <Eigen::MatrixXi> MCMC_sampling(int num_sample, int duration, double beta, Eigen::MatrixXd &temporal_prob,  Eigen::MatrixXd &spatial_hamiltonian, int num_burn_up = 1000){
        int max_int = 1E12;
	    std::random_device rd;  // a seed source for the random number engine
	    std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<> distrib(0, max_int);

		std::vector <Eigen::MatrixXi> samples(num_sample);
        Eigen::MatrixXi sample = Eigen::MatrixXi::Zero(temporal_prob.rows(), duration);
        Eigen::MatrixXd temporal_hamiltonian(temporal_prob.rows(), temporal_prob.cols());

        // Transform temporal conditional probabilities to hamiltonian
        for(int row_iter = 0; row_iter < sample.rows(); ++ row_iter){
        	Eigen::Vector2d hamiltonian = hamiltonian_get(temporal_prob.row(row_iter), beta);
        	temporal_hamiltonian.row(row_iter) = hamiltonian;
		}

        // Burn-up step
        for(int burn_up_iter = 0; burn_up_iter < num_burn_up; ++ burn_up_iter){
        	Gibbs_sampling(duration, beta, sample, temporal_hamiltonian, temporal_prob, spatial_hamiltonian);
		}

		// Actual sampling step (assuming stationary probability has been reached)
		#pragma omp parallel
		{
			#pragma omp for
			for(int sample_iter = 0; sample_iter < num_sample; ++ sample_iter){
				Gibbs_sampling(duration, beta, sample, temporal_hamiltonian, temporal_prob, spatial_hamiltonian);
				samples[sample_iter] = sample;
			}
		}

        return samples;
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

