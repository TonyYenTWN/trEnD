// Header file for contingency analysis
#pragma once

#include <cmath>
#include <filesystem>
#include <random>
#include "src/power_market/power_market.h"  // need to read flex_stat so power_market header is needed

namespace power_network{
    struct contingency_analysis_struct{
        // order of rows: transmission lines, large hydro power plants, large wind power plants, large PSPP

        // Input variables
        int num_burn_up = 10000; //default
        int num_component;
        int num_sample = 100000;
        int duration;
        double beta = 1.;
        Eigen::MatrixXd temporal_prob_0;
        Eigen::MatrixXd temporal_hamiltonian;         // translate probability to energy in the hamiltonian
        Eigen::MatrixXd temporal_energy_variation;  // add this for time-dependent failure energy (h) and correlation energy (j) in the Ising model
        Eigen::MatrixXd spatial_hamiltonian;             // add this for correlation of failure between components

        // Process variables
        std::vector <Eigen::MatrixXd> samples;
        alglib::minlpstate problem;

        // Output variable
        std::vector <Eigen::MatrixXd> energy_not_served;  // ENS  of each node during the entire time interval, for each simulation
        Eigen::MatrixXd energy_not_served_mean;               // Averaged ENS of each node at each time

        // Functions
        void samples_set(){
            this->samples = std::vector <Eigen::MatrixXd> (this->num_sample);
        }

        void hamiltonian_get(){
            this->temporal_hamiltonian = Eigen::MatrixXd(this->temporal_prob_0.rows(), this->temporal_prob_0.cols());
            for(int row_iter = 0; row_iter < temporal_hamiltonian.rows(); ++ row_iter){
                Eigen::Vector2d transition_prob = this->temporal_prob_0.row(row_iter);
                this->temporal_hamiltonian(row_iter, 0) = -std::log(transition_prob[0] * transition_prob[1] / std::pow(1 - transition_prob[0], 2)) / this->beta;
                this->temporal_hamiltonian(row_iter, 1) = std::log((1 - transition_prob[0]) * (1 - transition_prob[1]) / transition_prob[0] / transition_prob[1]) / this->beta;
            }
        }
    };

    // functions
    void flex_stat_input(power_market::market_whole_inform&, power_network::network_inform&, configuration::process_config&);
    void contingency_analysis_set(contingency_analysis_struct&, power_market::market_whole_inform&, configuration::process_config&);
    void contigency_sampling(contingency_analysis_struct&, int, int, configuration::process_config&);
    void contingency_analysis_solve(contingency_analysis_struct&, power_market::market_whole_inform&, power_network::network_inform&, configuration::process_config&);
    void contingency_analysis_print(contingency_analysis_struct&, power_market::market_whole_inform&, configuration::process_config&);
}

