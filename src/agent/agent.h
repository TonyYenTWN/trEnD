// Header file for agents
#pragma once
#include "src/alglib/optimization.h"
#include "src/basic/basic_definitions.h"
#include "src/basic/eigen_sparse.h"

namespace agent{
	namespace parameters{
		static inline int foresight_time(){
			int value = 24;
			return value;
		}

		static inline int sample_num(){
			int value = 3;
			return value;
		}

		static inline double residential_ratio(){
			//double value = 1.;
			double value = .4;
			return value;
		}
	}

	struct sorted_vector{
		Eigen::VectorXi id;
		Eigen::VectorXd value;
	};

	// Functions
	sorted_vector sort(Eigen::VectorXd);
}
