// Header file for operation and investment of industrial consumers
#include "src/alglib/optimization.h"
#include "src/basic/basic_definitions.h"
#include "src/basic/eigen_sparse.h"

namespace agent{
	namespace industrial{
		static inline double flexible_ratio(){
			double value = .05;
			return value;
		}
	}
}
