#pragma once

#ifndef ALGORITHM
#define ALGORITHM
#include <algorithm>
#endif

#ifndef CONSTANT
#define CONSTANT
#include <boost/math/constants/constants.hpp>
#endif

#ifndef EIGEN
#define EIGEN
#include <Eigen/Dense>
#include <Eigen/Core>
namespace Eigen{
	typedef Eigen::Matrix<ptrdiff_t, Eigen::Dynamic, 1> VectorXpd;
}
#endif

#ifndef MAP
#define MAP
#include <map>
#endif

#ifndef NUMERIC
#define NUMERIC
#include <numeric>
#endif

#ifndef OMP
#define OMP
#include <omp.h>
#endif

#ifndef VECTOR
#define VECTOR
#include <vector>
#endif
