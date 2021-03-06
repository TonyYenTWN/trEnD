#ifndef EIGEN_SPARSE
#define EIGEN_SPARSE

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>								// Use only when the matric is symmetric positive definite (ex. covariance or laplacian matrix)
namespace Eigen{
	typedef Eigen::Triplet <double> TripletXd;    			// Define a triplet object
	typedef Eigen::Triplet <std::complex <double>> TripletXcd;
}

#endif
