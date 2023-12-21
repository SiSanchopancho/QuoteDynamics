#pragma once

#ifndef KFSMLE_TYPES
#define KFSMLE_TYPES

// Include necessary Rcpp and Eigen libraries
#include <RcppCommon.h>
#include <Rcpp.h>
#include <Eigen/Eigen>

// Results structure to store the outcomes of the optimization process
struct Results {
    double min_val;           // Minimum value found by the optimization
    Eigen::VectorXd estimate; // Estimated parameters
    Eigen::MatrixXd hessian;  // Hessian matrix at the optimum
};

// Rcpp namespace to enable seamless integration with R
namespace Rcpp {
    // Specialize the wrap function for the Results struct
    template <>
    SEXP wrap(const Results& x);
}

#endif /* defined(KFSMLE_TYPES) */
