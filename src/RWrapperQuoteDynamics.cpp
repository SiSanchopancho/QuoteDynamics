#define _USE_MATH_DEFINES // If you need some math constants

// Defin PI when using g++
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define NLOPT_DLL

// Externakl includes

#include <stdlib.h>
#include <cfloat>
#include <iostream>
#include <Eigen/Eigen>
#include <math.h>
#include <nlopt/src/api/nlopt.h>
#include <RcppCommon.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "QuoteDynamics_types.h"

// Internal Incldues

#include "Internals/DataHandle.h" // Data Handlers (Mean, Var, Print Matrices, etc.)
#include "Internals/Filtering.h" // Univariate representataion of multivariate KFS according to Koopman and Durbin (2010)
#include "Internals/Orders.h" // Functions to infer on "orders" of the proces (number of factors, VAR order, etc.)
#include "Internals/LogLike.h" // Functions to infer on "orders" of the proces (number of factors, VAR order, etc.)

using namespace DataHandle;
using namespace Filtering;
using namespace Orders;
using namespace LogLike;
using namespace Rcpp;

struct OptimData {
    Eigen::MatrixXd X;
    Eigen::VectorXd tau;
    bool log;
};


namespace Rcpp {
    template <>
    SEXP wrap(const Results& x) {
        NumericVector estimate = Rcpp::wrap(x.estimate);
        NumericMatrix hessian = Rcpp::wrap(x.hessian);
        return List::create(Named("min_val") = x.min_val,
            Named("estimate") = estimate,
            Named("hessian") = hessian);
    }
}

double objective_function(unsigned n, const double* x, double* grad, void* data)
{

    // Load in the data
    OptimData* data_in = static_cast<OptimData*>(data);

    // Create an Eigen::VectorXd from the input array x
    Eigen::VectorXd parameters = Eigen::Map<const Eigen::VectorXd>(x, n);

    // Compute the LL value using parameters as the model parameters
    double ll = LL(parameters, data_in->X, data_in->tau);

    if (data_in->log) {

        Rcpp::Rcout << "Current objective value: " << ll << std::endl;

    }

    // Return the negative log-likelihood (as NLopt minimizes by default)
    return ll;
}

Eigen::MatrixXd ComputeHessian(const Eigen::VectorXd& parameters, void* data, double h) {

    // Load in the data
    OptimData* data_in = static_cast<OptimData*>(data);

    // Calculate the hessian

    Eigen::MatrixXd Hessian(parameters.size(), parameters.size());

    for (int i = 0; i < parameters.size(); ++i) {

        for (int j = i; j < parameters.size(); ++j) {

            Eigen::VectorXd parameters_plus = parameters, parameters_minus = parameters;

            parameters_plus(i) += h; parameters_plus(j) += h;
            parameters_minus(i) -= h; parameters_minus(j) -= h;

            double f1 = LL(parameters_plus, data_in->X, data_in->tau);
            double f2 = LL(parameters, data_in->X, data_in->tau);
            double f3 = LL(parameters_minus, data_in->X, data_in->tau);
            Hessian(i, j) = Hessian(j, i) = (f1 - 2 * f2 + f3) / (h * h);

        }

    }
    return Hessian;
}

//' FastOptim Function
//'
//' Fast KFS MLE implementation
//'
//' @param startR Starting values
//' @param XR Data Matrix
//' @param tauR Parameter vector
//' @param xtol Algorithm tolerance
//' @param stop_val Stopping rule
//' @return Returns minimum value (changes start in place)
//' @export
// [[Rcpp::export]]
Results FastOptim(NumericVector startR, NumericMatrix XR, NumericVector tauR, double xtol, double stop_val, int algorithm_id, bool hessian, double step_size, bool log) {
    
    // Convert R types to Eigen types
    
    Eigen::VectorXd estimate(Eigen::VectorXd::Map(startR.begin(), startR.size()));
    Eigen::Map<Eigen::MatrixXd> X(Eigen::MatrixXd::Map(XR.begin(), XR.rows(), XR.cols()));
    Eigen::Map<Eigen::VectorXd> tau(Eigen::VectorXd::Map(tauR.begin(), tauR.size()));

    OptimData data;
    data.X = X;
    data.tau = tau;
    data.log = log;

    nlopt_opt opt;

    opt = nlopt_create(static_cast<nlopt_algorithm>(algorithm_id), estimate.size()); /* algorithm and dimensionality */
    std::vector<double> initial_guess(estimate.data(), estimate.data() + estimate.size());

    nlopt_set_min_objective(opt, objective_function, &data);

    nlopt_set_xtol_rel(opt, xtol);

    nlopt_set_stopval(opt, stop_val);

    double minf; /* `*`the` `minimum` `objective` `value,` `upon` `return`*` */
    if (nlopt_optimize(opt, estimate.data(), &minf) < 0) {
        Rcpp::stop("NLopt failed!");
    }

    Results out;
    out.min_val = minf;
    out.estimate = estimate;

    if (hessian)
    {
        // Compute hessian

        out.hessian = ComputeHessian(out.estimate, &data, step_size);

    }

    nlopt_destroy(opt);

    return out;
}
