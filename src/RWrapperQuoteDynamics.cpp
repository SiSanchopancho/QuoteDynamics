//' SPDX-License-Identifier: GPL-3.0-or-later */
//'
//' Copyright (c) 2025 Domenic Franjic
//'
//' This file is part of QuoteDynamics.
//'
//' QuoteDynamics is free software: you can redistribute it
//' and/or modify it under the terms of the GNU General Public License
//' as published by the Free Software Foundation, either version 3
//' of the License, or (at your option) any later version.
//'
//' QuoteDynamics is distributed in the hope that it will be useful,
//' but WITHOUT ANY WARRANTY; without even the implied warranty
//' of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//' See the GNU General Public License for more details.
//'
//' You should have received a copy of the GNU General Public License
//' along with QuoteDynamics.  If not, see <https://www.gnu.org/licenses/>.
//' @docType package
//' @name QuoteDynamics

#define _USE_MATH_DEFINES

// Defin PI when using g++
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define NLOPT_DLL

#include <vector>
#include <climits>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Eigen>
#include <nlopt/src/api/nlopt.h>
#include "QuoteDynamics_types.h"
#include "Internals/LogLike.h"

using namespace LogLike;
using namespace Rcpp;

/** Data structure parsed to NLopt's optimisation routine */
struct OptimData {
    Eigen::MatrixXd X;
    Eigen::VectorXd tau;
    bool log;
    Eigen::MatrixXd ind_matrix;
};

/** Wrapper template to parse results back into R */
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


/**  NLopt objective function wrapper
*
* @param n Size of the paramter vector
* @param x Pointer to the data of the parameter vector
* @param grad Pointer to the data of the gradient
* @param data Pointer to the data structure that holds additional function parameters of the actual objective function
*
* @return Objective function value of the current parameter estimate
* */
double objective_function(unsigned n, const double* x, double* grad, void* data)
{

    // Load in the data
    OptimData* data_in = static_cast<OptimData*>(data);

    // Create an Eigen::VectorXd from the input array x
    Eigen::VectorXd parameters = Eigen::Map<const Eigen::VectorXd>(x, n);

    // Compute the LL value using parameters as the model parameters
    double ll = MVLL(parameters, data_in->X, data_in->tau, data_in->ind_matrix);

    if (data_in->log) {

        Rcpp::Rcout << "Current objective value: " << ll << std::endl;

    }

    // Return the negative log-likelihood (as NLopt minimizes by default)
    return ll;
}

/** Hessian computation
*
* @param parameters Reference to the Eigen::VectorXd that holds the (final) parameter estimates
* @param data Pointer to the data structure that holds additional function parameters of the actual objective function 
* @param h Step-size for computing the Hessian
* @return Eigen::MatrixXd Hessian of the objective function at the estimate
*/
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

            double f1 = MVLL(parameters_plus, data_in->X, data_in->tau, data_in->ind_matrix);
            double f2 = MVLL(parameters, data_in->X, data_in->tau, data_in->ind_matrix);
            double f3 = MVLL(parameters_minus, data_in->X, data_in->tau, data_in->ind_matrix);
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
//' @param ind_matrix_R Identifier matrix
//' @param xtol Algorithm tolerance
//' @param stop_val Stopping rule
//' @return Returns minimum value (changes start in place)
//' @export
// [[Rcpp::export]]
Results FastOptim(NumericVector startR, NumericMatrix XR, NumericVector tauR, NumericMatrix ind_matrix_R, double xtol, double stop_val, int max_eval, int algorithm_id, bool hessian, double step_size, bool log) {
    
    /* Initialisation */

    // Convert R types to Eigen types
    Eigen::VectorXd estimate(Eigen::VectorXd::Map(startR.begin(), startR.size()));
    Eigen::Map<Eigen::MatrixXd> X(Eigen::MatrixXd::Map(XR.begin(), XR.rows(), XR.cols()));
    Eigen::Map<Eigen::VectorXd> tau(Eigen::VectorXd::Map(tauR.begin(), tauR.size()));
    Eigen::Map<Eigen::MatrixXd> ind_matrix(Eigen::MatrixXd::Map(ind_matrix_R.begin(), ind_matrix_R.rows(), ind_matrix_R.cols()));

    // Handle the case where meax_eval has been disabled and set it to INT_MAX (DO NOT FULLY DISABLE max_eval AS SOME ALGORITHMS WILL NOT TERMINATE!)
    if (max_eval == -1) {
      max_eval = INT_MAX;
    }

    // Store all the data in an NLopt-convinient struct
    OptimData data;
    data.X = X;
    data.tau = tau;
    data.log = log;
    data.ind_matrix = ind_matrix;

    /* NLopt Set-Up */

    nlopt_opt opt; // Create object
    opt = nlopt_create(static_cast<nlopt_algorithm>(algorithm_id), estimate.size()); // Set the optimisation algorithm and the dimensionality (number of parameters) of the problem
    std::vector<double> initial_guess(estimate.data(), estimate.data() + estimate.size()); // Map the Eigen::VectorXd stored data into a std::vector for NLopt
    nlopt_set_min_objective(opt, objective_function, &data); // Creat the minimisation problem
    nlopt_set_xtol_rel(opt, xtol); // Set step-tolerance
    nlopt_set_stopval(opt, stop_val); // Set obj. function value tolerance
    nlopt_set_maxeval(opt, max_eval); // Set maximum number of obj. function evaluation
    double minf; // Minimum obj. value at return

    /* Optimisation */

    if (nlopt_optimize(opt, estimate.data(), &minf) < 0) { // Logic block for handling the case where NLopt fails
        Rcpp::stop("NLopt failed!");
    }

    /* Create an output object */

    Results out;
    out.min_val = minf;
    out.estimate = estimate;

    if (hessian) {// Compute hessian
        out.hessian = ComputeHessian(out.estimate, &data, step_size);
    }

    nlopt_destroy(opt); // Destroy NLopt object properly!

    return out;
}
