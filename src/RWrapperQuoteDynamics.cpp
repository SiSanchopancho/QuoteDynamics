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

#include <iomanip> 
#include <vector>
#include <climits>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <random>
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
	Eigen::VectorXd day_indicator;
	int eval_count;
	double initial_kf_var_cov;
};

/** Wrapper template to parse results back into R */
namespace Rcpp {
	template <>
	SEXP wrap(const Results& x) {
		NumericVector estimate = Rcpp::wrap(x.estimate);
		NumericMatrix hessian = Rcpp::wrap(x.hessian);
		NumericMatrix boot_straps = Rcpp::wrap(x.boot_straps);
		return List::create(Named("min_val") = x.min_val,
			Named("estimate") = estimate,
			Named("hessian") = hessian,
			Named("nlopt_return_code") = x.nlopt_return_code,
			Named("eval_count") = x.eval_count,
			Named("boot_straps") = boot_straps);
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
	double ll = MVLL(parameters, data_in->X, data_in->tau, data_in->ind_matrix, data_in->day_indicator, data_in->initial_kf_var_cov);
	data_in->eval_count++;

	if (data_in->log) {
		if (data_in->eval_count % 10 == 0) {
			Rcpp::Rcout << "\r\033[KCurrent objective value is " << std::scientific << std::setprecision(5) << ((ll < 0) ? "-" : " ") << std::abs(ll) << " after " << data_in->eval_count << " function evaluations.";
			Rcpp::Rcout.flush();
		}
	}

	// Return the negative log-likelihood (as NLopt minimizes by default)
	return ll;
}


//'  Wrapper to compute the neg-log-like ought to be called from within R for usage with Rs optim()
//'
//' @param n Size of the paramter vector
//' @param x Pointer to the data of the parameter vector
//' @param grad Pointer to the data of the gradient
//' @param data Pointer to the data structure that holds additional function parameters of the actual objective function
//' @return Objective function value of the current parameter estimate
//' @export
// [[Rcpp::export]]
double objFunctionCpp(NumericVector startR, NumericMatrix XR, NumericVector tauR, NumericMatrix ind_matrix_R, NumericMatrix day_indicator_R, const double& initial_kf_var_cov, bool log)
{

	Eigen::VectorXd parameters(Eigen::VectorXd::Map(startR.begin(), startR.size()));
	Eigen::Map<Eigen::MatrixXd> X(Eigen::MatrixXd::Map(XR.begin(), XR.rows(), XR.cols()));
	Eigen::Map<Eigen::VectorXd> tau(Eigen::VectorXd::Map(tauR.begin(), tauR.size()));
	Eigen::Map<Eigen::MatrixXd> ind_matrix(Eigen::MatrixXd::Map(ind_matrix_R.begin(), ind_matrix_R.rows(), ind_matrix_R.cols()));
	Eigen::Map<Eigen::VectorXd> day_indicator(Eigen::VectorXd::Map(day_indicator_R.begin(), day_indicator_R.size()));


	// Compute the LL value using parameters as the model parameters
	double ll = MVLL(parameters, X, tau, ind_matrix, day_indicator, initial_kf_var_cov);

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

			double f1 = MVLL(parameters_plus, data_in->X, data_in->tau, data_in->ind_matrix, data_in->ind_matrix, data_in->initial_kf_var_cov);
			double f2 = MVLL(parameters, data_in->X, data_in->tau, data_in->ind_matrix, data_in->ind_matrix, data_in->initial_kf_var_cov);
			double f3 = MVLL(parameters_minus, data_in->X, data_in->tau, data_in->ind_matrix, data_in->ind_matrix, data_in->initial_kf_var_cov);
			Hessian(i, j) = Hessian(j, i) = (f1 - 2 * f2 + f3) / (h * h);

		}

	}
	return Hessian;
}

/** Progress-bar */
void printProgress(unsigned current, unsigned total, int bar_width = 40) {

	// Compute progress
	double progress = double(current) / double(total);
	int pos = int(bar_width * progress);
	
	// Escape character used to flush the console without deleting earlier output
	Rcpp::Rcout << "\r\033[K[";

	// Set up the filled part
	for (int i = 0; i < bar_width; ++i) {
		if (i < pos) { 
			Rcpp::Rcout << "="; 
		}
		else if (i == pos) { 
			Rcpp::Rcout << ">"; 
		}
		else { Rcpp::Rcout << " "; }
	}

	// Percentage
	Rcpp::Rcout << "] " << std::setw(3) << std::fixed << std::setprecision(0) << (progress * 100) << "%";
	Rcpp::Rcout.flush();
}

/** Check whether the current algorithm is global */
bool isGlobal(nlopt_algorithm algo) {
	return algo == NLOPT_GN_DIRECT ||
		algo == NLOPT_GN_DIRECT_L ||
		algo == NLOPT_GN_DIRECT_L_RAND ||
		algo == NLOPT_GN_DIRECT_NOSCAL ||
		algo == NLOPT_GN_DIRECT_L_NOSCAL ||
		algo == NLOPT_GN_DIRECT_L_RAND_NOSCAL ||
		algo == NLOPT_GN_ORIG_DIRECT ||
		algo == NLOPT_GN_ORIG_DIRECT_L ||
		algo == NLOPT_GN_CRS2_LM ||
		algo == NLOPT_GN_ISRES ||
		algo == NLOPT_GN_ESCH ||
		algo == NLOPT_GN_MLSL ||
		algo == NLOPT_GN_MLSL_LDS;
}

//' FastOptim Function
//'
//' Fast KFS MLE implementation
//'
//' @param startR Starting values
//' @param XR Data Matrix
//' @param tauR Parameter vector
//' @param ind_matrix_R Identifier matrix
//' @param rel_xtol Algorithm tolerance
//' @param rel_ftol Stopping rule
//' @param algorithm_id NLopt algorithm ID
//' @param upper_bound_R Upper parameter bounds (only for global optimisers)
//' @param lower_bound_R Lower parameter bounds (only for global optimisers)
//' @param initial_kf_var_cov Initial Kalman filter state variance
//' @param compute_se Whether or not standard eroors should be computed
//' @param no_of_bootstraps Number of bootstrap repititions
//' @param length_of_bootstraps Length of each bootstrap block
//' @param seed Seed used for bootstrapping
//' @param hessian Whether or not to compute the hessian
//' @param step_size Step size for the computation of the hessian
//' @param log Whether or not to print output
//' @return Returns optimisation object
//' @export
// [[Rcpp::export]]
Results FastOptim(
	NumericVector startR,
	NumericMatrix XR,
	NumericVector tauR,
	NumericMatrix ind_matrix_R,
	NumericMatrix day_indicator_R,
	double rel_xtol,
	double rel_ftol,
	int max_eval,
	int algorithm_id,
	NumericVector upper_bound_R,
	NumericVector lower_bound_R,
	double initial_kf_var_cov,
	bool compute_se,
	int no_of_bootstraps,
	int length_of_bootstraps,
	int seed,
	bool hessian,
	double step_size,
	bool log
) {

	/* Initialisation */

	// Convert R types to Eigen types
	Eigen::VectorXd estimate(Eigen::VectorXd::Map(startR.begin(), startR.size()));
	Eigen::Map<Eigen::MatrixXd> X(Eigen::MatrixXd::Map(XR.begin(), XR.rows(), XR.cols()));
	Eigen::Map<Eigen::VectorXd> tau(Eigen::VectorXd::Map(tauR.begin(), tauR.size()));
	Eigen::Map<Eigen::MatrixXd> ind_matrix(Eigen::MatrixXd::Map(ind_matrix_R.begin(), ind_matrix_R.rows(), ind_matrix_R.cols()));
	Eigen::Map<Eigen::VectorXd> day_indicator(Eigen::VectorXd::Map(day_indicator_R.begin(), day_indicator_R.size()));


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
	data.eval_count = 0;
	data.day_indicator = day_indicator;
	data.initial_kf_var_cov = initial_kf_var_cov;

	/* NLopt Set-Up */

	nlopt_opt opt; // Create object
	opt = nlopt_create(static_cast<nlopt_algorithm>(algorithm_id), estimate.size()); // Set the optimisation algorithm and the dimensionality (number of parameters) of the problem
	std::vector<double> initial_guess(estimate.data(), estimate.data() + estimate.size()); // Map the Eigen::VectorXd stored data into a std::vector for NLopt
	nlopt_set_min_objective(opt, objective_function, &data); // Creat the minimisation problem
	nlopt_set_xtol_rel(opt, rel_xtol); // Set step-tolerance
	nlopt_set_ftol_rel(opt, rel_ftol); // Set obj. function value tolerance
	nlopt_set_maxeval(opt, max_eval); // Set maximum number of obj. function evaluation
	double minf; // Minimum obj. value at return

	if (isGlobal(static_cast<nlopt_algorithm>(algorithm_id))) {
		std::vector<double> upper_bound(upper_bound_R.begin(), upper_bound_R.begin() + upper_bound_R.size());
		std::vector<double> lower_bound(lower_bound_R.begin(), lower_bound_R.begin() + lower_bound_R.size());
		nlopt_set_upper_bounds(opt, upper_bound.data());
		nlopt_set_lower_bounds(opt, lower_bound.data());
	}

	/* Optimisation */

	if(log) { Rcpp::Rcout << "Start optimisation." << std::endl; }
	nlopt_result optim_outcome = nlopt_optimize(opt, estimate.data(), &minf);
	if(log) {
		Rcpp::Rcout << std::defaultfloat << std::setprecision(6);  // Reset from scientific notation back to default
		if (optim_outcome > 0) {
			Rcpp::Rcout << std::endl << "Optimisation exited succesfully." << std::endl;
		}
		else {
			Rcpp::Rcout << std::endl << "Optimisation did not exit succesfully." << std::endl;
		}
	}

	/* Create an output object */

	Results out;
	out.min_val = minf;
	out.estimate = estimate;
	out.nlopt_return_code = optim_outcome;
	out.eval_count = data.eval_count;

	if (hessian) {// Compute hessian
		out.hessian = ComputeHessian(out.estimate, &data, step_size);
	}

	nlopt_destroy(opt); // Destroy NLopt object properly!

	/* Start bootstrap standard error block */

	if (compute_se) {

		// Initialise progress bar
		printProgress(0, no_of_bootstraps);

		// Create output matrix
		unsigned no_of_bs_blocks = X.rows() - length_of_bootstraps + 1;
		Eigen::MatrixXd bootstrapped_estimates = Eigen::MatrixXd::Zero(estimate.size(), no_of_bootstraps);

		// Set-up RNG
		std::mt19937 gen(seed);
		std::uniform_int_distribution<int> dist(0, no_of_bs_blocks - 1);

		/* Create bootstrap block indices */

		// Populate the dummie matrix holding the time indices of the corresponding blocks
		Eigen::MatrixXi Bootstrap_blocks = Eigen::MatrixXi::Zero(length_of_bootstraps, no_of_bs_blocks);
		for (unsigned l = 0; l < no_of_bs_blocks; ++l) {
			Bootstrap_blocks.col(l) = Eigen::VectorXi::LinSpaced(length_of_bootstraps, l, l + length_of_bootstraps - 1);
		}

		/* Start main bootstrap loop */

		for (unsigned boot = 0; boot < no_of_bootstraps; ++boot) {

			// Store an initial result vector
			Eigen::VectorXd current_guess = estimate;

			// Construct bootstrapped data
			unsigned bootstrap_length_cap = std::ceil(double(X.rows()) / double(length_of_bootstraps));
			Eigen::MatrixXd X_bootstrapped = Eigen::MatrixXd::Zero(length_of_bootstraps * bootstrap_length_cap, X.cols());
			for (unsigned m = 0; m < bootstrap_length_cap; ++m) {
				int draw = dist(gen);
				//X_bootstrapped.block(m * length_of_bootstraps, 0, length_of_bootstraps, X.cols()) = X(Bootstrap_blocks.col(draw), Eigen::all);
				for (int k = 0; k < length_of_bootstraps; ++k) {
					X_bootstrapped.row(m * length_of_bootstraps + k) = X.row(Bootstrap_blocks(k, draw));
				}
			}
			X_bootstrapped.conservativeResize(X.rows(), X.cols());

			/* NLopt Set-Up */

			OptimData data_bs;
			data_bs.X = X_bootstrapped;
			data_bs.tau = tau;
			data_bs.log = false;
			data_bs.ind_matrix = ind_matrix;
			data_bs.eval_count = 0;
			data_bs.day_indicator = day_indicator;

			nlopt_opt opt_bs; // Create object
			opt_bs = nlopt_create(NLOPT_LN_SBPLX, current_guess.size()); // Set the optimisation algorithm and the dimensionality (number of parameters) of the problem
			std::vector<double> initial_guess(current_guess.data(), current_guess.data() + current_guess.size()); // Map the Eigen::VectorXd stored data into a std::vector for NLopt
			nlopt_set_min_objective(opt_bs, objective_function, &data_bs); // Creat the minimisation problem
			nlopt_set_xtol_rel(opt_bs, rel_xtol); // Set step-tolerance
			nlopt_set_ftol_rel(opt_bs, rel_ftol); // Set obj. function value tolerance
			nlopt_set_maxeval(opt_bs, max_eval); // Set maximum number of obj. function evaluation
			double minf_bs; // Minimum obj. value at return
			nlopt_optimize(opt_bs, current_guess.data(), &minf_bs);
			bootstrapped_estimates.col(boot) = current_guess;
			nlopt_destroy(opt_bs); // Destroy NLopt object properly!

			printProgress(boot + 1, no_of_bootstraps);

		}

		/* End main bootstrap loop */

		// Store results
		out.boot_straps = bootstrapped_estimates;

	}

	/* End bootstrap standard error block */

	return out;
}
