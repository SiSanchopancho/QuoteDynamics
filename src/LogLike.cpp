/** SPDX-License-Identifier: GPL-3.0-or-later
*
* Copyright (c) 2025 Domenic Franjic
*
* This file is part of QuoteDynamics.
*
* QuoteDynamics is free software: you can redistribute it
* and/or modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3
* of the License, or (at your option) any later version.
*
* QuoteDynamics is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty
* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
* See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with QuoteDynamics.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "Internals/LogLike.h"

double LogLike::MVLL(
  Eigen::VectorXd param,
  Eigen::MatrixXd X,
  Eigen::VectorXd tau,
  Eigen::MatrixXd ind_matrix
) {

  /* Initialisation */

  // Dimensions
  const unsigned M = X.cols() / 2, T = X.rows(), N = X.cols(), d = 2;

  // Spread vector
  Eigen::VectorXd spread = param.head(M);
  Eigen::VectorXd obs_intercept(2 * M);
  obs_intercept(Eigen::seq(0, 2 * M - 1, 2)) = -0.5 * spread;
  obs_intercept(Eigen::seq(1, 2 * M - 1, 2)) = 0.5 * spread;

  // Loadings
  Eigen::VectorXd loading_base = param(Eigen::seq(M, (3 * M) - 1));

  /* Initialisation of the variance-covariance matrix of the data */

  Eigen::MatrixXd lower_chol_factor = Eigen::MatrixXd::Zero(2, 2); // Dummy matrix
  Eigen::MatrixXd obs_var_cov = Eigen::MatrixXd::Zero(N, N); // Actual var.-cov. matrix of the data (Block diagonal structure with each of the M diagonal 2x2 block representing to the var.-cov. matrix between bid and ask within the correpsonding market)
  int ind = 0;
  for (unsigned m = 0; m < M; ++m) {

    lower_chol_factor.setZero();

    for (int i = 0; i < lower_chol_factor.rows(); i++) {

      for (int j = 0; j <= i; j++) {

        lower_chol_factor(i, j) = param(3 * M + ind);
        ++ind;

      }
    }
    obs_var_cov.block(2 * m, 2 * m, 2, 2) = lower_chol_factor * lower_chol_factor.transpose();
  }


  Eigen::VectorXd state_intercept = Eigen::VectorXd::Zero(2);
  Eigen::MatrixXd state_trans = Eigen::MatrixXd::Zero(2, 2);
  state_trans(0, 0) = 1.0;
  Eigen::MatrixXd Tt_transpose = state_trans.transpose();


  const double sigma = param(Eigen::last - 2);
  const double delta1 = param(Eigen::last - 1);
  const double delta2 = param(Eigen::last);

  /* Initialisation of the Kalman filter */

  // State moments
  Eigen::VectorXd state_mean = Eigen::VectorXd::Zero(d);
  Eigen::MatrixXd state_var_cov = 1000 * Eigen::MatrixXd::Identity(d, d);

  // Measurement equation
  Eigen::MatrixXd loading_matrix(loading_base.size(), d);
  loading_matrix.col(0).setConstant(1.0);
  loading_matrix.col(1) = loading_base * (std::pow(tau(0), delta2) * sigma);
  Eigen::MatrixXd loading_matrix_transpose = loading_matrix.transpose();
  Eigen::VectorXd cond_obs_mean = obs_intercept + loading_matrix * state_mean;
  Eigen::MatrixXd cond_obs_var = loading_matrix * state_var_cov * loading_matrix_transpose + obs_var_cov;

  // Additional dummies
  Eigen::VectorXd fc_error = X.row(0).transpose() - cond_obs_mean;
  Eigen::VectorXd lnl = Eigen::VectorXd::Zero(T);
  Eigen::MatrixXd selection_matrix = Eigen::MatrixXd::Identity(M, M);
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(d, d);

  // Initialisation of the ll
  const double log_two_pi = std::log(2 * M_PI);
  Eigen::LLT<Eigen::MatrixXd> llt(cond_obs_var);
  Eigen::MatrixXd cond_obs_var_inv = llt.solve(Eigen::MatrixXd::Identity(cond_obs_var.rows(), cond_obs_var.rows()));
  Eigen::MatrixXd cholesky_factor = llt.matrixL();
  double cond_obs_var_log_det = 2.0 * cholesky_factor.diagonal().array().log().sum();
  lnl(0) = -0.5 * N * log_two_pi - 0.5 * cond_obs_var_log_det - 0.5 * (fc_error.transpose() * cond_obs_var_inv * fc_error)(0, 0);

  /* Update */

  Eigen::MatrixXd kalman_gain = state_var_cov * loading_matrix_transpose * cond_obs_var_inv;
  Eigen::VectorXd updt_state_mean = state_mean + kalman_gain * fc_error;
  Eigen::MatrixXd updt_state_var_cov = state_var_cov - state_var_cov * loading_matrix_transpose * kalman_gain.transpose();

  /* Start filter loop */

  for (unsigned t = 1; t < T; ++t) {

    /* Construct selection Matrix */

    selection_matrix.conservativeResize(2 * ind_matrix.row(t).sum(), 2 * M);
    selection_matrix.setZero();
    unsigned row = 0;
    for (unsigned m = 0; m < M; ++m) {
      if (ind_matrix(t, m) == 1) {
        selection_matrix(row, 2 * m) = 1.0; // Bid
        selection_matrix(row + 1, 2 * m + 1) = 1.0; // Ask
        row += 2;
      }
    }

    /* Update the loading matrix */

    Eigen::MatrixXd loading_matrix_selected = loading_matrix;
    loading_matrix_selected.col(0).setConstant(1.0);
    loading_matrix_selected.col(1) = loading_base * (std::pow(tau(t), delta2) * sigma);
    loading_matrix_selected = (selection_matrix * loading_matrix_selected);
    Eigen::MatrixXd loading_matrix_selected_transpose = loading_matrix_selected.transpose();
    H.setZero();
    H(0, 0) = sigma * std::pow(tau(t), delta1);
    H(d - 1, d - 1) = 1;

    /* Predict */

    state_mean = state_intercept + state_trans * updt_state_mean;
    state_var_cov = state_trans * updt_state_var_cov * Tt_transpose + H * H.transpose();

    /* Evaluate */

    cond_obs_mean = selection_matrix * obs_intercept + loading_matrix_selected * state_mean;
    cond_obs_var = loading_matrix_selected * state_var_cov * loading_matrix_selected_transpose + selection_matrix * obs_var_cov * selection_matrix.transpose();
    fc_error = selection_matrix * X.row(t).transpose() - cond_obs_mean;

    /* Compute ll */;

    Eigen::LLT<Eigen::MatrixXd> llt_temp(cond_obs_var);
    cond_obs_var_inv = llt_temp.solve(Eigen::MatrixXd::Identity(cond_obs_var.rows(), cond_obs_var.rows()));
    cholesky_factor = llt_temp.matrixL();
    cond_obs_var_log_det = 2.0 * cholesky_factor.diagonal().array().log().sum();
    lnl(t) = -0.5 * N * log_two_pi - 0.5 * cond_obs_var_log_det - 0.5 * (fc_error.transpose() * cond_obs_var_inv * fc_error)(0, 0);


    /* Update */

    kalman_gain = state_var_cov * loading_matrix_selected_transpose * cond_obs_var_inv;
    updt_state_mean = state_mean + kalman_gain * fc_error;
    updt_state_var_cov = state_var_cov - state_var_cov * loading_matrix_selected_transpose * kalman_gain.transpose();

  }


  /* End filter loop*/

  return -lnl.mean();
}
