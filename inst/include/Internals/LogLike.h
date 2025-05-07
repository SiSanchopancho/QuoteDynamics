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

#pragma once

#ifndef LOGLIKE
#define LOGLIKE

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#include <stdlib.h>
#include <cfloat>
#include <iostream>
#include <RcppEigen.h>

namespace LogLike {

  /** nNegative log-ikelihood function of the quote dynamics model using a multivariate Kalman filter for the estimation of the latent state
  *
  * @param param Parameter vector
  * @param X Observation matrix
  * @param tau Vectzor of durations between the quotes
  * @return ind_matrix Tick indicator matrix
  */
  double MVLL(Eigen::VectorXd param, Eigen::MatrixXd X, Eigen::VectorXd tau, Eigen::MatrixXd ind_matrix);

};
#endif /* defined(LOGLIKE) */
