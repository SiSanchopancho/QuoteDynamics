#pragma once

#ifndef LOGLIKE
#define LOGLIKE

// Define PI for compatibility with various compilers
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// Include standard and Eigen libraries
#include <stdlib.h>
#include <cfloat>
#include <iostream>
#include <Eigen/Eigen>

// Include internal library headers
#include "DataHandle.h"
#include "Filtering.h"

// Use namespaces to simplify code readability
using namespace Eigen;
using namespace DataHandle;
using namespace Filtering;

// The LogLike namespace encapsulates functions related to log-likelihood calculations
namespace LogLike {

    // Calculate log-likelihood based on provided parameters and data
    double LL(VectorXd param, MatrixXd X, VectorXd tau, double m, unsigned d);

};
#endif /* defined(LOGLIKE) */
