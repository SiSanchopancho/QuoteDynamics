#pragma once

#ifndef ORDERS
#define ORDERS

#define _USE_MATH_DEFINES

// Include standard, Eigen, and other necessary libraries
#include <stdlib.h>
#include <random>
#include <cfloat>
#include <iostream>
#include <Eigen/Eigen>
#include <math.h>

// Include internal library headers
#include "DataHandle.h"

// Use namespaces for convenience
using namespace Eigen;
using namespace DataHandle;

// The Orders namespace contains functions for determining orders in statistical models
namespace Orders {

    // Determine the number of degrees of freedom for a given dataset and parameters
    extern double NoOfDF(const MatrixXd& X_in, const unsigned& k0, const unsigned& k1, const double& w0, const VectorXi& s);

    // Calculate the number of significant factors in a dataset
    double NoOfSF(const MatrixXd& X_in, const unsigned& k0, const unsigned& k1);

    // Overloaded functions to determine the number of factors in a dataset
    int NoOfFactors(const MatrixXd& X, const unsigned& k0, const unsigned& k1, const double& w0, const VectorXi& s, const double& alpha = 0.01, const bool& log = 1);
    int NoOfFactors(const MatrixXd& X, const unsigned& k0, const unsigned& k1, const double& alpha = 0.01, const bool& log = 1);

    // Determine the order of a VAR process based on given criteria
    int VARorder(const MatrixXd& F, const int& O, const double& comp_null = 10e-15, const char* crit = "BIC");
};
#endif /* defined(ORDERS) */
