#pragma once

#ifndef DATA_HANDLING
#define DATA_HANDLING

// Including necessary external and internal libraries
#include <stdlib.h>
#include <cfloat>
#include <iostream>
#include <Eigen/Eigen>
#include "DataHandle.h"

using namespace Eigen;
using namespace std;

namespace DataHandle {

    extern
        // Function to compute empirical variance-covariance matrix of the data
        MatrixXd Cov(const MatrixXd&);

    // Function to compute empirical correlation matrix of the data
    MatrixXd Corr(const MatrixXd&);

    // Function to compute empirical correlation between two vectors
    double Corr(const VectorXd&, const VectorXd&);

    // Overloaded function for empirical correlation between two matrix columns
    double Corr(const Block<MatrixXd, -1, 1, true>&, const Block<MatrixXd, -1, 1, true>&);

    // Function to demean each column in a matrix (in-place modification)
    void Demean(MatrixXd&);

    // Function to demean a vector (in-place modification)
    void Demean(VectorXd&);

    // Function to standardize each column in a matrix (in-place modification)
    void Standardise(MatrixXd&);

    // Function to standardize a vector (in-place modification)
    void Standardise(VectorXd&);

    // Overload for standardizing a matrix block (in-place modification)
    void Standardise(Block<MatrixXd, -1, -1, false>&);

    // Function to estimate the variance of a vector
    double Var(const VectorXd&);

    // Function for Kronecker product computation
    MatrixXd KroneckerProd(const MatrixXd&, const MatrixXd&);

    // Function to remove a specific row from a matrix (in-place modification)
    void removeRow(MatrixXd&, const int&, const bool& conservative = 1);

    // Function to remove a specific column from a matrix (in-place modification)
    void removeCol(MatrixXd&, const int&, const bool& conservative = 1);

    // Function to remove a specific element from a vector (in-place modification)
    void removeElement(VectorXd&, const int&, const bool& conservative = 1);
};
#endif /* defined(DATA_HANDLING) */
