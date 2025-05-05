#pragma once

#ifndef FILTERING
#define FILTERING

// Defin PI when using g++
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// Including external libraries
#include <stdlib.h>
#include <cfloat>
#include <iostream>
#include <Eigen/Eigen>

// Including internal libraries
#include "DataHandle.h"

using namespace Eigen;
using namespace DataHandle;

namespace Filtering {

    class KFS_fit
    {
    public:
        Eigen::MatrixXd Lambda_hat; // Matrix of estimated factor loadings
        Eigen::MatrixXd Pt; // Conditional variance of factors at t+1 given observations up to t
        Eigen::MatrixXd F; // Smoothed estimates of latent factors
        Eigen::MatrixXd Wt; // Conditional variance of factors at time t given all observations
        Eigen::MatrixXd C; // Inverse of the Lower Cholesky-Deco of Var-Cov
        int order; // Order of VAR process for factors in the transition equation
        bool conv; // Whether the filter converged

        // Default Constructor
        KFS_fit()
        {
            Lambda_hat = MatrixXd::Zero(0, 0);
            Pt = MatrixXd::Zero(0, 0);
            F = MatrixXd::Zero(0, 0);
            Wt = MatrixXd::Zero(0, 0);
            C = MatrixXd::Zero(0, 0);
            order = 0.;
            conv = 1;
        }

        // Constructor
        KFS_fit(int K, int T, int N) :
            Lambda_hat(N, T),
            Pt(K, K * T),
            F(K, T),
            Wt(K, K * T),
            C(N, N),
            order(0),
            conv(0)
        {}

        // Copy constructor
        KFS_fit(const KFS_fit& other) :
            Lambda_hat(other.Lambda_hat),
            Pt(other.Pt),
            F(other.F),
            Wt(other.Wt),
            C(other.C),
            order(other.order),
            conv(other.conv)
        {}

        // Move constructor
        KFS_fit(KFS_fit&& other) noexcept :
            Lambda_hat(std::move(other.Lambda_hat)),
            Pt(std::move(other.Pt)),
            F(std::move(other.F)),
            Wt(std::move(other.Wt)),
            C(std::move(other.C)),
            order(other.order),
            conv(other.conv)
        {}

        // Copy assignment operator
        KFS_fit& operator=(const KFS_fit& other) {
            if (this != &other) {
                Lambda_hat = other.Lambda_hat;
                Pt = other.Pt;
                F = other.F;
                Wt = other.Wt;
                C = other.C;
                order = other.order;
                conv = other.conv;
            }
            return *this;
        }

        // Move assignment operator
        KFS_fit& operator=(KFS_fit&& other) noexcept {
            if (this != &other) {
                Lambda_hat = std::move(other.Lambda_hat);
                Pt = std::move(other.Pt);
                F = std::move(other.F);
                Wt = std::move(other.Wt);
                C = std::move(other.C);
                order = other.order;
                conv = other.conv;
            }
            return *this;
        }

        // Destructor
        ~KFS_fit()
        {
            // Deallocate any dynamically allocated memory
            Lambda_hat.resize(0, 0);
            Pt.resize(0, 0);
            F.resize(0, 0);
            Wt.resize(0, 0);
        }
    };

    extern

        // Univariate representation of the multivariate Kalman Filter
        void UVMVKalmanFilter(
            Filtering::KFS_fit& results,
            MatrixXd& Vt_inv,
            MatrixXd& Kt,
            MatrixXd& e,
            const int& N,
            const int& T,
            const int& K,
            const MatrixXd& X_in,
            MatrixXd& Lambda,
            const MatrixXd& Phi,
            const MatrixXd& Sigma_e,
            const MatrixXd& Sigma_epsilon,
            const MatrixXd& Pt0,
            const VectorXd& Ft0,
            const VectorXi& missings,
            const VectorXd& tau,
            const VectorXd& alpha,
            const VectorXd& ct,
            const double delta2,
            const double delta1,
            const double sigma,
            const int& h = 0
        );

    // Univariate representation of the multivariate Kalman Smoother
    void UVMVKalmanSmoother(
        Filtering::KFS_fit& results,
        MatrixXd& Kt,
        MatrixXd& Vt_inv,
        MatrixXd& e,
        const int& N,
        const int& T,
        const int& K,
        const MatrixXd& X,
        const MatrixXd& Lambda,
        const MatrixXd& Phi,
        const MatrixXd& Sigma_e,
        const MatrixXd& Sigma_epsilon,
        const VectorXi& missings
    );
};
#endif /* defined(FILTERING) */

