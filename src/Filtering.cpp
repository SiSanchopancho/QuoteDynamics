#include "Internals/Filtering.h"

void Filtering::UVMVKalmanFilter(
	Filtering::KFS_fit& results, // Object for storing the results
	MatrixXd& Vt_inv, // Inverse of the variance matrix V
	MatrixXd& Kt, // Matrix storing the Kalman gain at time t
	MatrixXd& e, //Matrix storing the errors
	const int& N, // Number of variables
	const int& T, // Number of observations
	const int& K, // Number of factors
	const MatrixXd& X_in, // Data matrix of dimensions (NxT) !
	MatrixXd& Lambda, // Factor Loading Matrix
	const MatrixXd& Phi, // VAR coefficient matrix for factor process
	const MatrixXd& Sigma_e, // Idiosyncratic error variance-covariance-matrix
	const MatrixXd& Sigma_epsilon, // Variance-covariance matrix of the factor VAR process
	const MatrixXd& Pt0, // Initial conditional variance-covariance matrix of the factors 
	const VectorXd& Ft0, // Initial value of the factors
	const VectorXi& missings, // Vector indicating the obseravtion index for missing observations for each variable x_n in X
	const VectorXd& tau,
	const VectorXd& alpha,
	const VectorXd& ct,
	const double delta2,
	const double delta1,
	const double sigma,
	const int& h // Forecasting horizon
)
{
	// Kalman Filter

	// Dummies

	// Reals
	double Vt = 0.;

	// Vectors
	VectorXd Ftt = VectorXd::Zero(K);

	// Matrices
	MatrixXd Ptt = MatrixXd::Zero(K, K), Lambda_T = Lambda.transpose(), X = X_in, R = MatrixXd::Zero(K, K);

	// Initialisation

	results.F.setZero();
	results.Pt.setZero();

	results.F.col(0) = Ft0;

	results.Pt.block(0, 0, K, K) = Pt0;

	for (int t = 0; t < T; ++t)
	{

		Lambda.col(0).setConstant(1);
		Lambda.col(1) = alpha * (std::pow(tau(t), delta2) * sigma);

		Lambda_T = Lambda.transpose();

		Ftt = results.F.col(t);
		Ptt = results.Pt.block(0, (t)*K, K, K);

		// Only for the log-like save Vt by abusing the space saved for Wt in the smoothing step which is not required

		for (int n = 0; n < N; ++n)
		{
			if (T - missings(n) <= t)
			{
				// Skip missing observations

				continue;
			}

			// Calculate the one-step ahead filter estimates of the factors and their conditional variances

			Vt = Lambda.row(n) * Ptt * Lambda_T.col(n) + Sigma_e.diagonal()(n);
			results.Wt.block(0, t * N, N, N)(n, n) = Vt;
			e(n, t) = X(n, t) - ct(n) - Lambda.row(n) * Ftt;
			Vt_inv(n, t) = 1. / Vt;
			Kt.col(t * N + n) = Ptt * Lambda_T.col(n);
			Ftt += (Vt_inv(n, t)) * Kt.col(t * N + n) * e(n, t);
			Ptt -= Kt.col(t * N + n) * (Vt_inv(n, t)) * Kt.col(t * N + n).transpose();

		}

		// Store the results

		results.F.col(t) = Ftt;
		results.Pt.block(0, (t)*K, K, K) = Ptt;

		// Predict t + 1

		if (t < T - 1) {

			int ind = 0;
			for (int d1 = 0; d1 < K; ++d1) {

				for (int d2 = 0; d2 < K; ++d2) {

					if (ind % 4 == 0) { R(d1, d2) = 0; }
					else
						if (ind % 4 == 1) { R(d1, d2) = std::pow(tau(t + 1), delta1) * sigma; }
						else
							if (ind % 4 == 2) { R(d1, d2) = 0; }
							else { R(d1, d2) = 1; }

					++ind;

				}
			}

			// Update

			results.F.col(t + 1) = Phi * Ftt;
			results.Pt.block(0, (t + 1) * K, K, K) = Phi * Ptt * Phi.transpose() + R * Sigma_epsilon * R.transpose();
		}

	}

	if (0 < h)
	{
		for (int t = T; t < T + h; ++t)
		{
			results.F.col(t + 1) = Phi * results.F.col(t);
			results.Pt.block(0, (t + 1) * K, K, K) = Phi * results.Pt.block(0, t * K, K, K) * Phi.transpose() + Sigma_epsilon;
		}
	}

	return;
}

void Filtering::UVMVKalmanSmoother(
	Filtering::KFS_fit& results, // Object for storing the results
	MatrixXd& Kt, // Matrix storing the Kalman gain at time t
	MatrixXd& Vt_inv, // Inverse of the variance matrix V
	MatrixXd& e, //Matrix storing the errors
	const int& N, // Number of variables
	const int& T, // Number of observations
	const int& K, // Number of factors
	const MatrixXd& X, // Data matrix of dimensions (NxT) !
	const MatrixXd& Lambda, // Factor Loading Matrix
	const MatrixXd& Phi, // VAR coefficient matrix for factor process
	const MatrixXd& Sigma_e, // Idiosyncratic error variance-covariance-matrix
	const MatrixXd& Sigma_epsilon, // Variance-covariance matrix of the factor VAR process
	const VectorXi& missings
)
{
	// Kalman Smoother

	// Dummies

	// Vectors
	VectorXd rtt = VectorXd::Zero(K);

	// Matrices
	MatrixXd Lt = MatrixXd::Zero(K, K), Ntt = MatrixXd::Zero(K, K), Lambda_T = Lambda.transpose(), IdentK = MatrixXd::Identity(K, K);

	for (int t = T - 1; 0 <= t; --t)
	{

		for (int n = N - 1; 0 <= n; --n)
		{
			if (T - missings(n) <= t)
			{
				// Skip missing observations

				continue;

			}

			Lt = IdentK - Kt.col(t * N + n) * Lambda.row(n) * Vt_inv(n, t);
			rtt = Lambda_T.col(n) * Vt_inv(n, t) * e(n, t) + Lt.transpose() * rtt;
			Ntt = Lambda_T.col(n) * Vt_inv(n, t) * Lambda.row(n) + Lt.transpose() * Ntt * Lt;

		}

		// Calculate the smoothed estimates of the factors and their variance matrix and store the results

		results.F.col(t) += results.Pt.block(0, t * K, K, K) * rtt;
		results.Wt.block(0, t * K, K, K) = results.Pt.block(0, t * K, K, K) - results.Pt.block(0, t * K, K, K) * Ntt * results.Pt.block(0, t * K, K, K);
		rtt = Phi.transpose() * rtt;
		Ntt = Phi.transpose() * Ntt * Phi;

	}

	return;
}

