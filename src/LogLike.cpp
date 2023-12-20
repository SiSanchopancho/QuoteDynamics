#include "Internals/LogLike.h"

double LogLike::LL(
	VectorXd param,
	MatrixXd X,
	VectorXd tau,
	double m,
	unsigned d
) {

	VectorXd spread = param(seq(0, 1));
	VectorXd ct(4);
	ct(0) = -spread(0) / 2.;
	ct(1) = spread(0) / 2.;
	ct(2) = -spread(1) / 2.;
	ct(3) = spread(1) / 2.;
	VectorXd alpha = param(seq(2, 5));
	MatrixXd Cmat = MatrixXd::Zero(2, 2);

	int ind = 0;

	for (int i = 0; i < Cmat.rows(); i++) {

		for (int j = 0; j <= i; j++) {

			Cmat(i, j) = param(6 + ind);
			++ind;

		}
	}

	MatrixXd Omega1 = Cmat * Cmat.transpose();

	Cmat.setZero();

	ind = 0;

	for (int i = 0; i < Cmat.rows(); i++) {

		for (int j = 0; j <= i; j++) {

			Cmat(i, j) = param(9 + ind);
			++ind;

		}
	}

	MatrixXd Omega2 = Cmat * Cmat.transpose();

	MatrixXd GGt = MatrixXd::Zero(2 * Cmat.rows(), 2 * Cmat.rows());
	GGt.topLeftCorner(Cmat.rows(), Cmat.rows()) = Omega1;
	GGt.bottomRightCorner(Cmat.rows(), Cmat.rows()) = Omega2;

	VectorXd dt = VectorXd::Zero(2);
	MatrixXd Tt = MatrixXd::Zero(d, d);

	ind = 0;

	for (int d1 = 0; d1 < d; ++d1) {

		for (int d2 = 0; d2 < d; ++d2) {

			Tt(d2, d1) = ind % 4 == 0 ? 1. : 0.;
			++ind;

		}
	}

	double sigma = param(12);
	double delta1 = param(13);
	double delta2 = param(14);

	// Dummies

	// Integers

	int T = X.rows(), N = X.cols();

	// Matrices

	MatrixXd Vt_inv = MatrixXd::Zero(N, T), Kt = MatrixXd::Zero(d, N * T), e = MatrixXd::Zero(N, T), Pt0 = MatrixXd::Zero(d, d), Ft0 = MatrixXd::Zero(d, 1), Lambda = MatrixXd::Zero(N, d),
		Phi = Tt, Sigma_e = GGt, Sigma_epsilon = MatrixXd::Identity(d, d), R = MatrixXd::Zero(d, d * T), IN = MatrixXd::Identity(N, N);

	VectorXd lnl = VectorXd::Zero(T);

	// Set Pt0

	Pt0.diagonal().setConstant(1000);

	// Filtering

	// Filter

	KFS_fit fit(d, (T + 1), N);

	fit.Wt.conservativeResize(N, N * (T + 1));
	fit.Wt.setZero();
	
	UVMVKalmanFilter(fit, Vt_inv, Kt, e, N, T, d, X.transpose(), Lambda, Phi, Sigma_e, Sigma_epsilon, Pt0, Ft0, VectorXi::Zero(N), tau, alpha, ct, delta2, delta1, sigma);

	// calculate likelihood for each t



	if (N < 5) {

		for (int t = 0; t < T; ++t) {

			// Use exactr determinant and inverse for small matrices

			lnl(t) = -0.5 * N * std::log(2 * M_PI) - 0.5 * std::log(fit.Wt.block(0, t * N, N, N).determinant()) - 0.5 * (e.col(t).transpose() * fit.Wt.block(0, t * N, N, N).inverse() * e.col(t))(0, 0);

		}
	}
	else
	{
		for (int t = 0; t < T; ++t) {

			// Use approximate determinant and inverse for small matrices

			LLT<MatrixXd> llt(fit.Wt.block(0, t * N, N, N));

			MatrixXd C = llt.matrixL();

			lnl(t) = -0.5 * N * std::log(2 * M_PI) - 0.5 * std::log(C.diagonal().prod() * C.diagonal().prod()) - 0.5 * (e.col(t).transpose() * llt.solve(IN) * e.col(t))(0, 0);
		
		}
	}

	return -lnl.mean();
}