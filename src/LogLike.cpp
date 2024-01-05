#include "Internals/LogLike.h"

double LogLike::LL(
	VectorXd param,
	MatrixXd X,
	VectorXd tau
) {

	// Set the parameters

	const int m = X.cols() / 2;

	VectorXd spread = param.head(m);

	VectorXd ct(2 * m);

	ct(seq(0, 2 * m - 1, 2)) = -0.5 * spread;
	ct(seq(1, 2 * m - 1, 2)) = 0.5 * spread;

	VectorXd alpha = param(seq(m, (3 * m) - 1));

	MatrixXd Cmat = MatrixXd::Zero(2, 2);

	MatrixXd GGt = MatrixXd::Zero(2, 2);

	int ind = 0;

	for (int i = 0; i < Cmat.rows(); i++) {

		for (int j = 0; j <= i; j++) {

			Cmat(i, j) = param(3 * m + ind);
			++ind;

		}
	}

	GGt.bottomRightCorner(Cmat.rows(), Cmat.cols()) = Cmat * Cmat.transpose();


	for (int mm = 1; mm < m; ++mm) {

		Cmat.setZero();

		for (int i = 0; i < Cmat.rows(); i++) {

			for (int j = 0; j <= i; j++) {

				Cmat(i, j) = param(3 * m + ind);
				++ind;

			}
		}

		GGt.conservativeResize(GGt.rows() + Cmat.rows(), GGt.cols() + Cmat.cols());
		GGt.topRightCorner(Cmat.rows() * mm, Cmat.cols()).setZero();
		GGt.bottomLeftCorner(Cmat.rows(), Cmat.cols() * mm).setZero();
		GGt.bottomRightCorner(Cmat.rows(), Cmat.cols()) = Cmat * Cmat.transpose();

	}

	VectorXd dt = VectorXd::Zero(2);
	MatrixXd Tt = MatrixXd::Zero(2, 2);

	Tt(0, 0) = 1.;

	double sigma = param(last - 2);
	double delta1 = param(last - 1);
	double delta2 = param(last);

	// Dummies

	// Integers

	int T = X.rows(), N = X.cols(), d = 2;

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