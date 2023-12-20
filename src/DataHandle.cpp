#include "Internals/DataHandle.h"

MatrixXd DataHandle::Cov(const MatrixXd& X_in)
{
    // Dummies 

    // Matices

    MatrixXd X = X_in;

    // Demeaning the data

    VectorXd mu = X.colwise().mean();

    for (int n = 0; n < X.cols(); ++n)
    {
        for (int t = 0; t < X.rows(); ++t)
        {
            X(t, n) -= mu(n, 0);
        }
    }

    // Claculate the empirical variance-covariance matrix

    MatrixXd cov = (X.transpose() * X) / double(X.rows() - 1);

    return cov;
}

MatrixXd DataHandle::Corr(const MatrixXd& X_in)
{

    // Dummies

    //Matrices

    MatrixXd X = X_in;

    // Demeaning the data

    VectorXd mu = X.colwise().mean();

    for (int n = 0; n < X.cols(); ++n)
    {
        for (int t = 0; t < X.rows(); ++t)
        {
            X(t, n) -= mu(n, 0);
        }
    }

    // Claculate the empirical variance-covariance matrix

    MatrixXd cov = (X.transpose() * X) / double(X.rows() - 1);

    // Claculate the empirical correlation matrix using Cor = diag(Cov)^-1/2 * Cov * diag(Cov)^-1/2

    VectorXd ss = cov.diagonal();
    MatrixXd S_sqrt_inv = MatrixXd::Zero(X.cols(), X.cols());

    for (int n = 0; n < X.cols(); ++n)
    {
        S_sqrt_inv(n, n) = 1 / sqrt(ss(n));
    }

    return S_sqrt_inv * cov * S_sqrt_inv;
}

double DataHandle::Corr(const VectorXd& x, const VectorXd& y)
{

    // Dummies

    // Vectors

    VectorXd x_in = x;
    VectorXd y_in = y;

    // Demeaning the data

    VectorXd I = VectorXd::Ones(x.size());
    Demean(x_in);
    Demean(y_in);

    // Claculate the empirical correlation

    double corr = (x_in.dot(y_in)) / (x_in.norm() * y_in.norm());

    return corr;
}

double DataHandle::Corr(const Block<MatrixXd, -1, 1, true>& x, const Block<MatrixXd, -1, 1, true>& y)
{

    // Dummies

    // Vectors

    VectorXd x_in = x;
    VectorXd y_in = y;

    // Demeaning the data

    VectorXd I = VectorXd::Ones(x.size());
    Demean(x_in);
    Demean(y_in);

    // Claculate the empirical correlation

    double corr = (x_in.dot(y_in)) / (x_in.norm() * y_in.norm());

    return corr;
}

void DataHandle::Demean(MatrixXd& X)
{

    // Demean the data

    VectorXd mu = X.colwise().mean();

    for (int n = 0; n < X.cols(); ++n)
    {
        for (int t = 0; t < X.rows(); ++t)
        {
            X(t, n) -= mu(n, 0);
        }
    }
    return;
}

void DataHandle::Demean(VectorXd& v)
{
    v -= v.col(0).mean() * VectorXd::Ones(v.size());

    return;
}

void DataHandle::Standardise(MatrixXd& X)
{

    // standardise the data by demeaning and then normalising

    VectorXd mu = X.colwise().mean();
    VectorXd sigma = Cov(X).diagonal();

    for (int n = 0; n < X.cols(); ++n)
    {
        for (int t = 0; t < X.rows(); ++t)
        {
            X(t, n) -= mu(n);
            X(t, n) /= std::sqrt(sigma(n));
        }
    }

    return;
}

void DataHandle::Standardise(VectorXd& v)
{

    v -= v.col(0).mean() * VectorXd::Ones(v.size());

    double s = std::sqrt(Var(v));

    for (int n = 0; n < v.size(); ++n)
    {
        v(n) /= s;
    }

    return;
}

void DataHandle::Standardise(Block<MatrixXd, -1, -1, false>& v_in)
{
    VectorXd v = v_in;
    v -= v.col(0).mean() * VectorXd::Ones(v.size());

    double s = std::sqrt(Var(v));

    for (int n = 0; n < v.size(); ++n)
    {
        v(n) /= s;
    }

    v_in = v;

    return;
}

double DataHandle::Var(const VectorXd& v_in)
{
    VectorXd v = v_in;
    VectorXd mu_vec = VectorXd::Ones(v.size());
    mu_vec.setConstant(v.mean());
    return (1. / double(v.size())) * ((v - mu_vec).transpose() * (v - mu_vec))(0);
}



MatrixXd DataHandle::KroneckerProd(const MatrixXd& X, const MatrixXd& Y)
{
    const int M = X.rows();
    const int N = X.cols();
    const int P = Y.rows();
    const int R = Y.cols();

    MatrixXd K(M * P, N * R);

    for (int n = 0; n < N; ++n)
    {
        for (int m = 0; m < M; ++m)
        {
            K.block(n * P, m * R, P, R) = X(m, n) * Y;
        }
    }
    return K;
}

void DataHandle::removeRow(MatrixXd& X, const int& t, const bool& conservative)
{
    int T = X.rows() - 1;
    int N = X.cols();

    if (t < T)
    {
        X.block(t, 0, T - t, N) = X.block(t + 1, 0, T - t, N);
    }

    if (conservative)
    {
        X.row(T) = VectorXd::Zero(N);
    }
    else
    {
        X.conservativeResize(T, N);
    }

}

void DataHandle::removeCol(MatrixXd& X, const int& n, const bool& conservative)
{
    int T = X.rows();
    int N = X.cols() - 1;

    if (n < N)
    {
        X.block(0, n, T, N - n) = X.block(0, n + 1, T, N - n);
    }

    if (conservative)
    {
        X.col(N) = VectorXd::Zero(T);
    }
    else
    {
        X.conservativeResize(T, N);
    }
}

void DataHandle::removeElement(VectorXd& x, const int& t, const bool& conservative)
{
    int T = x.size() - 1;

    if (t < T)
    {
        x.block(t, 0, T - t, 1) = x.block(t + 1, 0, T - t, 1);
    }

    if (conservative)
    {
        x(T - 1) = 0.;
    }
    else
    {
        x.conservativeResize(T, 1);
    }
}
