# QuoteDynamics Package

## Overview
**QuoteDynamics** is an R package that provides **fast maximum‑likelihood estimation for high‑frequency bid/ask quote models**. The core is a C++ implementation of a multivariate Kalman filter combined with the NLopt optimisation library. Heavy numerical work is done in C++17 with Eigen for linear algebra and NLopt for (derivative‑free) optimisation. Rcpp and RcppEigen handle the interface to R.

### Key Features

- Block‑diagonal state‑space structure exploited for speed
- Choice of 13 NLopt algorithms (PrintAlgorithms() lists them)
- Optional numerical Hessian (central differences)
- Fully open source – GPL ≥ 3, links LGPL NLopt dynamically
- Works on Windows, macOS, and Linux

## Installation

The library is not on CRAN yet and must be installed from GitHub.

### System requirements

| Platform    | NLopt note                                                                   |
| ----------- | ---------------------------------------------------------------------------- |
| **Windows** | Pre‑compiled `nlopt.dll` is shipped inside the package.      		     |
| **Linux**   | Install dev package beforehand, e.g. <br>`sudo apt-get install libnlopt-dev` |
| **macOS**   | Install via Homebrew: <br>`brew install nlopt`                               |

No manual Eigen/RcppEigen steps are required as CRAN’s RcppEigen is used.

## Quick start

Generate dummy data with two markets (`M = 2`, i.e., four series) and fit the model:

```R
library(QuoteDynamics)

set.seed(42)
T  <- 100L           # time points
M  <- 2L             # markets
N  <- 2 * M          # bid + ask per market

X   <- matrix(rnorm(T * N), nrow = T, ncol = N)
tau <- rexp(T, rate = 5)                  # durations
mask <- matrix(1, nrow = T, ncol = M)     # all quotes active

# naive starting values
theta0 <- c(
  spread  = rep(0.02, M),
  alpha   = rep( 0.9, 2 * M),
  Clower  = rep( 0.01, 3 * M),
  sigma   = 0.1,
  delta1  = 0.3,
  delta2  = 0.7
)

## run optimiser
fit <- quoteDynOptim(
  start      = theta0,
  X          = X,
  tau        = tau,
  ident_mat  = mask,
  algorithm  = "LN_NELDERMEAD",
  xtol       = 1e-6,
  stop_val   = 1e-12,
  verbose    = TRUE,
  hessian    = TRUE
)

str(fit, max.level = 1)

```

For a list of available NLopt algorithms together with their integer IDs, call:

```R
PrintAlgorithms()
```

## Function reference

| Function            | Purpose                                                  |
| ------------------- | -------------------------------------------------------- |
| `quoteDynOptim()`   | High‑level R wrapper with argument checks                |
| `FastOptim()`       | Low‑level C++ routine            			 |
| `PrintAlgorithms()` | Prints NLopt algorithm table                             |
| `nlopt_algorithms`  | Named vector: `"ALG_NAME" = id`                          |

## License

- Package code: GPL‑3 (or later)
- NLopt (embedded/shared): LGPL 2.1+
- Eigen (via RcppEigen): MPL 2.0

See the `LICENSE` file and the original license texts in `inst/licenses/`.

## References

To acknowledge the core libraries that QuoteDynamics builds on, please cite the following when appropriate:

| Library       | Reference                                                                                                                                                            |
| ------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Eigen**     | Guennebaud, G., Jacob, B. *et al.* (2024). **Eigen 3.4** – A C++ template library for linear algebra. <[http://eigen.tuxfamily.org\\>](http://eigen.tuxfamily.org\>) |
| **NLopt**     | Johnson, S.G. (2014). **The NLopt nonlinear‑optimization package**. <[https://nlopt.readthedocs.io\\>](https://nlopt.readthedocs.io\>)                               |
| **Rcpp**      | Eddelbuettel, D., François, R. (2011). *Rcpp: Seamless R and C++ Integration.* *Journal of Statistical Software*, **40**(8), 1‑18.                                   |
| **RcppEigen** | Bates, D., Eddelbuettel, D. (2013). *Fast and Elegant Numerical Linear Algebra Using the RcppEigen Package.* *Journal of Statistical Software*, **52**(5), 1‑24.     |


## Citation

Coming soon

## Contact

For any queries or feedback, please contact me.
