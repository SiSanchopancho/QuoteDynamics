#' Available Optimization Algorithms for FastOptim
#'
#' A named vector of available optimization algorithms for use with the `FastOptim` function.
#' @name nlopt_algorithms
#' @export
nlopt_algorithms <- c(
  "GN_DIRECT" = 0,
  "GN_DIRECT_L" = 1,
  "GN_ORIG_DIRECT" = 6,
  "GN_CRS2_LM" = 20,
  "GN_MLSL" = 21,
  "GN_ESCH" = 37,
  "LN_COBYLA" = 26,
  "LN_NEWUOA_BOUND" = 28,
  "LN_NELDERMEAD" = 29,
  "LN_SBPLX" = 30,
  "LN_AUGLAG" = 31,
  "LN_BOBYQA" = 34,
  "LN_PRAXIS" = 24
)

#' Print Available NLopt Algorithms and Descriptions
#'
#' Lists the names of available optimization algorithms in `nlopt_algorithms`
#' along with a short description of each.
#'
#' @export
#' @examples
#' PrintAlgorithms()
PrintAlgorithms <- function() {
  cat("Available NLopt Algorithms:\n")

  descriptions <- list(
    "GN_DIRECT" = "Global optimization via the DIviding RECTangles algorithm.",
    "GN_DIRECT_L" = "A variant of GN_DIRECT with potentially faster convergence.",
    "GN_ORIG_DIRECT" = "The original version of the DIRECT algorithm.",
    "GN_CRS2_LM" = "Controlled Random Search algorithm with local mutation.",
    "GN_MLSL" = "Multi-Level Single-Linkage, global optimization using a sequence of local searches.",
    "GN_ESCH" = "Global optimization using the ESCH evolutionary algorithm.",
    "LN_COBYLA" = "Local optimization using the COBYLA (Constrained Optimization BY Linear Approximations) algorithm.",
    "LN_NEWUOA_BOUND" = "Local derivative-free optimization using the NEWUOA algorithm, bound constrained.",
    "LN_NELDERMEAD" = "Local optimization using the Nelder-Mead simplex algorithm.",
    "LN_SBPLX" = "Local optimization using a variant of the Nelder-Mead simplex algorithm.",
    "LN_AUGLAG" = "Local optimization using the AUGmented LAGrangian algorithm.",
    "LN_BOBYQA" = "Local optimization using the BOBYQA (Bound Optimization BY Quadratic Approximation) algorithm.",
    "LN_PRAXIS" = "Local optimization using the PRAXIS (PRincipal AXIS) algorithm."
    # Add other algorithm descriptions here
  )

  for (alg in names(nlopt_algorithms)) {
    desc <- descriptions[[alg]]
    cat(sprintf(" - %s (%s): %s\n", alg, nlopt_algorithms[alg], desc))
  }
}


#' Fast Kalman Filter and Smoother Maximum Likelihood Estimation
#'
#' This function provides a fast implementation of the Kalman Filter and Smoother
#' Maximum Likelihood Estimation using the NLopt optimization library.
#'
#' @param start Numeric vector of starting values.
#' @param X Numeric matrix of data.
#' @param tau Numeric vector of parameter values.
#' @param xtol Numeric, the algorithm tolerance.
#' @param stop_val Numeric, the stopping value.
#' @param algorithm Character string, the optimization algorithm to use (default: "LN_NELDERMEAD").
#' @param hessian Logical, whether to compute the hessian matrix.
#' @param step_size Numeric, step size for hessian computation.
#' @param verbose Logical, whether to print verbose output.
#' @return A list containing the minimum value found, estimated parameters, and optionally the hessian matrix.
#' @export
#' @examples
#' start <- c(1, 1)
#' data <- matrix(rnorm(100), ncol = 10)
#' tau <- runif(10)
#' result <- NLoptKFS(start, data, tau, 2, 2)
quoteDynOptim <- function(start, X, tau, xtol = 10e-5, stop_val = 10e-10, algorithm = "LN_NELDERMEAD", hessian = TRUE, step_size = 10e-5, verbose = FALSE) {

  # Convert inputs to matrices

  errorOccurred <- FALSE

  ## If start is a list of named parameters, convert start to a vector

  XR <- try(as.matrix(X), silent = TRUE)
  if (inherits(XR, "try-error")) {
    errorOccurred <- TRUE
  }

  if(is.list(start) && !is.data.frame(start)){

    # Check the names of the parameter list

    if (sum(names(start) == c("spread", "alpha", "Cmat.lowerTriangular", "sigma", "delta1", "delta2")) == 6) {
      cat("An error occurred in the conversion of the starting parameters list. The elements of the starting parameter",
          "list must follow this naming convention: \n names(start_list) = c(\"spread\", \"alpha\", \"Cmat.lowerTriangular\",",
          "\"sigma\", \"delta1\", \"delta2\")")
      return()
    }

    # Check the dimensions of the parameters stored in the list

    M <- dim(X)[2] / 2

    if (length(start$spread) != M ||
        length(start$alpha) != 2 * M ||
        length(start$Cmat.lowerTriangular) == 3 * M ||
        length(start$sigma) != 1 ||
        length(start$delta1) != 1 ||
        length(start$delta1) != 1) {
      cat("An error occurred in the conversion of the starting parameters list. The parameters appear to have false dimensions.",
          "Given the data you provided, the elements of start_list should have the following dimensions: |start_list$spread| = ", M,
          ", |alpha| = ", 2*M, ", |Cmat.lowerTriangular| = ", 3*M, ", |sigma| = ", 1, ", |delta1| = ", 1, ", |delta2| = ", 1)
      return()
    }

    startR <- matrix(NaN, length(unlist(start)), 1)
    startR[1:M, 1] <- start$spread
    startR[(M + 1):(3*M), 1] <- start$alpha
    startR[(3*M + 1):(3*M + 9), 1] <- start$Cmat.lowerTriangular
    startR[3*M + 10, ] <- start$sigma
    startR[3*M + 11, ] <- start$delta1
    startR[3*M + 12, ] <- start$delta2
  }else{
    startR <- if (!errorOccurred) try(as.matrix(start), silent = TRUE)
    if (inherits(startR, "try-error")) {
      errorOccurred <- TRUE
    }
  }

  tauR <- if (!errorOccurred) try(as.matrix(tau), silent = TRUE) else NULL
  if (inherits(tauR, "try-error")) {
    errorOccurred <- TRUE
  }

  # If an error occurred, terminate
  if (errorOccurred) {
    cat("An error occurred in matrix conversion. \"start\", \"X\", and \"tau\" must be convertable to matrices. \n")
    return()
  }

  if(!is.matrix(XR) || !is.matrix(startR) || !is.matrix(tauR)){

    print("start, X, and tau must be convertible to matrices!")

    return()

  }

  algorithm_id <- nlopt_algorithms[algorithm]
  .Call('_KFSMLE_FastOptim', PACKAGE = 'KFSMLE', startR, XR, tauR, xtol, stop_val, algorithm_id, hessian, step_size, verbose)
}
