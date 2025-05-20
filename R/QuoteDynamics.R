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
#' @param start Vector of starting values.
#' @param X Matrix of data.
#' @param tau Vector of parameter values.
#' @param ident_mat Selection matrix.
#' @param rel_xtol Relative parameter change tolerance.
#' @param rel_ftol Relative obj. func. tolerance.
#' @param max_eval Maximum number of obj. function evaluations.
#' @param algorithm Character string, the optimization algorithm to use.
#' @param hessian Whether to compute the hessian matrix.
#' @param step_size Step size for hessian computation.
#' @param verbose Whether to print verbose output.
#' @return A list containing the minimum value found, estimated parameters, and optionally the hessian matrix.
#' @export
#' @examples
#' start <- c(1, 1)
#' data <- matrix(rnorm(100), ncol = 10)
#' tau <- runif(10)
#' ident_mat <- matrix(round(runif(100 * 10, -0.5, 1.5), 0), 100, 10)
#' result <- NLoptKFS(start, data, tau, ident_mat)
quoteDynOptim <- function(start, X, tau, ident_mat, rel_xtol = 1e-8, rel_ftol = 1e-8, max_eval = 100000L, algorithm = "LN_NELDERMEAD", hessian = TRUE, step_size = 1e-8, verbose = FALSE) {

  # Start error handling block #

  N <- dim(X)[2]
  M <- N / 2
  T <- dim(X)[1]

  # Convert inputs to matrices #

  errorOccurred <- FALSE

  XR <- try(as.matrix(X), silent = TRUE)
  if (inherits(XR, "try-error")) {
    errorOccurred <- TRUE
  }

  # Logic block for handling different input types for the start vector

  if(is.list(start) && !is.data.frame(start)){ # Logic for converting the start values into a matrix if its provided as a named list

    # Check the names of the parameter list
    if (!setequal(names(start), c("spread","alpha","Cmat.lowerTriangular","sigma","delta1","delta2"))) {
      stop("An error occurred in the conversion of the starting parameters list. The elements of the starting parameter",
           "list must follow this naming convention: \n names(start_list) = c(\"spread\", \"alpha\", \"Cmat.lowerTriangular\",",
           "\"sigma\", \"delta1\", \"delta2\")")
    }

    # Check the dimensions of the parameters stored in the list

    if (length(start$spread) != M ||
        length(start$alpha) != 2 * M ||
        length(start$Cmat.lowerTriangular) != 3 * M ||
        length(start$sigma) != 1 ||
        length(start$delta1) != 1 ||
        length(start$delta2) != 1) {
      stop("An error occurred in the conversion of the starting parameters list. The parameters appear to have false dimensions.",
           "Given the data you provided, the elements of start_list should have the following dimensions: |start_list$spread| = ", M,
           ", |alpha| = ", 2*M, ", |Cmat.lowerTriangular| = ", 3*M, ", |sigma| = ", 1, ", |delta1| = ", 1, ", |delta2| = ", 1)
    }

    startR <- matrix(NaN, length(unlist(start)), 1)
    startR[1:M, 1] <- start$spread
    startR[(M + 1):(3*M), 1] <- start$alpha
    startR[(3*M + 1):(3*M + 9), 1] <- start$Cmat.lowerTriangular
    startR[3*M + 10, ] <- start$sigma
    startR[3*M + 11, ] <- start$delta1
    startR[3*M + 12, ] <- start$delta2

  }else{ # Logic for converting the start values into a matrix if its provided as a vector/matrix

    startR <- if (!errorOccurred) try(as.matrix(start), silent = TRUE)
    if (inherits(startR, "try-error")) {
      errorOccurred <- TRUE
    }

  }

  tauR <- if (!errorOccurred) try(as.matrix(tau), silent = TRUE) else NULL
  if (inherits(tauR, "try-error")) {
    errorOccurred <- TRUE
  }

  ident_mat_R <- try(matrix(as.numeric(ident_mat), dim(ident_mat)[1], dim(ident_mat)[2]), silent = TRUE)
  if (inherits(ident_mat_R, "try-error")) {
    errorOccurred <- TRUE
  }

  # If an error occurred, terminate
  if (errorOccurred) {
    stop("An error occurred in matrix conversion. \"start\", \"X\", \"tau\", and \"ident_mat_R\" must be convertable to matrices. \n")
  }

  # Dimension checks #

  if(dim(tauR)[1] != T){
    stop(paste0("tau must be of length ", T, " but is ", dim(tauR)[1]))
  }

  if(dim(ident_mat_R)[2] != M || dim(ident_mat_R)[1] != T){
    stop(paste0("ident_mat_R must be of dimensions ", T, "x", M, " but is ", dim(ident_mat_R)[1], "x", dim(ident_mat_R)[2]))
  }

  # Check of the elements of the indicator matrix

  if(!all((ident_mat_R - floor(ident_mat_R)) == 0)){
    stop("ident_mat must be a matrix of integer entries only!")
  }

  if(any(ident_mat_R < 0) || any(ident_mat_R > 1)){
    stop("The elements of ident_mat must be either 1 or zero!")
  }

  # Bounds check for the NLopt algorithm coefficients #

  if(!is.null(rel_xtol) && !is.na(rel_xtol) && !is.infinite(rel_xtol)){

    if(rel_xtol < 0){
      stop("rel_xtol must be strictly non-negative!")
    }

  }else{
    stop("rel_xtol cannot be set to NULL/NaN/Inf!")
  }

  if(!is.null(rel_ftol) && !is.na(rel_ftol) && !is.infinite(rel_ftol)){

    if(rel_ftol < 0){
      stop("rel_ftol must be non-negative!")
    }

  }else{
    stop("rel_ftol cannot be set to NULL/NaN/Inf!")
  }

  if(!is.null(max_eval) && !is.na(max_eval) && !is.infinite(max_eval)){

    if(max_eval - floor(max_eval) != 0){
      stop("max_eval must be an integer!")
    }
    if(max_eval <= 0){
      stop("max_eval must be strictly positive!")
    }

  }else{
    cat("Warning! max_eval is disabled.")
    max_eval <- -1
  }

  # Algorithm-type check #

  if(!(algorithm %in% names(nlopt_algorithms))){
    stop(paste0(algorithm, " is not a valid algorithm. See PrintAlgorithms() for a list of the available optimisation algorithms."))
  }
  algorithm_id <- nlopt_algorithms[algorithm]

  # Bounds check for the parameters for the computation of the hessian #

  if(hessian == 1){
    hessian <- TRUE
  }else if(hessian == 0){
    hessian <- FALSE
  }else if(!is.logical(hessian)){
    stop("hessian must be boolean.")
  }

  if(hessian){
    if(!is.null(step_size) && !is.na(step_size) && !is.infinite(step_size)){

      if(step_size <= 0){
        stop("step_size must be strictly positive!")
      }
    }else{
      stop("step_size cannot be set to NULL/NaN/Inf!")
    }
  }

  # Check misc. parameters #

  if(verbose == 1){
    verbose <- TRUE
  }else if(verbose == 0){
    verbose <- FALSE
  }else if(!is.logical(verbose)){
    stop("verbose must be boolean.")
  }

  # End error handling block #

  out <- .Call('_QuoteDynamics_FastOptim', PACKAGE = 'QuoteDynamics', startR, XR, tauR, ident_mat_R, rel_xtol, rel_ftol, max_eval, algorithm_id, hessian, step_size, verbose)

  optimisation_status <- ""
  if(out$nlopt_return_code == 1){
    optimisation_status <- "NLopt exited succesfully."
  }else if(out$nlopt_return_code == 2){
    optimisation_status <- "Objective function stopping value is reached."
  }else if(out$nlopt_return_code == 3){
    optimisation_status <- "Objective function relative/absolute stopping criterion is reached."
  }else if(out$nlopt_return_code == 4){
    optimisation_status <- "Parameter change relative/absolute stopping criterion is reached."
  }else if(out$nlopt_return_code == 5){
    optimisation_status <- "Maximum number of obj. function calls is reached."
  }else if(out$nlopt_return_code == 6){
    optimisation_status <- "Maximum computation time is reached."
  }else if(out$nlopt_return_code == -1){
    optimisation_status <- "Error: NLopt did not exit succesfully."
  }else if(out$nlopt_return_code == -2){
    optimisation_status <- "Error: Invalid arguments. NLopt did not exit succesfully."
  }else if(out$nlopt_return_code == -3){
    optimisation_status <- "Error: Insufficient memory. NLopt did not exit succesfully."
  }else if(out$nlopt_return_code == -4){
    optimisation_status <- "Warning: NLopt exits due to roundoff errors limited progress. NLopt may not have converged."
  }else if(out$nlopt_return_code == -5){
    optimisation_status <- "Warning: Userforced exit. NLopt did not exit succesfully."
  }

  out_sthree <- structure(list(
    call = match.call(),
    minf = out$min_val,
    estimate = out$estimate,
    hessian = out$hessian,
    eval_count = out$eval_count,
    nlopt_code = out$nlopt_return_code,
    optimisation_status = optimisation_status
  ), class = "quoteDynFit")

  return(out_sthree)

}

#'  Summary function for the return object of the optimisation call
#'
#' @export
summary.quoteDynFit <- function(object, ...) {
  cat("------ QuoteDynamics Estimation Summary ------\n\n")
  cat("Call: ")
  print(object$call)
  cat("\n")
  cat("##############################################\n")
  cat("NLopt Return Code:", object$optimisation_status, "\n")
  cat("Minimum Value of Obj.Func.:", object$min_val, "\n")
  cat("Number of Evaluations:", object$eval_count, "\n")
  cat("Estimated Parameters:\n")
  print(object$estimate)
  cat("\n")
  if(object$call$hessian){
    cat("Hessian:\n")
    print(object$hessian)
  }
  cat("\n")
  cat("##############################################\n")
  invisible(object)
}


#'  Compute the negative log-likelihood of the quote dynamics model
#'
#' @param start Vector of starting values.
#' @param X Matrix of data.
#' @param tau Vector of parameter values.
#' @param ident_mat Selection matrix.
#' @param verbose Whether to print verbose output.
#' @return Negative log-likelihood
#' @export
#' @examples
#' start <- c(1, 1)
#' data <- matrix(rnorm(100), ncol = 10)
#' tau <- runif(10)
#' ident_mat <- matrix(round(runif(100 * 10, -0.5, 1.5), 0), 100, 10)
#' quoteDynNegLL(start, data, tau, ident_mat)
quoteDynNegLL <- function(start, X, tau, ident_mat, verbose = FALSE) {

  # Start error handling block #

  N <- dim(X)[2]
  M <- N / 2
  T <- dim(X)[1]

  # Convert inputs to matrices #

  errorOccurred <- FALSE

  XR <- try(as.matrix(X), silent = TRUE)
  if (inherits(XR, "try-error")) {
    errorOccurred <- TRUE
  }

  # Logic block for handling different input types for the start vector

  if(is.list(start) && !is.data.frame(start)){ # Logic for converting the start values into a matrix if its provided as a named list

    # Check the names of the parameter list
    if (!setequal(names(start), c("spread","alpha","Cmat.lowerTriangular","sigma","delta1","delta2"))) {
      stop("An error occurred in the conversion of the starting parameters list. The elements of the starting parameter",
           "list must follow this naming convention: \n names(start_list) = c(\"spread\", \"alpha\", \"Cmat.lowerTriangular\",",
           "\"sigma\", \"delta1\", \"delta2\")")
    }

    # Check the dimensions of the parameters stored in the list

    if (length(start$spread) != M ||
        length(start$alpha) != 2 * M ||
        length(start$Cmat.lowerTriangular) != 3 * M ||
        length(start$sigma) != 1 ||
        length(start$delta1) != 1 ||
        length(start$delta2) != 1) {
      stop("An error occurred in the conversion of the starting parameters list. The parameters appear to have false dimensions.",
           "Given the data you provided, the elements of start_list should have the following dimensions: |start_list$spread| = ", M,
           ", |alpha| = ", 2*M, ", |Cmat.lowerTriangular| = ", 3*M, ", |sigma| = ", 1, ", |delta1| = ", 1, ", |delta2| = ", 1)
    }

    startR <- matrix(NaN, length(unlist(start)), 1)
    startR[1:M, 1] <- start$spread
    startR[(M + 1):(3*M), 1] <- start$alpha
    startR[(3*M + 1):(3*M + 9), 1] <- start$Cmat.lowerTriangular
    startR[3*M + 10, ] <- start$sigma
    startR[3*M + 11, ] <- start$delta1
    startR[3*M + 12, ] <- start$delta2

  }else{ # Logic for converting the start values into a matrix if its provided as a vector/matrix

    startR <- if (!errorOccurred) try(as.matrix(start), silent = TRUE)
    if (inherits(startR, "try-error")) {
      errorOccurred <- TRUE
    }

  }

  tauR <- if (!errorOccurred) try(as.matrix(tau), silent = TRUE) else NULL
  if (inherits(tauR, "try-error")) {
    errorOccurred <- TRUE
  }

  ident_mat_R <- try(matrix(as.numeric(ident_mat), dim(ident_mat)[1], dim(ident_mat)[2]), silent = TRUE)
  if (inherits(ident_mat_R, "try-error")) {
    errorOccurred <- TRUE
  }

  # If an error occurred, terminate
  if (errorOccurred) {
    stop("An error occurred in matrix conversion. \"start\", \"X\", \"tau\", and \"ident_mat_R\" must be convertable to matrices. \n")
  }

  # Dimension checks #

  if(dim(tauR)[1] != T){
    stop(paste0("tau must be of length ", T, " but is ", dim(tauR)[1]))
  }

  if(dim(ident_mat_R)[2] != M || dim(ident_mat_R)[1] != T){
    stop(paste0("ident_mat_R must be of dimensions ", T, "x", M, " but is ", dim(ident_mat_R)[1], "x", dim(ident_mat_R)[2]))
  }

  # Check of the elements of the indicator matrix

  if(!all((ident_mat_R - floor(ident_mat_R)) == 0)){
    stop("ident_mat must be a matrix of integer entries only!")
  }

  if(any(ident_mat_R < 0) || any(ident_mat_R > 1)){
    stop("The elements of ident_mat must be either 1 or zero!")
  }

  # End error handling block #

  .Call('_QuoteDynamics_objFunctionCpp', PACKAGE = 'QuoteDynamics', startR, XR, tauR, ident_mat_R, verbose)

}
