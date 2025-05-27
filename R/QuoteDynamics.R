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

checkParameterVector <- function(param, M){

  if(is.list(param) && !is.data.frame(param)){ # Logic for converting the param values into a matrix if its provided as a named list

    # Check the names of the parameter list
    if (!setequal(names(param), c("spread","alpha","Cmat.lowerTriangular","sigma","delta1","delta2"))) {
      stop("An error occurred in the conversion of the ", substitute(param), " parameters list. The elements of the ", substitute(param), " parameter",
           "list must follow this naming convention: \n names(start_list) = c(\"spread\", \"alpha\", \"Cmat.lowerTriangular\",",
           "\"sigma\", \"delta1\", \"delta2\")")
    }

    # Check the dimensions of the parameters stored in the list

    if (length(param$spread) != M ||
        length(param$alpha) != 2 * M ||
        length(param$Cmat.lowerTriangular) != 3 * M ||
        length(param$sigma) != 1 ||
        length(param$delta1) != 1 ||
        length(param$delta2) != 1) {
      stop("An error occurred in the conversion of the ", substitute(param), " parameters list. The parameters appear to have false dimensions.",
           "Given the data provided, the elements of ", substitute(param), " should have the following dimensions: |spread| = ", M,
           ", |alpha| = ", 2*M, ", |Cmat.lowerTriangular| = ", 3*M, ", |sigma| = ", 1, ", |delta1| = ", 1, ", |delta2| = ", 1)
    }

    param_out <- matrix(NaN, length(unlist(param)), 1)
    param_out[1:M, 1] <- param$spread
    param_out[(M + 1):(3*M), 1] <- param$alpha
    param_out[(3*M + 1):(6*M), 1] <- param$Cmat.lowerTriangular
    param_out[6*M + 1, ] <- param$sigma
    param_out[6*M + 2, ] <- param$delta1
    param_out[6*M + 3, ] <- param$delta2

  }else{ # Logic for converting the param values into a matrix if its provided as a vector/matrix

    param_out <- try(as.matrix(param), silent = TRUE)
    if(inherits(param_out, "try-error")){
      return(param_out)
    }

    if(((dim(param_out)[1] != (6*M + 3)) && (dim(param_out)[2] != (6*M + 3)))
       || (dim(param_out)[1] == (6*M + 3) && dim(param_out)[2] != 1)
       || (dim(param_out)[2] == (6*M + 3) && dim(param_out)[1] != 1)
       ) {
      stop(paste0(substitute(param), " must be of dimensions ", (6*M + 3), "x1."))
    }

    if(dim(param_out)[2] == (6*M + 3)){
      param_out <- t(param_out)
    }

  }

  return(param_out)
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
#' @param upper_bound Upper parameter bounds (only for global optimisers)
#' @param lower_bound Lower parameter bounds (only for global optimisers)
#' @param initial_filter_variance Initial Kalman filter state variance
#' @param compute_se Whether or not standard eroors should be computed
#' @param no_of_bootstraps Number of bootstrap repititions
#' @param length_of_bootstraps Length of each bootstrap block
#' @param seed Seed used for bootstrapping
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
quoteDynOptim <- function(
    start,
    X,
    tau,
    ident_mat,
    rel_xtol = 1e-8,
    rel_ftol = 1e-8,
    max_eval = 100000L,
    algorithm = "LN_NELDERMEAD",
    upper_bound = NULL,
    lower_bound = NULL,
    initial_filter_variance = 1000,
    compute_se = TRUE,
    no_of_bootstraps = 100,
    length_of_bootstraps = NULL,
    seed = 21052025,
    hessian = TRUE,
    step_size = 1e-8,
    verbose = FALSE
) {

  # Start error handling block #

  N <- dim(X)[2]
  M <- N / 2

  if (N %% 2 != 0) stop("Number of columns of X must be even.")

  T <- dim(X)[1]

  # Convert inputs to matrices #

  errorOccurred <- FALSE

  if(!is.zoo(X) && !is.xts(X)){
    stop("X must be a zoo/xts object")
  }

  # Create date change indicator
  day_date <- day(time(X))
  day_change_ind <- c(0, as.integer(day_date[-1] == day_date[-length(day_date)]))
  date_change_index <- unlist(lapply(which(day_change_ind == 0), function(i) seq(i, i + 49)))
  if(max(date_change_index) >= length(day_date)){
    stop("Not enough observations for the last day of the data set. Each day must contain at least 50 Observations.")
  }
  day_change_ind[date_change_index] <- 0
  day_change_ind <- as.matrix(day_change_ind)

  XR <- try(coredata(X), silent = TRUE)
  if (inherits(XR, "try-error")) {
    errorOccurred <- TRUE
  }

  # Logic block for handling different input types for the start vector #

  startR <- checkParameterVector(start, M)
  if (inherits(startR, "try-error")) {
    errorOccurred <- TRUE
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

  # Check bounds #

  if(is.null(lower_bound) && is.null(upper_bound)){ # Logic block for unbound optimisation

    # Check whether bounds are disabled but global optimisation should be used
    if(substr(algorithm, 1, 2) == "GN"){
      stop(paste0(algorithm, " corresponds to a global optimisation procedure. Upper and lower parameter bounds must be set."))
    }

    # Create empty dummy matrices so C++ doesn't cry
    upper_bound_R <- as.matrix(-1, dim(startR)[1], 1)
    lower_bound_R <- as.matrix(-1, dim(startR)[1], 1)

  }else if((is.null(lower_bound) && !is.null(upper_bound))
           || (!is.null(lower_bound) && is.null(upper_bound))){ # Logic block for handling misspecification of the bounds

    stop("Either both or no bounds must be provided. Currently, one of the two is parsed as NULL.")

  }else if(!is.null(lower_bound) && !is.null(upper_bound)){ # Logic block for the case where bounds are correctly provided

    # Check whether bounds are enabled but local optimisation should be used
    if(substr(algorithm, 1, 2) == "LN"){
      warning(paste0("Warning!", algorithm, " corresponds to a local optimisation procedure. Upper and lower parameter bounds will be ignored."))
    }

    # Check the validity of the provided upper bound #

    upper_bound_R <- checkParameterVector(upper_bound, M)
    if (inherits(upper_bound_R, "try-error")) {
      errorOccurred <- TRUE
    }

    lower_bound_R <- checkParameterVector(lower_bound, M)
    if (inherits(lower_bound_R, "try-error")) {
      errorOccurred <- TRUE
    }

    # If an error occurred in the transformation of the bounds, terminate
    if (errorOccurred) {
      stop("An error occurred in matrix conversion. \"upper_bound\" and \"lower_bound\" must be convertable to matrices. \n")
    }

  }

  # Check initial KF state variance

  if(initial_filter_variance <= 0){
    stop("The initial Kalman filter variance must be strictly positive.")
  }

  # Bounds check for the parameters for the computation of the bootstrap parameters #

  if(compute_se == 1){
    compute_se <- TRUE
  }else if(compute_se == 0){
    compute_se <- FALSE
  }else if(!is.logical(compute_se)){
    stop("compute_se must be boolean.")
  }

  if(compute_se){

    if(!is.null(no_of_bootstraps) && !is.na(no_of_bootstraps) && !is.infinite(no_of_bootstraps)){
      if(no_of_bootstraps <= 0){
        stop("no_of_bootstraps must be strictly positive!")
      }
    }else{
      stop("no_of_bootstraps cannot be set to NULL/NaN/Inf!")
    }

    if(is.null(length_of_bootstraps)){
      length_of_bootstraps <- floor(T^(1/3))
      cat("length_of_bootstraps has been set to default floor(T^(1/3)).")
    }

    if(!is.na(length_of_bootstraps) && !is.infinite(length_of_bootstraps)){
      if(length_of_bootstraps <= 0){
        stop("length_of_bootstraps must be strictly positive!")
      }else if(length_of_bootstraps >= T){
        stop("length_of_bootstraps must be smaller dim(X)[1]!")
      }
    }else{
      stop("length_of_bootstraps cannot be set to NULL/NaN/Inf!")
    }

    if(!is.null(seed) && !is.na(seed) && !is.infinite(seed)){
      if(seed < 0){
        stop("seed must be strictly positive!")
      }else if(seed >= 2^32 - 1){
        stop("seed must be smaller then 2^32 - 1")
      }
    }else{
      stop("seed cannot be set to NULL/NaN/Inf!")
    }
  }

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

  out <- .Call('_QuoteDynamics_FastOptim',
               PACKAGE = 'QuoteDynamics',
               startR,
               XR,
               tauR,
               ident_mat_R,
               day_change_ind,
               rel_xtol,
               rel_ftol,
               max_eval,
               algorithm_id,
               upper_bound_R,
               lower_bound_R,
               initial_filter_variance,
               compute_se,
               no_of_bootstraps,
               length_of_bootstraps,
               seed,
               hessian,
               step_size,
               verbose)

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

  # Compute standard errors, 95% confidence intervals, and p-test statistics
  boot_straps_sorted <- out$boot_straps
  ses <- c()
  p_vals <- c()
  cis <- matrix(NaN, nrow(boot_straps_sorted), 2)

  for(n in seq_len(nrow(boot_straps_sorted))) {
    boot_straps_sorted[n, ] <- sort(boot_straps_sorted[n, ]) # Sort bootstrap values for confidence intervals
    ses[n] <- sqrt(var(boot_straps_sorted[n, ])) # Compute standard errors
    cis[n, 1] <- quantile(boot_straps_sorted[n, ], probs = 0.025) # Extract lower CI
    cis[n, 2] <- quantile(boot_straps_sorted[n, ], probs = 0.975) # Extract upper CI

    # Compute the p-value
    p_plus  <- mean(boot_straps_sorted[n, ] >= 0)
    p_minus <- mean(boot_straps_sorted[n, ] <= 0)
    p_vals[n] <- max(min(2 * min(p_plus, p_minus), 1), 0)
  }

  # Naming the output #

  # Construct name vector
  names <- c(paste0("spread", 1:M),
             paste0("alpha", 1:(2*M)),
             paste0("Cmat", 1:(3*M)),
             "sigma",
             "delta1",
             "delta2")

  # Naming the objects
  names(out$estimate) <- names
  if(compute_se){
    names(ses) <- names
    rownames(cis) <- names
    colnames(cis) <- c("lower 95%-CI", "upper 95%-CI")
    names(p_vals) <- names
    rownames(out$boot_straps) <- names
  }

  out_sthree <- structure(list(
    call = match.call(),
    minf = out$min_val,
    estimate = out$estimate,
    ses = ses,
    CIs = cis,
    p_vals = p_vals,
    boot_straps = out$boot_straps,
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
  cat("\n")
  cat("Call: ")
  print(object$call)
  cat("\n")
  cat("##############################################\n")
  cat("\n")
  cat("NLopt Return Code:", object$optimisation_status, "\n")
  cat("Minimum Value of Obj.Func.:", object$min_val, "\n")
  cat("Number of Evaluations:", object$eval_count, "\n")
  cat("Estimated Parameters:\n\n")
  if (isTRUE(object$call$compute_se)) { # Logical block for when ses have been computed
    # Construct significance stars
    stars <- symnum(object$p_vals, corr = FALSE, cutpoints = c(0, .001, .01, .05, .1, 1),
                    symbols = c("***", "**", "*", ".", " "))

    table <- data.frame(Estimate = format(object$estimate, scientific = TRUE, digits = 4),
                        `95%-CI`= paste0("[", format(object$CIs[, 1], scientific = TRUE, digits = 4),
                                         ", ", format(object$CIs[, 2], scientific = TRUE, digits = 4),
                                         "]"),
                        Sign. = stars,
                        SE = format(object$ses, scientific = TRUE, digits = 4),
                        check.names = FALSE)
    colnames(table)[3] <- ""
    print(table)
  } else {
    print(data.frame(Estimate = object$estimate))
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

  if (N %% 2 != 0) stop("Number of columns of X must be even.")

  T <- dim(X)[1]

  # Convert inputs to matrices #

  errorOccurred <- FALSE

  if(!is.zoo(X) && !is.xts(X)){
    stop("X must be a zoo/xts object")
  }

  # Create date change indicator
  day_date <- day(time(X))
  day_change_ind <- c(0, as.integer(day_date[-1] == day_date[-length(day_date)]))
  date_change_index <- unlist(lapply(which(day_change_ind == 0), function(i) seq(i, i + 49)))
  if(max(date_change_index) >= length(day_date)){
    stop("Not enough observations for the last day of the data set. Each day must contain at least 50 Observations.")
  }
  day_change_ind[date_change_index] <- 0
  day_change_ind <- as.matrix(day_change_ind)

  XR <- try(coredata(X), silent = TRUE)
  if (inherits(XR, "try-error")) {
    errorOccurred <- TRUE
  }

  # Logic block for handling different input types for the start vector #

  startR <- checkParameterVector(start, M)
  if (inherits(startR, "try-error")) {
    errorOccurred <- TRUE
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

  .Call('_QuoteDynamics_objFunctionCpp', PACKAGE = 'QuoteDynamics', startR, XR, tauR, ident_mat_R, day_change_ind, verbose)

}
