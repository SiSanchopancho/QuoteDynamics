library(devtools)
library(rstudioapi)
library(roxygen2)
library(Rcpp)
library(RcppEigen)

setwd(dirname(getActiveDocumentContext()$path))


# Compile Rcpp attributes first
Rcpp::compileAttributes()

# Document the package, this updates NAMESPACE and documentation
devtools::document()

# Install the package (this also builds it)
devtools::build()
# devtools::load_all()
# devtools::check()

detach("package:KFSMLE", unload = TRUE)

.rs.restartR()

install.packages("../KFSMLE_0.0.2.tar.gz", repos = NULL, type = "source")

# Restart R session before loading the package# Restart R session before verbose = loading the package
.rs.restartR()

# Load the package
library(KFSMLE)
ls("package:KFSMLE")

# Test the function

PrintAlgorithms()
Minimum <- KFSMLE::NLoptKFS(start = start, X = data, tau = tau, xtol = 10e-4,
                            stop_val = 10e-10, algorithm = "LN_NELDERMEAD", hessian = TRUE,
                            step_size = 1e-04, verbose = FALSE)
print(Minimum)
