library(devtools)
library(rstudioapi)
library(roxygen2)
library(Rcpp)
library(RcppEigen)

setwd(dirname(getActiveDocumentContext()$path))

# Rcpp::sourceCpp("version_info.cpp")
# build_info()

# Document the package, this updates NAMESPACE and documentation
unlink("src/*.o")
#    – alle plattformspezifischen Shared-Libs
unlink("src/*.so")    # Unix / macOS
unlink("src/*.dll")   # Windows
#    – zusätzlich evtl. vorher gebaute DLLs in .Rproj.user
devtools::clean_dll()

# 3) (optional) Vignetten-Build-Artefakte entfernen
devtools::clean_vignettes()

# Compile Rcpp attributes first
Rcpp::compileAttributes()

# 4) Roxygen → Rd + NAMESPACE neu erzeugen
devtools::document()

# Install the package (this also builds it)
devtools::build()
# devtools::load_all()
# devtools::check()

detach("package:QuoteDynamics", unload = TRUE)

.rs.restartR()

install.packages("../QuoteDynamics_0.0.4.tar.gz", repos = NULL, type = "source")

# Restart R session before loading the package# Restart R session before verbose = loading the package
.rs.restartR()

# Load the package
library(QuoteDynamics)
ls("package:QuoteDynamics")

# Test the function

PrintAlgorithms()
Minimum <- QuoteDynamics::quoteDynOptim(start = start, X = data, tau = tau, xtol = 10e-4,
                            stop_val = 10e-10, algorithm = "LN_NELDERMEAD", hessian = TRUE,
                            step_size = 1e-04, verbose = TRUE)
print(Minimum)
