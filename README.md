

# QuoteDynamics Package

## Overview
QuoteDynamics is an R package providing fast Kalman Filter and Smoother Maximum Likelihood Estimation using the NLopt optimization library.

## Installation

### Windows Users
Windows users can directly install the package as it includes the necessary NLopt `.dll` file.

### Linux and macOS Users
Linux and macOS users need to install NLopt separately.

#### Installing NLopt on Linux
You can install NLopt on Linux using the package manager. For example, on Ubuntu:

```bash
sudo apt-get install libnlopt-dev
```

#### Installing NLopt on macOS
On macOS, NLopt can be installed using Homebrew:

```bash
brew install nlopt
```

## Usage
To use the QuoteDynamics package, first load it into your R session:

```R
library(QuoteDynamics)
```

You can then call the main function NLoptKFS with appropriate parameters. Here is a working example using the simulated data that is included in the package.

```R
Minimum <-     QuoteDynamics::QuoteDynamics(start = start, X = data, tau = tau, xtol = 10e-4,
                            stop_val = 10e-10, algorithm = "LN_NELDERMEAD", hessian = TRUE,
                            step_size = 1e-04, verbose = FALSE)
print(Minimum)
```

Refer to the package documentation for detailed usage instructions.

## Contact

For any queries or feedback, please contact me.
