#' rOpenSci statistical software standards
#'
#' This file records compliance annotations for the rOpenSci statistical
#' software review standards that apply to dyadicMarkov. Standards that are not
#' applicable are listed separately with brief justifications.
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstats {G1.4} All exported functions are documented using roxygen2, and devtools::document() generates the corresponding Rd files.
#' @srrstats {G1.4a} All internal helper functions are documented with roxygen2 blocks and @noRd tags, including validation helpers, likelihood ratio helpers, AIC helpers, and S3 print helper functions.
#' @srrstats {G3.0} Numeric equality checks are used only for integer coded categorical states and scalar integer validation. dyadicMarkov does not compare floating point statistical estimates for exact equality.
#' @srrstats {G5.1} The package created example data sets are exported as package data and are used in shipped tests to verify their structure, state coding, and expected workflow classifications.
#' @srrstats {G5.3} Shipped tests check that returned empirical counts and MLE probability matrices contain finite, non negative values and valid probabilities after row normalization.
#' @srrstats {G5.4} Correctness tests use fixed hand constructed dyadic sequences and empirical count matrices with known outputs for countEmp(), countEmpBivariate(), and mleEstimation().
#' @srrstats {G5.4a} Because the dyadic transition matrix method is package specific, implementation correctness is tested against hand computed cases, including exact univariate counts, exact bivariate counts, zero row MLE behavior, and all identical valid chains.
#' @srrstats {G5.6} Parameter recovery is tested for mleEstimation() using empirical transition count matrices with known probabilities obtained by row normalization.
#' @srrstats {G5.6a} Parameter recovery tests compare recovered transition probabilities with known expected probabilities using an explicit numerical tolerance.
#' @srrstats {G5.7} Scaling behaviour tests demonstrate that mleEstimation() returns stable transition probabilities when empirical transition counts are multiplied by a constant factor.
#' @noRd
NULL

#' NA_standards
#'
#' The following standards are recorded as not applicable to the current package
#' scope. Each annotation gives the reason for non applicability.
#'
#' @srrstatsNA {G1.5} dyadicMarkov does not currently make algorithmic performance claims requiring reproduction code. The package examples and tests focus on correctness of transition counting, estimation, and pattern identification.
#' @srrstatsNA {G1.6} dyadicMarkov does not currently make comparative performance claims against alternative R implementations.
#' @srrstatsNA {G2.3} dyadicMarkov does not use free form univariate character parameters to control statistical algorithms.
#' @srrstatsNA {G2.3a} Not applicable because dyadicMarkov does not use free form univariate character parameters requiring match.arg().
#' @srrstatsNA {G2.3b} Not applicable because dyadicMarkov does not use case sensitive free form character parameters.
#' @srrstatsNA {G2.4c} dyadicMarkov does not require conversion of inputs to character for its statistical algorithms.
#' @srrstatsNA {G2.4d} dyadicMarkov does not require conversion of inputs to factor for its statistical algorithms.
#' @srrstatsNA {G2.4e} dyadicMarkov does not require conversion from factor inputs for its statistical algorithms.
#' @srrstatsNA {G2.5} dyadicMarkov does not expect factor inputs; categorical states are represented through integer coded state vectors.
#' @srrstatsNA {G2.7} dyadicMarkov operates on state vectors and transition count matrices, not on general tabular time series containers.
#' @srrstatsNA {G2.9} dyadicMarkov does not perform lossy type conversions or metadata altering conversions such as factor to character or spatial metadata removal.
#' @srrstatsNA {G2.10} dyadicMarkov does not extract single columns from tabular inputs as part of its statistical algorithms.
#' @srrstatsNA {G2.11} dyadicMarkov does not accept data.frame like tabular inputs with non standard column classes as primary statistical inputs.
#' @srrstatsNA {G2.12} dyadicMarkov does not accept data.frame like tabular inputs with list columns as primary statistical inputs.
#' @srrstatsNA {G2.14b} dyadicMarkov does not ignore missing data because doing so would alter empirical transition counts.
#' @srrstatsNA {G2.14c} dyadicMarkov does not impute missing categorical states because imputation would alter empirical transition counts and inferred transition patterns.
#' @srrstatsNA {G3.1} dyadicMarkov does not rely on covariance calculations.
#' @srrstatsNA {G3.1a} Not applicable because dyadicMarkov does not use covariance calculations.
#' @srrstatsNA {G4.0} dyadicMarkov does not write statistical outputs to local files.
#' @srrstatsNA {G5.5} dyadicMarkov's core algorithms are deterministic and do not require random seeds for correctness tests.
#' @srrstatsNA {G5.6b} dyadicMarkov's core algorithms are deterministic and do not involve random seeds or random initial conditions.
#' @srrstatsNA {G5.9} Noise susceptibility tests are not applicable to dyadicMarkov's categorical state inputs because adding numerical noise to integer coded categories would produce invalid states rather than meaningful perturbed data.
#' @srrstatsNA {G5.9a} Not applicable because adding floating point noise to categorical state labels is not meaningful for dyadic transition sequence analysis.
#' @srrstatsNA {G5.9b} Not applicable because dyadicMarkov's core algorithms are deterministic and do not depend on random seeds or random initial conditions.
#' @srrstatsNA {G5.10} dyadicMarkov currently does not require extended tests beyond the regular testthat suite because the core examples and correctness checks are lightweight.
#' @srrstatsNA {G5.11} dyadicMarkov does not require large external data sets or external assets for extended tests.
#' @srrstatsNA {G5.11a} Not applicable because dyadicMarkov does not download external data for tests.
#' @srrstatsNA {G5.12} Not applicable because dyadicMarkov does not currently require extended tests with special platform, memory, runtime, or external artifact conditions.
#' @srrstatsNA {G2.4b} dyadicMarkov does not require conversion of inputs to continuous numeric values. Categorical states are integer coded, empirical transition counts are checked as numeric matrices, and continuous valued statistical inputs are not part of the package API.
#' @srrstatsNA {G5.0} Not applicable because no external standard reference data set exists for the package specific dyadic categorical transition matrix methods implemented in dyadicMarkov; correctness is instead tested using hand constructed data with known properties.
#' @srrstatsNA {G5.4b} Not applicable because dyadicMarkov is not a new implementation of an existing external software method; the package implements a package specific dyadic transition matrix workflow.
#' @srrstatsNA {G5.4c} Not applicable because the shipped correctness tests use hand computed expected values and exported example data rather than relying on external published output snapshots.
#' @srrstatsNA {EA2.0} dyadicMarkov does not accept standard tabular data requiring extensive filtering or joins, so an index column system is not applicable.
#' @srrstatsNA {EA2.1} Not applicable because dyadicMarkov does not use index columns for table filtering or joining operations.
#' @srrstatsNA {EA2.2} Not applicable because dyadicMarkov does not use index columns or table join workflows.
#' @srrstatsNA {EA2.2a} Not applicable because no index-column class system is used or required.
#' @srrstatsNA {EA2.2b} Not applicable because no index-column attribute is used or required.
#' @srrstatsNA {EA2.3} dyadicMarkov does not perform table join operations.
#' @srrstatsNA {EA2.4} dyadicMarkov does not accept multi-tabular input.
#' @srrstatsNA {EA2.5} Not applicable because dyadicMarkov does not accept multi-tabular input or use index columns.
#' @srrstatsNA {EA5.0} dyadicMarkov does not currently implement graphical output functions as part of the package API.
#' @srrstatsNA {EA5.0a} Not applicable because dyadicMarkov does not currently implement graphical output functions.
#' @srrstatsNA {EA5.0b} Not applicable because dyadicMarkov does not currently implement graphical output functions or default colour schemes.
#' @srrstatsNA {EA5.1} Not applicable because dyadicMarkov does not specify typefaces in graphical output.
#' @srrstatsNA {EA5.4} Not applicable because dyadicMarkov does not currently implement graphical output functions.
#' @srrstatsNA {EA5.5} Not applicable because dyadicMarkov does not currently implement graphical output functions with axes or units.
#' @srrstatsNA {EA5.6} Not applicable because dyadicMarkov does not bundle dynamic visualization libraries.
#' @srrstatsNA {EA6.0d} dyadicMarkov does not return data.frame-type tabular objects as primary statistical outputs.
#' @srrstatsNA {EA6.1} Not applicable because dyadicMarkov does not currently implement graphical output functions.
#' @noRd
NULL
