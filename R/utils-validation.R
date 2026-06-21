#' Validate significance level
#'
#' Internal helper for validating the `alpha` argument.
#'
#' @param alpha A single number in (0, 1) giving the significance level.
#'
#' @return Invisibly returns the validated value.
#'
#' @srrstats {G2.0} Input-validation helpers assert lengths of scalar arguments and sequence inputs. For example, alpha and states must have length one, dyadic chains must have equal lengths, and chains must have length at least two.
#' @srrstats {G2.0a} Function documentation specifies length expectations for inputs, including equal-length dyadic chains, scalar alpha, scalar states, and fixed empirical matrix dimensions.
#' @srrstats {G2.1} Input-validation helpers assert expected input types, including numeric state vectors, numeric empirical count matrices, scalar numeric alpha, and scalar numeric states.
#' @srrstats {G2.1a} Function documentation specifies expected data types for vector and matrix inputs, including numeric/integer-coded state vectors and numeric empirical transition count matrices.
#' @srrstats {G2.13} Input-validation helpers explicitly check for missing data before analytical routines are applied, using anyNA() for dyadic chains and empirical transition count matrices.
#' @srrstats {G2.14} dyadicMarkov handles missing data by rejecting missing categorical states or missing empirical counts before transition count estimation or pattern identification.
#' @srrstats {G2.14a} dyadicMarkov implements the error-on-missing-data option: inputs containing NA values trigger informative errors before analysis proceeds.
#' @srrstats {G2.15} Functions do not assume non-missingness. Missing values are checked in validation helpers before data are passed to the transition count, estimation, or testing routines.
#' @srrstats {G2.16} Undefined numeric values are handled during validation: empirical matrices are checked for finite values, and scalar arguments such as alpha and states must be finite before analysis proceeds.
#'
#' @srrstats {G2.2} Functions restrict inputs to the expected dimensionality. Chain inputs are checked as vectors through length and type validation, univariate empirical inputs must be matrices with states^2 rows and states columns, and bivariate empirical inputs must be 16x2 matrices.
#' @srrstats {G2.4} dyadicMarkov uses explicit validation and limited conversion mechanisms where appropriate. The state-space size is converted to integer after validation, while categorical chain values and empirical counts are checked rather than silently coerced.
#' @srrstats {G2.4a} The states argument is explicitly validated and converted with as.integer() after checking that it is a single finite integer-like value.
#' @srrstats {G2.6} One-dimensional chain inputs are pre-processed through validation helpers that check type, length, missingness, integer coding, and supported state ranges before empirical transition counts are computed.
#' @srrstats {G2.8} The package uses validation and preprocessing routines to ensure that analytical routines receive standardized inputs: integer-coded state vectors for counting functions and numeric empirical transition count matrices for estimation and testing functions.
#' @srrstats {EA2.6} Validation tests demonstrate that vector inputs are checked for length, missing values, finite integer-coded states, and admissible state ranges before analytic routines are applied.
#' @noRd
.validate_alpha <- function(alpha) {
  if (!(is.numeric(alpha) && length(alpha) == 1L &&
        is.finite(alpha) && alpha > 0 && alpha < 1)) {
    stop("alpha must be a single number in (0, 1).", call. = FALSE)
  }

  invisible(alpha)
}


#' Validate number of states
#'
#' Internal helper for validating the number of states.
#'
#' @param states Integer number of categorical states.
#'
#' @return Returns the validated number of states as an integer.
#' @noRd
.validate_states <- function(states) {
  ok_states <- !missing(states) &&
    is.numeric(states) &&
    length(states) == 1L &&
    is.finite(states) &&
    states == floor(states) &&
    states >= 2L &&
    states <= .Machine$integer.max

  if (!ok_states) {
    stop("states must be provided as a single integer >= 2.", call. = FALSE)
  }

  as.integer(states)
}


#' Validate univariate dyadic chains
#'
#' Internal helper for validating categorical dyadic input vectors.
#'
#' @param chainFM Vector of states for the first member.
#' @param chainSM Vector of states for the second member.
#' @param states Integer number of categorical states.
#'
#' @return Invisibly returns `TRUE` if validation succeeds;
#'   otherwise raises an error.
#' @noRd
.validate_univariate_chains <- function(chainFM, chainSM, states) {
  n <- length(chainFM)

  if (!is.numeric(chainFM) || !is.numeric(chainSM)) {
    stop(
      "chain values must be numeric vectors of integer-coded states.",
      call. = FALSE
    )
  }
  if (n != length(chainSM)) {
    stop("chainFM and chainSM must have the same length.", call. = FALSE)
  }
  if (n < 2L) {
    stop("chains must have length >= 2.", call. = FALSE)
  }
  if (anyNA(chainFM) || anyNA(chainSM)) {
    stop("chains must not contain NA.", call. = FALSE)
  }
  if (!all(is.finite(chainFM)) || !all(is.finite(chainSM))) {
    stop("chains must contain finite integer-coded states.", call. = FALSE)
  }

  bad <- any(chainFM != as.integer(chainFM)) ||
    any(chainSM != as.integer(chainSM)) ||
    any(chainFM < 1L | chainFM > states) ||
    any(chainSM < 1L | chainSM > states)

  if (bad) {
    stop("chain values must be integers in 1:states.", call. = FALSE)
  }

  invisible(TRUE)
}


#' Validate empirical transition count matrix
#'
#' Internal helper for validating univariate empirical transition
#' count matrices.
#'
#' @param empirical Empirical transition count matrix.
#'
#' @return Invisibly returns `TRUE` if validation succeeds;
#'   otherwise raises an error.
#' @noRd
.validate_empirical_matrix <- function(empirical) {
  if (!is.matrix(empirical)) {
    stop("empirical transition count input must be a matrix.", call. = FALSE)
  }
  if (!is.numeric(empirical) || anyNA(empirical) ||
      !all(is.finite(empirical)) || any(empirical < 0)) {
    stop(
      "empirical transition count input must contain finite non-negative ",
      "counts with no NA.",
      call. = FALSE
    )
  }

  states <- ncol(empirical)
  if (states < 2L) {
    stop("empirical must have at least 2 columns (states >= 2).", call. = FALSE)
  }
  if (nrow(empirical) != states * states) {
    stop("empirical must have states^2 rows and states columns.", call. = FALSE)
  }

  invisible(TRUE)
}


#' Validate empirical bivariate transition count matrix
#'
#' Internal helper for validating bivariate empirical transition count matrices.
#'
#' @param empirical Empirical bivariate transition count matrix.
#'
#' @return Invisibly returns `TRUE` if validation succeeds;
#'   otherwise raises an error.
#' @noRd
.validate_bivariate_empirical_matrix <- function(empirical) {
  if (!is.matrix(empirical) ||
      nrow(empirical) != 16L ||
      ncol(empirical) != 2L) {
    stop(
      "bivariate functions currently support states = 2 only ",
      "(empirical must be a 16x2 matrix).",
      call. = FALSE
    )
  }
  if (!is.numeric(empirical) || anyNA(empirical) ||
      !all(is.finite(empirical)) || any(empirical < 0)) {
    stop(
      "bivariate empirical matrix must contain finite non-negative ",
      "counts with no NA.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
