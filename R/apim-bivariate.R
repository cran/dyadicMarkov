#' Construct bivariate dyadic transition dimnames
#'
#' Internal helper used for bivariate empirical count matrices.
#'
#' @param states Number of categorical states.
#'
#' @return A list suitable for use as matrix dimnames.
#' @noRd
.dyadic_bivariate_dimnames <- function(states) {
  grid <- expand.grid(
    second_second = seq_len(states),
    second_first = seq_len(states),
    main_second = seq_len(states),
    main_first = seq_len(states)
  )

  rows <- paste0(
    "mainFM", grid$main_first,
    "_mainSM", grid$main_second,
    "_secondFM", grid$second_first,
    "_secondSM", grid$second_second
  )

  list(rows, paste0("next_", seq_len(states)))
}


# Validate that bivariate counting currently uses binary states.
.validate_bivariate_states <- function(states) {
  ok_states <- !missing(states) &&
    is.numeric(states) &&
    length(states) == 1L &&
    is.finite(states) &&
    states == floor(states) &&
    states == 2L

  if (!ok_states) {
    stop(
      "bivariate functions currently support states = 2 only.",
      call. = FALSE
    )
  }

  2L
}


# Validate that all bivariate chains are numeric.
.validate_bivariate_chain_types <- function(chains) {
  is_num <- vapply(chains, is.numeric, logical(1))

  if (!all(is_num)) {
    stop("chain values must be numeric integers in {1, 2}.", call. = FALSE)
  }

  invisible(chains)
}


# Validate that bivariate chains do not contain missing values.
.validate_bivariate_chain_missing <- function(chains) {
  has_na <- vapply(chains, anyNA, logical(1))

  if (any(has_na)) {
    stop("bivariate chains must not contain NA.", call. = FALSE)
  }

  invisible(chains)
}


# Validate that bivariate chains contain finite numeric values.
.validate_bivariate_chain_finite <- function(chains) {
  finite <- vapply(
    chains,
    function(x) all(is.finite(x)),
    logical(1)
  )

  if (!all(finite)) {
    stop(
      "bivariate chains must contain finite integer-coded states.",
      call. = FALSE
    )
  }

  invisible(chains)
}


# Validate that all bivariate chains have a common usable length.
.validate_bivariate_chain_lengths <- function(chains) {
  lengths <- vapply(chains, length, integer(1))

  if (length(unique(lengths)) != 1L) {
    stop("chains must have the same length.", call. = FALSE)
  }

  n <- lengths[[1L]]

  if (n < 2L) {
    stop("bivariate chains must have length >= 2.", call. = FALSE)
  }

  n
}


# Validate that bivariate chains use integer-coded states in {1, 2}.
.validate_bivariate_chain_values <- function(chains) {
  bad_chain <- function(x) {
    any(x != as.integer(x)) || any(x < 1L | x > 2L)
  }

  bad <- vapply(chains, bad_chain, logical(1))

  if (any(bad)) {
    stop(
      "with states = 2, chain values must be integers in {1, 2}.",
      call. = FALSE
    )
  }

  invisible(chains)
}


#' Empirical transition counts for bivariate dyadic sequences
#'
#' Computes empirical transition counts for bivariate categorical dyadic
#' sequences with two variables. This function currently supports
#' \code{states = 2} only.
#'
#' @param chainFM_V1,chainSM_V1 Vectors of observed states for variable 1
#'   for the first and second member.
#' @param chainFM_V2,chainSM_V2 Vectors of observed states for variable 2
#'   for the first and second member.
#' @param states A single integer. Currently only \code{2} is supported.
#'   Default is 2.
#' @details The bivariate counter currently supports \code{states = 2} only.
#'   Rows represent the previous dyadic states of variable 1 and variable 2.
#'   The implementation uses the row mapping
#'   \code{states^2 * (states * (FM_V1,t - 1) + (SM_V1,t - 1)) +}
#'   \code{states * (FM_V2,t - 1) + (SM_V2,t - 1) + 1}. Columns correspond to
#'   the state of the first member on variable 1 at the next time point,
#'   \eqn{FM_{V1,t+1}}.
#' @returns An integer matrix with class
#'   \code{c("dyadic_counts", "matrix", "array")} with 16 rows and 2 columns
#'   when \code{states = 2}. It remains usable as an ordinary matrix.
#' @examples
#' chainFM_V1 <- c(1L, 2L, 1L, 2L, 2L, 1L)
#' chainSM_V1 <- c(2L, 1L, 2L, 1L, 1L, 2L)
#' chainFM_V2 <- c(1L, 1L, 2L, 2L, 1L, 2L)
#' chainSM_V2 <- c(2L, 2L, 1L, 1L, 2L, 1L)
#' emp <- countEmpBivariate(
#'   chainFM_V1, chainSM_V1, chainFM_V2, chainSM_V2,
#'   states = 2L
#' )
#' dim(emp)
#' @export
countEmpBivariate <- function(
  chainFM_V1, chainSM_V1, chainFM_V2, chainSM_V2,
  states = 2L
) {

  states <- .validate_bivariate_states(states)

  chains <- list(
    chainFM_V1 = chainFM_V1,
    chainSM_V1 = chainSM_V1,
    chainFM_V2 = chainFM_V2,
    chainSM_V2 = chainSM_V2
  )

  .validate_bivariate_chain_types(chains)
  .validate_bivariate_chain_missing(chains)
  .validate_bivariate_chain_finite(chains)

  n <- .validate_bivariate_chain_lengths(chains)

  .validate_bivariate_chain_values(chains)

  # Number of transitions and next FM state for variable 1 (columns)
  chainCount <- n - 1L
  col <- chainFM_V1[2L:(chainCount + 1L)]

  # Current dyad states for both variables
  bfm1 <- chainFM_V1[1L:chainCount]
  bsm1 <- chainSM_V1[1L:chainCount]
  bfm2 <- chainFM_V2[1L:chainCount]
  bsm2 <- chainSM_V2[1L:chainCount]

  # Legacy dimensions: 16 rows (states = 2) and 2 columns
  nrow_out <- 4L * states * states

  # Map (V1 state, V2 state) to a row index in 1:16, then flatten
  # (row, col) for tabulate()
  row <- states^2 * (states * (bfm1 - 1L) + (bsm1 - 1L)) +
    states * (bfm2 - 1L) + (bsm2 - 1L) + 1L

  idx <- row + (col - 1L) * nrow_out

  tab <- tabulate(idx, nbins = nrow_out * states)
  count <- matrix(tab, nrow = nrow_out, ncol = states)
  dimnames(count) <- .dyadic_bivariate_dimnames(states)

  class(count) <- c("dyadic_counts", "matrix", "array")
  count
}


#' Bivariate case identification for dyadic Markov chains
#'
#' Identifies the bivariate case as \code{"trivial"}, \code{"univariate"},
#' \code{"partial"}, or \code{"complete"} using two likelihood-ratio tests
#' against constrained bivariate structures.
#'
#' @param empirical An empirical bivariate count matrix with 16 rows and 2
#'   columns, as returned by \code{\link{countEmpBivariate}}.
#' @srrstats {EA3.1} The package provides standardized comparison of restricted univariate and bivariate transition structures that would otherwise require manual construction of theoretical transition matrices and separate test statistics.
#' @srrstats {EA3.0} The package automates extraction and reporting of dyadic transition counts, MLE transition probabilities, likelihood-ratio comparisons, AIC comparisons, and selected interaction patterns.
#' @param alpha A single number in (0, 1) giving the significance level.
#'   Default is 0.05.
#' @details The returned case corresponds to the global approach of the
#'   bivariate method. It determines whether the sequence analyzed is treated
#'   as trivial, univariate, partial bivariate, or complete bivariate before
#'   the local identification of the pattern of interaction.
#' @returns A list with class \code{c("dyadic_case", "list")} containing
#'   components \code{testUnivariate}, \code{testPartial}, \code{case}, and
#'   metadata fields \code{alpha} and \code{call}. It remains usable as an
#'   ordinary list.
#' @examples
#' chainFM_V1 <- c(1L, 2L, 1L, 2L, 2L, 1L)
#' chainSM_V1 <- c(2L, 1L, 2L, 1L, 1L, 2L)
#' chainFM_V2 <- c(1L, 1L, 2L, 2L, 1L, 2L)
#' chainSM_V2 <- c(2L, 2L, 1L, 1L, 2L, 1L)
#' emp <- countEmpBivariate(
#'   chainFM_V1, chainSM_V1, chainFM_V2, chainSM_V2,
#'   states = 2L
#' )
#' bivariateCase(emp, alpha = 0.05)
#' @export
bivariateCase <- function(empirical, alpha = 0.05) {

  .validate_alpha(alpha)

  .validate_bivariate_empirical_matrix(empirical)

  # Build constrained theoretical count matrices (G-family)
  g <- countTheoBivariateG(empirical)
  theoA1 <- g[[1]]
  theoB1 <- g[[2]]

  # Two likelihood-ratio tests used for case identification
  testUnivariate <- bivariateTest(population = theoA1, empirical = empirical)
  testPartial    <- bivariateTest(population = theoB1, empirical = empirical)

  pA1 <- testUnivariate[["p.value"]]
  pB1 <- testPartial[["p.value"]]

  # Map test outcomes to the final case
  if (is.na(pA1) || is.na(pB1)) {
    case <- NA_character_
  } else if (pA1 > alpha && pB1 > alpha) {
    case <- "trivial"
  } else if (pA1 > alpha && pB1 < alpha) {
    case <- "univariate"
  } else if (pA1 < alpha && pB1 > alpha) {
    case <- "partial"
  } else {
    case <- "complete"
  }

  out <- list(
    testUnivariate = testUnivariate,
    testPartial = testPartial,
    case = case,
    alpha = alpha,
    call = match.call()
  )
  class(out) <- c("dyadic_case", "list")
  out
}


#' Partial bivariate pattern identification by AIC
#'
#' Compares the partial bivariate patterns B1, B2, and B3 using AIC and returns
#' the selected pattern.
#'
#' @param empirical An empirical bivariate count matrix with 16 rows and 2
#'   columns, as returned by \code{\link{countEmpBivariate}}.
#' @details Conditional on the partial bivariate case, AIC is used to select
#'   among the B1, B2, and B3 structures.
#' @returns A list with class \code{c("dyadic_pattern", "list")} containing
#'   components \code{aic} (a data frame with candidate patterns and AIC
#'   values), \code{pattern} (the selected pattern label), and \code{call}. It
#'   remains usable as an ordinary list.
#' @examples
#' chainFM_V1 <- c(1L, 2L, 1L, 2L, 2L, 1L)
#' chainSM_V1 <- c(2L, 1L, 2L, 1L, 1L, 2L)
#' chainFM_V2 <- c(1L, 1L, 2L, 2L, 1L, 2L)
#' chainSM_V2 <- c(2L, 2L, 1L, 1L, 2L, 1L)
#' emp <- countEmpBivariate(
#'   chainFM_V1, chainSM_V1, chainFM_V2, chainSM_V2,
#'   states = 2L
#' )
#' partialPattern(emp)
#' @export
partialPattern <- function(empirical) {

  .validate_bivariate_empirical_matrix(empirical)

  # Constrained theoretical matrices for partial patterns B1/B2/B3
  g <- countTheoBivariateG(empirical)  # list(theoA1, theoB1)
  p <- countTheoBivariateP(empirical)  # list(theoB2, theoB3)

  theoB1 <- g[[2]]
  theoB2 <- p[[1]]
  theoB3 <- p[[2]]

  # AIC values for each candidate pattern (k differs by constraint structure)
  aic_vals <- c(
    "partial actor partner (B1)" = aicBivariate(
      population = theoB1,
      empirical = empirical,
      test = "duo"
    ),
    "partial actor (B2)" = aicBivariate(
      population = theoB2,
      empirical = empirical,
      test = "single"
    ),
    "partial partner (B3)" = aicBivariate(
      population = theoB3,
      empirical = empirical,
      test = "single"
    )
  )

  # Summary table
  aic <- data.frame(
    pattern = names(aic_vals),
    matrix  = c("B1", "B2", "B3"),
    aic     = as.numeric(aic_vals),
    row.names = NULL
  )

  # Select best (smallest AIC)
  pattern <- names(aic_vals)[which.min(aic_vals)]

  out <- list(
    aic = aic,
    pattern = pattern,
    call = match.call()
  )
  class(out) <- c("dyadic_pattern", "list")
  out
}

#' Complete bivariate pattern identification by AIC
#'
#' Compares the complete bivariate patterns C, D1--D4, and E1--E4 using AIC and
#' returns the selected pattern.
#'
#' @param empirical An empirical bivariate count matrix with 16 rows and 2
#'   columns, as returned by \code{\link{countEmpBivariate}}.
#' @details Conditional on the complete bivariate case, AIC is used to select
#'   among the C, D1--D4, and E1--E4 structures.
#' @returns A list with class \code{c("dyadic_pattern", "list")} containing
#'   components \code{aic} (a data frame with columns \code{pattern},
#'   \code{matrix}, and \code{aic}), \code{pattern} (the selected pattern
#'   label), and \code{call}. It remains usable as an ordinary list.
#' @examples
#' chainFM_V1 <- c(1L, 2L, 1L, 2L, 2L, 1L)
#' chainSM_V1 <- c(2L, 1L, 2L, 1L, 1L, 2L)
#' chainFM_V2 <- c(1L, 1L, 2L, 2L, 1L, 2L)
#' chainSM_V2 <- c(2L, 2L, 1L, 1L, 2L, 1L)
#' emp <- countEmpBivariate(
#'   chainFM_V1, chainSM_V1, chainFM_V2, chainSM_V2,
#'   states = 2L
#' )
#' completePattern(emp)
#' @export
completePattern <- function(empirical) {

  .validate_bivariate_empirical_matrix(empirical)

  # Constrained theoretical matrices for complete patterns
  c3 <- countTheoBivariateC3(empirical)  # D1..D4
  c2 <- countTheoBivariateC2(empirical)  # E1..E4

  pops <- c(
    list(C = empirical),
    stats::setNames(c3, c("D1","D2","D3","D4")),
    stats::setNames(c2, c("E1","E2","E3","E4"))
  )

  # Number of free parameters for each structure
  ks <- c(
    C  = 16L,
    D1 = 8L, D2 = 8L, D3 = 8L, D4 = 8L,
    E1 = 4L, E2 = 4L, E3 = 4L, E4 = 4L
  )

  # AIC for each candidate
  aics <- .aicBivariate_many(empirical = empirical, populations = pops, ks = ks)

  # Human-readable labels (in the same order as pops/ks)
  labels <- c(
    "complete actor partner on both (C)",
    "complete partner on the main, actor partner on the second (D1)",
    "complete actor on the main, actor partner on the second (D2)",
    "complete actor partner on the main, partner on the second (D3)",
    "complete actor partner on the main, actor on the second (D4)",
    "complete partner on both (E1)",
    "complete partner on the main, actor on the second (E2)",
    "complete actor on the main, partner on the second (E3)",
    "complete actor on both (E4)"
  )

  codes <- names(ks)
  aic_vals <- as.numeric(aics[codes])

  aic <- data.frame(
    pattern = labels,
    matrix  = codes,
    aic     = aic_vals,
    row.names = NULL
  )

  # Select best (smallest AIC)
  best <- which.min(aic_vals)
  pattern <- labels[best]

  out <- list(
    aic = aic,
    pattern = pattern,
    call = match.call()
  )
  class(out) <- c("dyadic_pattern", "list")
  out
}
