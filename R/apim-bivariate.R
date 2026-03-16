#' Empirical transition counts for the bivariate dyadic model
#'
#' Computes empirical transition counts for the bivariate dyadic model (two variables).
#' The current implementation supports \code{states = 2} only.
#'
#' @param chainFM_V1,chainSM_V1 Vectors of observed states for variable 1 (FM and SM).
#' @param chainFM_V2,chainSM_V2 Vectors of observed states for variable 2 (FM and SM).
#' @param states A single integer. Currently only \code{2} is supported.
#' @returns An integer matrix of counts with 16 rows and 2 columns (when \code{states = 2}).
#' @examples
#' chainFM_V1 <- c(1L, 2L, 1L, 2L, 2L, 1L)
#' chainSM_V1 <- c(2L, 1L, 2L, 1L, 1L, 2L)
#' chainFM_V2 <- c(1L, 1L, 2L, 2L, 1L, 2L)
#' chainSM_V2 <- c(2L, 2L, 1L, 1L, 2L, 1L)
#' emp <- countEmpBivariate(chainFM_V1, chainSM_V1, chainFM_V2, chainSM_V2, states = 2L)
#' dim(emp)
#' @export
countEmpBivariate <- function(chainFM_V1, chainSM_V1, chainFM_V2, chainSM_V2, states = 2L) {

  # Bivariate theory currently implemented for states = 2 only
  ok_states <- length(states) == 1L && is.finite(states) &&
    states == as.integer(states) && as.integer(states) == 2L
  if (!ok_states) stop("bivariate functions currently support states = 2 only.")
  states <- 2L

  # No missing values allowed
  if (anyNA(chainFM_V1) || anyNA(chainSM_V1) || anyNA(chainFM_V2) || anyNA(chainSM_V2)) {
    stop("chains must not contain NA.")
  }

  # All four sequences must have the same length (and allow at least one transition)
  n <- length(chainFM_V1)
  if (n != length(chainSM_V1) || n != length(chainFM_V2) || n != length(chainSM_V2)) {
    stop("chains must have the same length.")
  }
  if (n < 2L) stop("chains must have length >= 2.")

  # With states = 2, values must be in {1, 2} (and integer-like)
  bad <- any(chainFM_V1 != as.integer(chainFM_V1)) || any(chainSM_V1 != as.integer(chainSM_V1)) ||
    any(chainFM_V2 != as.integer(chainFM_V2)) || any(chainSM_V2 != as.integer(chainSM_V2)) ||
    any(chainFM_V1 < 1L | chainFM_V1 > 2L) || any(chainSM_V1 < 1L | chainSM_V1 > 2L) ||
    any(chainFM_V2 < 1L | chainFM_V2 > 2L) || any(chainSM_V2 < 1L | chainSM_V2 > 2L)
  if (bad) stop("with states = 2, chain values must be integers in {1, 2}.")

  # Number of transitions and next FM state for variable 1 (columns)
  chainCount <- n - 1L
  col <- chainFM_V1[2L:(chainCount + 1L)]

  # Current dyad states for both variables
  bfm1 <- chainFM_V1[1L:chainCount]; bsm1 <- chainSM_V1[1L:chainCount]
  bfm2 <- chainFM_V2[1L:chainCount]; bsm2 <- chainSM_V2[1L:chainCount]

  # Legacy dimensions: 16 rows (states = 2) and 2 columns
  nrow_out <- 4L * states * states

  # Map (V1 state, V2 state) to a row index in 1:16, then flatten (row, col) for tabulate()
  row <- states^2 * (states * (bfm1 - 1L) + (bsm1 - 1L)) +
    states * (bfm2 - 1L) + (bsm2 - 1L) + 1L
  idx <- row + (col - 1L) * nrow_out

  tab <- tabulate(idx, nbins = nrow_out * states)
  count <- matrix(tab, nrow = nrow_out, ncol = states)

  count
}


#' Classify the bivariate dependence case
#'
#' Classifies the bivariate case as \code{"trivial"}, \code{"univariate"},
#' \code{"partial"}, or \code{"complete"} using two chi-squared tests against
#' constrained models (states = 2 only).
#'
#' @param empirical An empirical bivariate count matrix (must be 16x2; states = 2).
#' @param alpha A single number in (0, 1) giving the significance level.
#' @returns A list with components \code{testUnivariate}, \code{testPartial}, and \code{case}.
#' @examples
#' chainFM_V1 <- c(1L, 2L, 1L, 2L, 2L, 1L)
#' chainSM_V1 <- c(2L, 1L, 2L, 1L, 1L, 2L)
#' chainFM_V2 <- c(1L, 1L, 2L, 2L, 1L, 2L)
#' chainSM_V2 <- c(2L, 2L, 1L, 1L, 2L, 1L)
#' emp <- countEmpBivariate(chainFM_V1, chainSM_V1, chainFM_V2, chainSM_V2, states = 2L)
#' bivariateCase(emp, alpha = 0.05)
#' @export
bivariateCase <- function(empirical, alpha = 0.05){

  # Significance level
  if (!(length(alpha) == 1L && is.finite(alpha) && alpha > 0 && alpha < 1)) {
    stop("alpha must be a single number in (0, 1).")
  }

  # Bivariate case currently supports states = 2 only (16x2 empirical counts)
  if (!is.matrix(empirical) || nrow(empirical) != 16L || ncol(empirical) != 2L) {
    stop("bivariate functions currently support states = 2 only (empirical must be a 16x2 matrix).")
  }
  if (anyNA(empirical) || any(empirical < 0)) {
    stop("empirical must contain non-negative counts with no NA.")
  }

  # Build constrained theoretical count matrices (G-family)
  g <- countTheoBivariateG(empirical)
  theoA1 <- g[[1]]
  theoB1 <- g[[2]]

  # Two chi-squared tests used for case classification
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

  list(testUnivariate = testUnivariate, testPartial = testPartial, case = case)
}


#' Select the best partial bivariate pattern by AIC
#'
#' Compares the partial bivariate patterns B1/B2/B3 using AIC and returns the
#' selected pattern.
#'
#' @param empirical An empirical bivariate count matrix (must be 16x2; states = 2).
#' @details Requires a bivariate empirical count matrix for \code{states = 2}
#'   (output of \code{\link{countEmpBivariate}}).
#' @returns A list with components \code{aic} (a data frame) and \code{pattern}
#'   (the selected pattern label).
#' @examples
#' chainFM_V1 <- c(1L, 2L, 1L, 2L, 2L, 1L)
#' chainSM_V1 <- c(2L, 1L, 2L, 1L, 1L, 2L)
#' chainFM_V2 <- c(1L, 1L, 2L, 2L, 1L, 2L)
#' chainSM_V2 <- c(2L, 2L, 1L, 1L, 2L, 1L)
#' emp <- countEmpBivariate(chainFM_V1, chainSM_V1, chainFM_V2, chainSM_V2, states = 2L)
#' partialPattern(emp)
#' @export
partialPattern <- function(empirical){

  # Bivariate case currently supports states = 2 only (16x2 empirical counts)
  if (!is.matrix(empirical) || nrow(empirical) != 16L || ncol(empirical) != 2L) {
    stop("bivariate functions currently support states = 2 only (empirical must be a 16x2 matrix).")
  }
  if (anyNA(empirical) || any(empirical < 0)) {
    stop("empirical must contain non-negative counts with no NA.")
  }

  # Constrained theoretical matrices for partial patterns B1/B2/B3
  g <- countTheoBivariateG(empirical)  # list(theoA1, theoB1)
  p <- countTheoBivariateP(empirical)  # list(theoB2, theoB3)

  theoB1 <- g[[2]]
  theoB2 <- p[[1]]
  theoB3 <- p[[2]]

  # AIC values for each candidate pattern (k differs by constraint structure)
  aic_vals <- c(
    "partial actor partner (B1)" = aicBivariate(population = theoB1, empirical = empirical, test = "duo"),
    "partial actor (B2)"         = aicBivariate(population = theoB2, empirical = empirical, test = "single"),
    "partial partner (B3)"       = aicBivariate(population = theoB3, empirical = empirical, test = "single")
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

  list(aic = aic, pattern = pattern)
}

#' Select the best complete bivariate pattern by AIC
#'
#' Compares complete bivariate patterns (C, D1–D4, E1–E4) using AIC and returns
#' the selected pattern.
#'
#' @param empirical An empirical bivariate count matrix (must be 16x2; states = 2).
#' @details Requires a bivariate empirical count matrix for \code{states = 2}
#'   (output of \code{\link{countEmpBivariate}}).
#' @returns A list with components \code{aic} (a data frame with columns
#'   \code{pattern}, \code{matrix}, \code{aic}) and \code{pattern} (the selected
#'   pattern label).
#' @examples
#' chainFM_V1 <- c(1L, 2L, 1L, 2L, 2L, 1L)
#' chainSM_V1 <- c(2L, 1L, 2L, 1L, 1L, 2L)
#' chainFM_V2 <- c(1L, 1L, 2L, 2L, 1L, 2L)
#' chainSM_V2 <- c(2L, 2L, 1L, 1L, 2L, 1L)
#' emp <- countEmpBivariate(chainFM_V1, chainSM_V1, chainFM_V2, chainSM_V2, states = 2L)
#' completePattern(emp)
#' @export
completePattern <- function(empirical){

  # Bivariate case currently supports states = 2 only (16x2 empirical counts)
  if (!is.matrix(empirical) || nrow(empirical) != 16L || ncol(empirical) != 2L) {
    stop("bivariate functions currently support states = 2 only (empirical must be a 16x2 matrix).")
  }
  if (anyNA(empirical) || any(empirical < 0)) {
    stop("empirical must contain non-negative counts with no NA.")
  }

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

  list(aic = aic, pattern = pattern)
}
