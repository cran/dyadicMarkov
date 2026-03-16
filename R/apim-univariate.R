#' Univariate pattern classification for dyadic Markov chains
#'
#' Computes empirical transition counts, fits the unrestricted model by maximum
#' likelihood, and performs chi-squared goodness-of-fit tests against Actor-only
#' (AM) and Partner-only (PM) constrained models to classify the univariate dyadic
#' pattern.
#'
#' @param chainFM Vector of observed states for the first member (FM).
#' @param chainSM Vector of observed states for the second member (SM).
#' @param states A single integer >= 2 giving the number of states.
#' @param alpha A single number in (0, 1) giving the significance level.
#' @returns A list with two \code{htest} objects (\code{TEST.AM}, \code{TEST.PM})
#'   and a string \code{pattern}.
#' @examples
#' chainFM <- c(1L, 2L, 1L, 2L, 2L, 1L)
#' chainSM <- c(2L, 1L, 2L, 1L, 1L, 2L)
#' univariatePattern(chainFM, chainSM, states = 2L, alpha = 0.05)
#' @export
univariatePattern <- function(chainFM, chainSM, states, alpha = 0.05) {

  # Significance level
  if (!(length(alpha) == 1L && is.finite(alpha) && alpha > 0 && alpha < 1)) {
    stop("alpha must be a single number in (0, 1).")
  }

  # States must be provided explicitly
  if (missing(states) ||
      !(length(states) == 1L && is.finite(states) &&
        states == as.integer(states) && states >= 2L)) {
    stop("states must be provided as a single integer >= 2.")
  }
  states <- as.integer(states)

  # No missing values in the observed sequences
  if (anyNA(chainFM) || anyNA(chainSM)) {
    stop("chains must not contain NA.")
  }

  # Empirical counts under the unrestricted model
  emp <- countEmp(chainFM = chainFM, chainSM = chainSM, states = states)

  # Theoretical counts under AM and PM constraints
  theoAM <- countTheo(empirical = emp, pattern = "AM")
  theoPM <- countTheo(empirical = emp, pattern = "PM")

  # Local chi-square tests against each constrained model
  lrtAM <- lrtLocal(population = theoAM, empirical = emp)
  lrtPM <- lrtLocal(population = theoPM, empirical = emp)

  pvalueAM <- lrtAM[["p.value"]]
  pvaluePM <- lrtPM[["p.value"]]

  # Pattern classification based on which constraints are rejected
  if (is.na(pvalueAM) || is.na(pvaluePM)) {
    type <- NA_character_
  } else if (pvalueAM > alpha && pvaluePM > alpha) {
    type <- "IM (A0)"
  } else if (pvalueAM < alpha && pvaluePM > alpha) {
    type <- "PM (A3)"
  } else if (pvalueAM > alpha && pvaluePM < alpha) {
    type <- "AM (A2)"
  } else {
    type <- "APM (A1)"
  }

  # Wrap results as htest objects with clear method labels
  TEST1 <- .make_htest(
    method = "Chi-squared test, Actor-only model",
    dataName = "Observed vs Estimated",
    alternative = "The unrestricted model fits the data better",
    statistic = lrtAM[["statistic"]],
    df = lrtAM[["parameter"]],
    pValue = lrtAM[["p.value"]]
  )

  TEST2 <- .make_htest(
    method = "Chi-squared test, Partner-only model",
    dataName = "Observed vs Estimated",
    alternative = "The unrestricted model fits the data better",
    statistic = lrtPM[["statistic"]],
    df = lrtPM[["parameter"]],
    pValue = lrtPM[["p.value"]]
  )

  list(TEST.AM = TEST1, TEST.PM = TEST2, pattern = type)
}
