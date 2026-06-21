#' Univariate pattern identification for dyadic Markov chains
#'
#' Computes empirical transition counts, estimates transition probabilities by
#' maximum likelihood, and performs likelihood-ratio tests against the
#' actor-only
#' and partner-only constrained models to identify the univariate pattern of
#' interaction.
#'
#' @param chainFM Vector of observed states for the first member (FM).
#' @param chainSM Vector of observed states for the second member (SM).
#' @param states A single integer >= 2 giving the number of states.
#' @param alpha A single number in (0, 1) giving the significance level.
#'   Default is 0.05.
#' @details Pattern labels summarize which structure is retained by the tests:
#'   \code{IM (A0)} denotes an independence pattern, \code{APM (A1)} an
#'   actor-partner pattern, \code{AM (A2)} an actor-only pattern, and
#'   \code{PM (A3)} a partner-only pattern.
#' @returns A list with class \code{c("dyadic_pattern", "list")} containing two
#'   \code{htest} objects (\code{TEST.AM}, \code{TEST.PM}), a string
#'   \code{pattern}, and metadata fields \code{alpha}, \code{states}, and
#'   \code{call}. It remains usable as an ordinary list.
#' @examples
#' chainFM <- c(1L, 2L, 1L, 2L, 2L, 1L)
#' chainSM <- c(2L, 1L, 2L, 1L, 1L, 2L)
#' univariatePattern(chainFM, chainSM, states = 2L, alpha = 0.05)
#' @export
univariatePattern <- function(chainFM, chainSM, states, alpha = 0.05) {

  .validate_alpha(alpha)

  states <- .validate_states(states)

  .validate_univariate_chains(chainFM, chainSM, states)

  # Empirical counts under the unrestricted model
  emp <- countEmp(chainFM = chainFM, chainSM = chainSM, states = states)

  # Theoretical counts under AM and PM constraints
  theoAM <- countTheo(empirical = emp, pattern = "AM")
  theoPM <- countTheo(empirical = emp, pattern = "PM")

  # Local likelihood-ratio tests against each constrained model
  lrtAM <- lrtLocal(population = theoAM, empirical = emp)
  lrtPM <- lrtLocal(population = theoPM, empirical = emp)

  pvalueAM <- lrtAM[["p.value"]]
  pvaluePM <- lrtPM[["p.value"]]

  # Pattern identification based on which constraints are rejected
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
    method = "Likelihood-ratio test, Actor-only model",
    dataName = "Observed vs Estimated",
    alternative = "The unrestricted model fits the data better",
    statistic = lrtAM[["statistic"]],
    df = lrtAM[["parameter"]],
    pValue = lrtAM[["p.value"]]
  )

  TEST2 <- .make_htest(
    method = "Likelihood-ratio test, Partner-only model",
    dataName = "Observed vs Estimated",
    alternative = "The unrestricted model fits the data better",
    statistic = lrtPM[["statistic"]],
    df = lrtPM[["parameter"]],
    pValue = lrtPM[["p.value"]]
  )

  out <- list(
    TEST.AM = TEST1,
    TEST.PM = TEST2,
    pattern = type,
    alpha = alpha,
    states = states,
    call = match.call()
  )
  class(out) <- c("dyadic_pattern", "list")
  out
}
