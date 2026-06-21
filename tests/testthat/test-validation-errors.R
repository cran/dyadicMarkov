#' @srrstats {G5.2} Shipped validation tests explicitly trigger representative error conditions for exported univariate, bivariate, empirical-matrix, and pattern-identification functions.
#' @srrstats {G5.2a} Validation errors use explicit, context-specific messages; shipped validation tests compare representative user-facing errors against expected messages.
#' @srrstats {G5.2b} Shipped tests use expect_error(..., fixed = TRUE) to compare validation failures against expected error messages for unsupported state spaces, invalid chain values, missing values, length mismatches, invalid empirical matrices, and invalid alpha values.
#' @srrstats {G5.8} Shipped tests cover edge and validation conditions, including zero-length chains, unsupported input types, missing values, invalid state values, length mismatches, invalid empirical matrices, and valid all-identical chains.
#' @srrstats {G5.8a} Zero-length chain inputs are tested and produce the expected explicit error.
#' @srrstats {G5.8b} Unsupported character chain inputs are tested and produce the expected explicit error.
#' @srrstats {G5.8c} Missing-value chains are tested and rejected explicitly; all-identical valid chains are tested to produce expected counts and probabilities.
#' @srrstats {G5.8d} Inputs outside the supported algorithmic scope are tested, including unsupported bivariate states, invalid state values, invalid empirical dimensions, and invalid alpha values.

test_that("univariate input validation errors are explicit", {
  expect_error(
    dyadicMarkov::countEmp(numeric(0), numeric(0), states = 2L),
    "chains must have length >= 2.",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::countEmp(c("1", "2"), c("1", "2"), states = 2L),
    "chain values must be numeric vectors of integer-coded states.",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::countEmp(c(1L, NA_integer_), c(1L, 2L), states = 2L),
    "chains must not contain NA.",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::countEmp(c(1L, 2L, 1L), c(1L, 2L), states = 2L),
    "chainFM and chainSM must have the same length.",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::countEmp(c(1L, 3L), c(1L, 2L), states = 2L),
    "chain values must be integers in 1:states.",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::countEmp(c(1L, 2L), c(1L, 2L), states = 1L),
    "states must be provided as a single integer >= 2.",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::countEmp(c(1L, 2L), c(1L, 2L), states = 1e10),
    "states must be provided as a single integer >= 2.",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::countEmp(c(1, Inf), c(1L, 2L), states = 2L),
    "chains must contain finite integer-coded states.",
    fixed = TRUE
  )
})


test_that("empirical matrix validation errors are explicit", {
  expect_error(
    dyadicMarkov::mleEstimation(c(1, 2, 3)),
    "empirical must be a matrix.",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::mleEstimation(matrix(numeric(0), nrow = 2L, ncol = 0L)),
    "empirical must have at least one row and at least two columns.",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::mleEstimation(matrix(c(1, NA, 2, 3), nrow = 2L)),
    "empirical must contain finite non-negative counts with no NA.",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::mleEstimation(matrix(c(1, -1, 2, 3), nrow = 2L)),
    "empirical must contain finite non-negative counts with no NA.",
    fixed = TRUE
  )
})


test_that("bivariate validation errors are explicit", {
  expect_error(
    dyadicMarkov::countEmpBivariate(
      c(1L, 2L), c(1L, 2L), c(1L, 2L), c(1L, 2L), states = 3L
    ),
    "bivariate functions currently support states = 2 only.",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::countEmpBivariate(
      c(1L, 2L), c(1L, 2L), c(1L, 2L), c(1L, 2L), states = 1e10
    ),
    "bivariate functions currently support states = 2 only.",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::countEmpBivariate(
      c("1", "2"), c(1L, 2L), c(1L, 2L), c(1L, 2L),
      states = 2L
    ),
    "chain values must be numeric integers in {1, 2}.",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::countEmpBivariate(
      c(1L, Inf), c(1L, 2L), c(1L, 2L), c(1L, 2L), states = 2L
    ),
    "bivariate chains must contain finite integer-coded states.",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::countEmpBivariate(
      c(1L, 2L), c(1L, 2L), c(1L, 2L), c(1L, 2L, 1L), states = 2L
    ),
    "chains must have the same length.",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::countEmpBivariate(
      c(1L, NA_integer_), c(1L, 2L), c(1L, 2L), c(1L, 2L), states = 2L
    ),
    "bivariate chains must not contain NA.",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::countEmpBivariate(
      c(1L, 3L), c(1L, 2L), c(1L, 2L), c(1L, 2L), states = 2L
    ),
    "with states = 2, chain values must be integers in {1, 2}.",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::bivariateCase(matrix(1, nrow = 4L, ncol = 2L), alpha = 0.05),
    "bivariate functions currently support states = 2 only (empirical must be a 16x2 matrix).",
    fixed = TRUE
  )
})


test_that("all-identical valid chains produce expected counts and probabilities", {
  chainFM <- rep(1L, 5L)
  chainSM <- rep(1L, 5L)

  counts <- dyadicMarkov::countEmp(chainFM, chainSM, states = 2L)

  expected_counts <- matrix(0L, nrow = 4L, ncol = 2L)
  expected_counts[1L, 1L] <- 4L

  expect_identical(unname(unclass(counts)), expected_counts)
  expect_false(is.null(dimnames(counts)))

  probs <- dyadicMarkov::mleEstimation(counts)

  expected_probs <- matrix(
    c(
      1.0, 0.0,
      0.5, 0.5,
      0.5, 0.5,
      0.5, 0.5
    ),
    nrow = 4L, ncol = 2L, byrow = TRUE
  )

  expect_equal(unname(unclass(probs)), expected_probs)
  expect_false(is.null(dimnames(probs)))
  expect_true(all(is.finite(probs)))
  expect_true(all(abs(rowSums(probs) - 1) < 1e-10))
})


test_that("pattern functions validate alpha and bivariate empirical inputs explicitly", {
  chainFM <- c(1L, 2L, 1L)
  chainSM <- c(2L, 1L, 1L)

  expect_error(
    dyadicMarkov::univariatePattern(chainFM, chainSM, states = 2L, alpha = 1),
    "alpha must be a single number in (0, 1).",
    fixed = TRUE
  )

  empirical_bad <- matrix(1, nrow = 4L, ncol = 2L)

  expect_error(
    dyadicMarkov::partialPattern(empirical_bad),
    "bivariate functions currently support states = 2 only (empirical must be a 16x2 matrix).",
    fixed = TRUE
  )

  expect_error(
    dyadicMarkov::completePattern(empirical_bad),
    "bivariate functions currently support states = 2 only (empirical must be a 16x2 matrix).",
    fixed = TRUE
  )
})
