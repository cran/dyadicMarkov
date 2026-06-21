test_that("example datasets have expected shape and state values", {
  utils::data("dyadic_univariate_example", package = "dyadicMarkov", envir = environment())
  utils::data("dyadic_bivariate_example", package = "dyadicMarkov", envir = environment())

  expect_identical(nrow(dyadic_univariate_example), 90L)
  expect_identical(names(dyadic_univariate_example), c("time", "FM", "SM"))

  expect_identical(nrow(dyadic_bivariate_example), 90L)
  expect_identical(
    names(dyadic_bivariate_example),
    c("time", "FM_V1", "SM_V1", "FM_V2", "SM_V2")
  )

  univariate_states <- dyadic_univariate_example[c("FM", "SM")]
  bivariate_states <- dyadic_bivariate_example[c("FM_V1", "SM_V1", "FM_V2", "SM_V2")]

  expect_true(all(vapply(univariate_states, function(x) {
    all(x == as.integer(x) & x %in% c(1L, 2L))
  }, logical(1))))

  expect_true(all(vapply(bivariate_states, function(x) {
    all(x == as.integer(x) & x %in% c(1L, 2L))
  }, logical(1))))
})


test_that("example datasets reproduce expected workflow classifications", {
  utils::data("dyadic_univariate_example", package = "dyadicMarkov", envir = environment())
  utils::data("dyadic_bivariate_example", package = "dyadicMarkov", envir = environment())

  univariate_result <- dyadicMarkov::univariatePattern(
    dyadic_univariate_example$FM,
    dyadic_univariate_example$SM,
    states = 2L,
    alpha = 0.05
  )
  expect_identical(univariate_result$pattern, "PM (A3)")

  bivariate_counts <- dyadicMarkov::countEmpBivariate(
    chainFM_V1 = dyadic_bivariate_example$FM_V1,
    chainSM_V1 = dyadic_bivariate_example$SM_V1,
    chainFM_V2 = dyadic_bivariate_example$FM_V2,
    chainSM_V2 = dyadic_bivariate_example$SM_V2,
    states = 2L
  )

  bivariate_case <- dyadicMarkov::bivariateCase(bivariate_counts, alpha = 0.05)
  expect_identical(bivariate_case$case, "complete")

  complete_pattern <- dyadicMarkov::completePattern(bivariate_counts)
  expect_identical(
    complete_pattern$pattern,
    "complete actor on the main, actor partner on the second (D2)"
  )
})
