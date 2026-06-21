test_that("pattern and case functions return S3 list objects", {
  chainFM <- c(1L, 2L, 1L, 2L, 2L, 1L)
  chainSM <- c(2L, 1L, 2L, 1L, 1L, 2L)

  uni <- dyadicMarkov::univariatePattern(chainFM, chainSM, states = 2L, alpha = 0.05)

  expect_s3_class(uni, "dyadic_pattern")
  expect_true("pattern" %in% names(uni))
  expect_true(is.character(uni$pattern) || is.na(uni$pattern))

  chainFM_V1 <- c(1L, 2L, 1L, 2L, 2L, 1L)
  chainSM_V1 <- c(2L, 1L, 2L, 1L, 1L, 2L)
  chainFM_V2 <- c(1L, 1L, 2L, 2L, 1L, 2L)
  chainSM_V2 <- c(2L, 2L, 1L, 1L, 2L, 1L)

  emp <- dyadicMarkov::countEmpBivariate(
    chainFM_V1, chainSM_V1, chainFM_V2, chainSM_V2, states = 2L
  )

  case <- dyadicMarkov::bivariateCase(emp, alpha = 0.05)
  partial <- dyadicMarkov::partialPattern(emp)
  complete <- dyadicMarkov::completePattern(emp)

  expect_s3_class(case, "dyadic_case")
  expect_true("case" %in% names(case))
  expect_true(is.character(case$case) || is.na(case$case))

  expect_s3_class(partial, "dyadic_pattern")
  expect_true("pattern" %in% names(partial))
  expect_type(partial$pattern, "character")

  expect_s3_class(complete, "dyadic_pattern")
  expect_true("pattern" %in% names(complete))
  expect_type(complete$pattern, "character")

  uni_summary <- summary(uni)
  expect_s3_class(uni_summary, "summary_dyadic_pattern")
  expect_true(all(c("pattern", "alpha", "states", "call") %in% names(uni_summary)))

  partial_summary <- summary(partial)
  expect_s3_class(partial_summary, "summary_dyadic_pattern")
  expect_true(all(c("pattern", "aic", "call") %in% names(partial_summary)))

  complete_summary <- summary(complete)
  expect_s3_class(complete_summary, "summary_dyadic_pattern")
  expect_true(all(c("pattern", "aic", "call") %in% names(complete_summary)))

  case_summary <- summary(case)
  expect_s3_class(case_summary, "summary_dyadic_case")
  expect_true(all(c("case", "alpha", "call") %in% names(case_summary)))

  expect_output(print(uni), "Dyadic interaction pattern")
  expect_output(print(case), "Bivariate dyadic case")
  expect_output(print(partial), "Dyadic interaction pattern")
  expect_output(print(complete), "Dyadic interaction pattern")
})
