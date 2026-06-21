test_that("mleEstimation returns a dyadic_mle matrix with unchanged values", {
#' @srrstats {EA6.0e} Deterministic tests compare exact count matrices and known transition probability values using explicit expected values and numerical tolerances.
#' @srrstats {EA6.0c} Tests explicitly check equivalent row and column names for empirical count and MLE matrix outputs through their matrix dimnames.
#' @srrstats {EA6.0b} Tests explicitly check the dimensions of empirical transition count matrices and exported workflow outputs.
#' @srrstats {EA6.0a} Tests explicitly check S3 classes and object types returned by pattern, case, empirical count, and MLE functions.
#' @srrstats {EA6.0} Tests cover return values from exported workflows, including classes, dimensions, finite values, valid probabilities, and known numerical outputs.
  empirical <- matrix(
    c(
      1, 1,
      0, 0,
      0, 2,
      3, 1
    ),
    nrow = 4L, ncol = 2L, byrow = TRUE
  )

  fit <- dyadicMarkov::mleEstimation(empirical)
  expected <- matrix(
    c(
      0.5,  0.5,
      0.5,  0.5,
      0.0,  1.0,
      0.75, 0.25
    ),
    nrow = 4L, ncol = 2L, byrow = TRUE
  )

  expect_s3_class(fit, "dyadic_mle")
  expect_true(is.matrix(fit))
  expect_equal(unname(unclass(fit)), expected)
  expect_false(is.null(dimnames(fit)))
  expect_true(all(abs(rowSums(fit) - 1) < 1e-10))
})


test_that("empirical count functions return dyadic_counts matrices", {
  chainFM <- c(1L, 2L, 1L)
  chainSM <- c(2L, 1L, 1L)

  counts <- dyadicMarkov::countEmp(chainFM, chainSM, states = 2L)
  expected_counts <- matrix(
    c(
      0L, 0L,
      0L, 1L,
      1L, 0L,
      0L, 0L
    ),
    nrow = 4L, ncol = 2L, byrow = TRUE
  )

  expect_s3_class(counts, "dyadic_counts")
  expect_true(is.matrix(counts))
  expect_identical(unname(unclass(counts)), expected_counts)
  expect_identical(colnames(counts), c("next_1", "next_2"))
  expect_identical(dim(counts), c(4L, 2L))

  chainFM_V1 <- c(1L, 2L, 1L)
  chainSM_V1 <- c(1L, 1L, 2L)
  chainFM_V2 <- c(2L, 1L, 2L)
  chainSM_V2 <- c(1L, 2L, 2L)

  bivar_counts <- dyadicMarkov::countEmpBivariate(
    chainFM_V1, chainSM_V1, chainFM_V2, chainSM_V2, states = 2L
  )
  expected_bivar_counts <- matrix(0L, nrow = 16L, ncol = 2L)
  expected_bivar_counts[3L, 2L] <- 1L
  expected_bivar_counts[10L, 1L] <- 1L

  expect_s3_class(bivar_counts, "dyadic_counts")
  expect_true(is.matrix(bivar_counts))
  expect_identical(unname(unclass(bivar_counts)), expected_bivar_counts)
  expect_identical(colnames(bivar_counts), c("next_1", "next_2"))
  expect_identical(length(rownames(bivar_counts)), 16L)
  expect_identical(dim(bivar_counts), c(16L, 2L))
})


test_that("matrix-like dyadic objects have formatted print and summary methods", {
  counts <- dyadicMarkov::countEmp(
    chainFM = c(1L, 1L, 2L, 2L),
    chainSM = c(1L, 2L, 1L, 2L),
    states = 2L
  )
  probs <- dyadicMarkov::mleEstimation(counts)

  counts_print <- capture.output(print(counts, digits = 0L))
  probs_print <- capture.output(print(probs, digits = 2L))

  expect_true(any(grepl("next_1", counts_print, fixed = TRUE)))
  expect_true(any(grepl("FM1_SM1", counts_print, fixed = TRUE)))
  expect_true(any(grepl("next_1", probs_print, fixed = TRUE)))
  expect_true(any(grepl("0.50", probs_print, fixed = TRUE)))

  counts_summary <- summary(counts)
  probs_summary <- summary(probs)

  expect_s3_class(counts_summary, "summary_dyadic_counts")
  expect_s3_class(probs_summary, "summary_dyadic_mle")

  expect_identical(counts_summary$object_class, class(counts))
  expect_identical(probs_summary$object_class, class(probs))
  expect_identical(counts_summary$storage_mode, storage.mode(counts))
  expect_identical(probs_summary$storage_mode, storage.mode(probs))
  expect_identical(counts_summary$dimensions, dim(counts))
  expect_identical(probs_summary$dimensions, dim(probs))
  expect_identical(counts_summary$column_names, colnames(counts))
  expect_identical(probs_summary$column_names, colnames(probs))
})
