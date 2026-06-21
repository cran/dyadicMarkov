test_that("mleEstimation recovers known transition probabilities within tolerance", {
  empirical <- matrix(
    c(
      30, 70,
      25, 75,
      60, 40,
      10, 90
    ),
    nrow = 4L, ncol = 2L, byrow = TRUE
  )

  expected <- matrix(
    c(
      0.30, 0.70,
      0.25, 0.75,
      0.60, 0.40,
      0.10, 0.90
    ),
    nrow = 4L, ncol = 2L, byrow = TRUE
  )

  recovered <- dyadicMarkov::mleEstimation(empirical)

  expect_equal(unname(unclass(recovered)), expected, tolerance = 1e-12)
  expect_true(all(is.finite(recovered)))
  expect_true(all(abs(rowSums(recovered) - 1) < 1e-12))
})


test_that("mleEstimation is stable when empirical counts are scaled", {
  empirical <- matrix(
    c(
      30, 70,
      25, 75,
      60, 40,
      10, 90
    ),
    nrow = 4L, ncol = 2L, byrow = TRUE
  )

  recovered_small <- dyadicMarkov::mleEstimation(empirical)
  recovered_large <- dyadicMarkov::mleEstimation(10L * empirical)

  expect_equal(
    unclass(recovered_large),
    unclass(recovered_small),
    tolerance = 1e-12
  )
})
