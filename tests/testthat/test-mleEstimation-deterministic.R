test_that("mleEstimation normalizes rows and uses uniform for zero rows", {
  empirical <- matrix(
    c(
      1, 1,  # sum 2 -> (0.5,0.5)
      0, 0,  # sum 0 -> uniform
      0, 2,  # sum 2 -> (0,1)
      3, 1   # sum 4 -> (0.75,0.25)
    ),
    nrow = 4L, ncol = 2L, byrow = TRUE
  )

  got <- mleEstimation(empirical)

  exp <- matrix(
    c(
      0.5,  0.5,
      0.5,  0.5,
      0.0,  1.0,
      0.75, 0.25
    ),
    nrow = 4L, ncol = 2L, byrow = TRUE
  )

  expect_identical(unname(unclass(got)), exp)
  expect_identical(rownames(got), c(
    "previous_1", "previous_2", "previous_3", "previous_4"
  ))
  expect_identical(colnames(got), c("next_1", "next_2"))
})
