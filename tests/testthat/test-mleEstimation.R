test_that("mleEstimation runs on countEmp output and returns numeric matrix", {
  set.seed(888)
  chainFM <- sample(1:2, 60, TRUE)
  chainSM <- sample(1:2, 60, TRUE)

  emp <- countEmp(chainFM, chainSM, states = 2L)
  fit <- mleEstimation(emp)

  expect_true(is.matrix(fit))
  expect_type(fit, "double")
})
