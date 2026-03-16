test_that("countEmp returns expected shape and nonnegative counts", {
  chainFM <- rep(c(1L,2L,1L,2L,2L), length.out = 200L)
  chainSM <- rep(c(2L,1L,2L,1L,1L), length.out = 200L)

  emp <- dyadicMarkov::countEmp(chainFM = chainFM, chainSM = chainSM, states = 2L)

  expect_true(is.matrix(emp))
  expect_true(all(emp >= 0))
  expect_equal(ncol(emp), 2L)
  expect_equal(nrow(emp), 2L^2)
})

test_that("mleEstimation returns valid probabilities", {
  chainFM <- rep(c(1L,2L,1L,2L,2L), length.out = 200L)
  chainSM <- rep(c(2L,1L,2L,1L,1L), length.out = 200L)

  emp  <- dyadicMarkov::countEmp(chainFM = chainFM, chainSM = chainSM, states = 2L)
  prob <- dyadicMarkov::mleEstimation(emp)

  expect_true(is.matrix(prob))
  expect_true(all(is.finite(prob)))
  expect_true(all(prob >= 0))
  expect_true(all(prob <= 1))
  expect_true(all(abs(rowSums(prob) - 1) < 1e-10))
})

test_that("univariatePattern runs", {
  chainFM <- rep(c(1L,2L,1L,2L,2L), length.out = 200L)
  chainSM <- rep(c(2L,1L,2L,1L,1L), length.out = 200L)

  res <- dyadicMarkov::univariatePattern(
    chainFM = chainFM, chainSM = chainSM, states = 2L, alpha = 0.05
  )

  expect_true(is.list(res))
})



test_that("bivariate pipeline runs for states=2 (exported API only)", {
  states <- 2L
  T <- 200L

  chainFM_V1 <- rep(c(1L,2L,1L,2L,2L), length.out = T)
  chainSM_V1 <- rep(c(2L,1L,2L,1L,1L), length.out = T)
  chainFM_V2 <- rep(c(1L,1L,2L,2L,1L), length.out = T)
  chainSM_V2 <- rep(c(2L,2L,1L,1L,2L), length.out = T)

  emp <- dyadicMarkov::countEmpBivariate(
    chainFM_V1, chainSM_V1, chainFM_V2, chainSM_V2, states = states
  )

  expect_true(is.matrix(emp))
  expect_equal(ncol(emp), states)
  expect_equal(nrow(emp), states^4)

  res_case <- dyadicMarkov::bivariateCase(empirical = emp, alpha = 0.05)
  expect_true(is.list(res_case))

  res_partial <- dyadicMarkov::partialPattern(emp)
  expect_true(is.list(res_partial))

  res_complete <- dyadicMarkov::completePattern(emp)
  expect_true(is.list(res_complete))
})

