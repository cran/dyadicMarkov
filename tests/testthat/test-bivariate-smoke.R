test_that("countEmpBivariate returns integer matrix with correct dimensions", {
  chainFM_V1 <- c(1L, 2L, 1L)
  chainSM_V1 <- c(1L, 1L, 2L)
  chainFM_V2 <- c(2L, 1L, 2L)
  chainSM_V2 <- c(1L, 2L, 2L)

  emp <- countEmpBivariate(chainFM_V1, chainSM_V1, chainFM_V2, chainSM_V2, states = 2L)

  expect_true(is.matrix(emp))
  expect_type(emp, "integer")
  expect_equal(dim(emp), c(16L, 2L))
})
