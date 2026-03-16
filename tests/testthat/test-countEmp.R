test_that("countEmp returns integer matrix with correct dimensions", {
  chainFM <- c(1L, 2L, 1L)
  chainSM <- c(2L, 1L, 1L)

  out <- countEmp(chainFM = chainFM, chainSM = chainSM, states = 2L)

  expect_true(is.matrix(out))
  expect_type(out, "integer")
  expect_equal(dim(out), c(4L, 2L))
})
