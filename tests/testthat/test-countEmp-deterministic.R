test_that("countEmp returns exact counts on a tiny hand-constructed chain", {
  chainFM <- c(1L, 2L, 1L)
  chainSM <- c(2L, 1L, 1L)

  got <- countEmp(chainFM, chainSM, states = 2L)

  # rows correspond to (bfm,bsm): (1,1),(1,2),(2,1),(2,2)
  exp <- matrix(
    c(
      0L, 0L,  # (1,1) -> next FM 1/2
      0L, 1L,  # (1,2)
      1L, 0L,  # (2,1)
      0L, 0L   # (2,2)
    ),
    nrow = 4L, ncol = 2L, byrow = TRUE
  )

  expect_identical(got, exp)
})
