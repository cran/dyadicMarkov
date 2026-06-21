test_that("countEmpBivariate returns exact counts on a tiny hand-constructed case", {
  chainFM_V1 <- c(1L, 2L, 1L)
  chainSM_V1 <- c(1L, 1L, 2L)
  chainFM_V2 <- c(2L, 1L, 2L)
  chainSM_V2 <- c(1L, 2L, 2L)

  got <- countEmpBivariate(chainFM_V1, chainSM_V1, chainFM_V2, chainSM_V2, states = 2L)

  exp <- matrix(0L, nrow = 16L, ncol = 2L)

  # Transition 1: row=3, col=2 gets +1
  exp[3L, 2L] <- 1L
  # Transition 2: row=10, col=1 gets +1
  exp[10L, 1L] <- 1L

  expect_identical(unname(unclass(got)), exp)
  expect_identical(length(rownames(got)), 16L)
  expect_identical(rownames(got)[1L], "mainFM1_mainSM1_secondFM1_secondSM1")
  expect_identical(rownames(got)[10L], "mainFM2_mainSM1_secondFM1_secondSM2")
  expect_identical(colnames(got), c("next_1", "next_2"))
})
