## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
chainFM <- c(
  1L,2L,1L,2L,2L,1L,  1L,1L,2L,2L,1L,2L,
  2L,1L,1L,2L,1L,2L,  1L,2L,2L,1L,2L,1L
)
chainSM <- c(
  2L,1L,2L,1L,1L,2L,  2L,2L,1L,1L,2L,1L,
  1L,2L,2L,1L,2L,1L,  2L,1L,1L,2L,1L,2L
)

length(chainFM)
head(chainFM)
length(chainSM)
head(chainSM)

## -----------------------------------------------------------------------------
states <- 2L
emp <- dyadicMarkov::countEmp(chainFM = chainFM, chainSM = chainSM, states = states)
emp

## -----------------------------------------------------------------------------
fit <- dyadicMarkov::mleEstimation(emp)
fit
rowSums(fit)

## -----------------------------------------------------------------------------
pat <- dyadicMarkov::univariatePattern(
  chainFM = chainFM,
  chainSM = chainSM,
  states  = 2L,
  alpha   = 0.05
)

pat[["pattern"]]
pat[["TEST.AM"]]
pat[["TEST.PM"]]

## -----------------------------------------------------------------------------
chainFM_V1 <- chainFM
chainSM_V1 <- chainSM

chainFM_V2 <- c(
  1L,1L,2L,2L,1L,2L,  2L,1L,1L,2L,1L,2L,
  1L,2L,1L,2L,2L,1L,  2L,2L,1L,1L,2L,1L
)
chainSM_V2 <- c(
  2L,2L,1L,1L,2L,1L,  1L,2L,2L,1L,2L,1L,
  2L,1L,2L,1L,1L,2L,  1L,1L,2L,2L,1L,2L
)

emp2 <- dyadicMarkov::countEmpBivariate(
  chainFM_V1, chainSM_V1,
  chainFM_V2, chainSM_V2,
  states = 2L
)

emp2

dyadicMarkov::bivariateCase(emp2, alpha = 0.05)

