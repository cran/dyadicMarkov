## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----univariate-data----------------------------------------------------------
utils::data("dyadic_univariate_example", package = "dyadicMarkov")

head(dyadic_univariate_example)
dim(dyadic_univariate_example)

## ----univariate-states--------------------------------------------------------
table(dyadic_univariate_example$FM)
table(dyadic_univariate_example$SM)

## ----univariate-counts--------------------------------------------------------
emp_uni <- dyadicMarkov::countEmp(
  chainFM = dyadic_univariate_example$FM,
  chainSM = dyadic_univariate_example$SM,
  states = 2L
)

emp_uni
class(emp_uni)

## ----univariate-counts-matrix-------------------------------------------------
dim(emp_uni)
rowSums(emp_uni)

## ----univariate-mle-----------------------------------------------------------
fit_uni <- dyadicMarkov::mleEstimation(emp_uni)

round(fit_uni, 3)
rowSums(fit_uni)
class(fit_uni)

## ----univariate-pattern-------------------------------------------------------
pat_uni <- dyadicMarkov::univariatePattern(
  chainFM = dyadic_univariate_example$FM,
  chainSM = dyadic_univariate_example$SM,
  states = 2L,
  alpha = 0.05
)

pat_uni

## ----univariate-pattern-details-----------------------------------------------
pat_uni$pattern
pat_uni$TEST.AM
pat_uni$TEST.PM
summary(pat_uni)

## ----reverse-perspective------------------------------------------------------
pat_uni_reverse <- dyadicMarkov::univariatePattern(
  chainFM = dyadic_univariate_example$SM,
  chainSM = dyadic_univariate_example$FM,
  states = 2L,
  alpha = 0.05
)

pat_uni_reverse

