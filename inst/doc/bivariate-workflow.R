## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----bivariate-data-----------------------------------------------------------
utils::data("dyadic_bivariate_example", package = "dyadicMarkov")

head(dyadic_bivariate_example)
dim(dyadic_bivariate_example)

## ----bivariate-states---------------------------------------------------------
table(dyadic_bivariate_example$FM_V1)
table(dyadic_bivariate_example$SM_V1)

## ----bivariate-counts---------------------------------------------------------
emp_bi <- dyadicMarkov::countEmpBivariate(
  chainFM_V1 = dyadic_bivariate_example$FM_V1,
  chainSM_V1 = dyadic_bivariate_example$SM_V1,
  chainFM_V2 = dyadic_bivariate_example$FM_V2,
  chainSM_V2 = dyadic_bivariate_example$SM_V2,
  states = 2L
)

emp_bi
class(emp_bi)
dim(emp_bi)

## ----bivariate-case-----------------------------------------------------------
case_bi <- dyadicMarkov::bivariateCase(emp_bi, alpha = 0.05)

case_bi
case_bi$case
summary(case_bi)

## ----complete-pattern---------------------------------------------------------
complete_bi <- dyadicMarkov::completePattern(emp_bi)

complete_bi
complete_bi$pattern
complete_bi$aic

## ----bivariate-role-rotation--------------------------------------------------
analyze_bivariate <- function(label, fm_v1, sm_v1, fm_v2, sm_v2) {
  emp <- dyadicMarkov::countEmpBivariate(
    chainFM_V1 = fm_v1,
    chainSM_V1 = sm_v1,
    chainFM_V2 = fm_v2,
    chainSM_V2 = sm_v2,
    states = 2L
  )

  case <- dyadicMarkov::bivariateCase(emp, alpha = 0.05)

  cat("\n", label, "\n", sep = "")
  print(case)

  if (identical(case$case, "complete")) {
    print(dyadicMarkov::completePattern(emp))
  }

  if (identical(case$case, "partial")) {
    print(dyadicMarkov::partialPattern(emp))
  }

  if (identical(case$case, "univariate")) {
    print(dyadicMarkov::univariatePattern(fm_v1, sm_v1, states = 2L, alpha = 0.05))
  }
}

d <- dyadic_bivariate_example

analyze_bivariate(
  "FM_V1 as analyzed sequence, V1 as main variable",
  d$FM_V1, d$SM_V1, d$FM_V2, d$SM_V2
)

analyze_bivariate(
  "SM_V1 as analyzed sequence, V1 as main variable",
  d$SM_V1, d$FM_V1, d$SM_V2, d$FM_V2
)

analyze_bivariate(
  "FM_V2 as analyzed sequence, V2 as main variable",
  d$FM_V2, d$SM_V2, d$FM_V1, d$SM_V1
)

analyze_bivariate(
  "SM_V2 as analyzed sequence, V2 as main variable",
  d$SM_V2, d$FM_V2, d$SM_V1, d$FM_V1
)

