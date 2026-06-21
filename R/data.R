#' Synthetic univariate dyadic sequence example
#'
#' A synthetic dyadic sequence with 90 observations, designed for package
#' workflow examples.
#'
#' @format A data frame with 90 rows and 3 columns:
#' \describe{
#'   \item{time}{Index of the measurement occasion.}
#'   \item{FM}{Integer state for the first member, taking values 1 or 2.}
#'   \item{SM}{Integer state for the second member, taking values 1 or 2.}
#' }
#' @details The package workflow classifies this example as \code{PM (A3)}
#'   using \code{\link{univariatePattern}} with \code{states = 2}.
#' @docType data
#' @name dyadic_univariate_example
"dyadic_univariate_example"


#' Synthetic bivariate dyadic sequence example
#'
#' A synthetic bivariate dyadic sequence with 90 observations, designed for
#' package workflow examples.
#'
#' @format A data frame with 90 rows and 5 columns:
#' \describe{
#'   \item{time}{Index of the measurement occasion.}
#'   \item{FM_V1}{Integer variable 1 state for the first member, taking values
#'     1 or 2.}
#'   \item{SM_V1}{Integer variable 1 state for the second member, taking values
#'     1 or 2.}
#'   \item{FM_V2}{Integer variable 2 state for the first member, taking values
#'     1 or 2.}
#'   \item{SM_V2}{Integer variable 2 state for the second member, taking values
#'     1 or 2.}
#' }
#' @details The bivariate workflow classifies this example as \code{complete}
#'   using \code{\link{bivariateCase}} with \code{alpha = 0.05}. The complete
#'   bivariate pattern selected by \code{\link{completePattern}} is \code{D2}.
#' @docType data
#' @name dyadic_bivariate_example
"dyadic_bivariate_example"
