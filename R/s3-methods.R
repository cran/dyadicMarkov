#' Print a line with newline termination
#'
#' Internal helper used by print methods.
#'
#' @param ... Objects passed to `cat()`.
#'
#' @return Invisibly returns `NULL`.
#' @noRd
#' @srrstats {EA5.3} Summary methods for matrix-like dyadic objects report object class, storage mode, dimensions, and row and column names.
#' @srrstats {EA5.2} Matrix-like dyadic print methods format numeric output with formatC() rather than relying on default numeric printing.
#' @srrstats {EA4.2} Primary returned pattern and case objects implement print and summary methods, while empirical counts and MLE probabilities preserve matrix-like behaviour.
#' @srrstats {EA4.1} Matrix-like dyadic print methods provide a digits argument to explicitly control numeric precision in screen output.
#' @srrstats {EA4.0} Return objects use documented S3 classes, including dyadic_counts, dyadic_mle, dyadic_pattern, and dyadic_case, with matrix-like or list-like structure matching the analytic output.
.cat_line <- function(...) {
  cat(..., "\n", sep = "")
  invisible(NULL)
}


#' @exportS3Method base::print
#' @noRd
print.dyadic_pattern <- function(x, ...) {
  .cat_line("Dyadic interaction pattern")

  if ("pattern" %in% names(x)) {
    .cat_line("Pattern: ", paste0(x[["pattern"]], collapse = ", "))
  }
  if ("alpha" %in% names(x)) {
    .cat_line("Alpha: ", paste0(x[["alpha"]], collapse = ", "))
  }
  if ("states" %in% names(x)) {
    .cat_line("States: ", paste0(x[["states"]], collapse = ", "))
  }

  invisible(x)
}


#' @exportS3Method base::print
#' @noRd
print.dyadic_case <- function(x, ...) {
  .cat_line("Bivariate dyadic case")

  if ("case" %in% names(x)) {
    .cat_line("Case: ", paste0(x[["case"]], collapse = ", "))
  }
  if ("alpha" %in% names(x)) {
    .cat_line("Alpha: ", paste0(x[["alpha"]], collapse = ", "))
  }

  invisible(x)
}


#' @exportS3Method base::summary
#' @noRd
summary.dyadic_pattern <- function(object, ...) {
  fields <- c("pattern", "aic", "alpha", "states", "call")
  out <- object[intersect(fields, names(object))]
  class(out) <- c("summary_dyadic_pattern", "list")
  out
}


#' @exportS3Method base::summary
#' @noRd
summary.dyadic_case <- function(object, ...) {
  fields <- c("case", "alpha", "call")
  out <- object[intersect(fields, names(object))]
  class(out) <- c("summary_dyadic_case", "list")
  out
}


#' Format a dyadic matrix for screen output
#'
#' Internal helper used by matrix-like print methods.
#'
#' @param x Matrix-like dyadic object.
#' @param digits Number of digits to display.
#'
#' @return A character matrix preserving dimnames.
#' @noRd
.format_dyadic_matrix <- function(x, digits) {
  digits <- as.integer(digits)[1L]
  if (is.na(digits) || digits < 0L) {
    digits <- 0L
  }

  out <- formatC(unclass(x), format = "f", digits = digits)
  dim(out) <- dim(x)
  dimnames(out) <- dimnames(x)
  out
}


#' Build a summary for matrix-like dyadic objects
#'
#' Internal helper used by summary methods.
#'
#' @param object Matrix-like dyadic object.
#' @param object_type Character description of the object type.
#'
#' @return A list with class, type, dimension, and name metadata.
#' @noRd
.summary_dyadic_matrix <- function(object, object_type) {
  list(
    object_type = object_type,
    object_class = class(object),
    storage_mode = storage.mode(object),
    dimensions = dim(object),
    row_names = rownames(object),
    column_names = colnames(object)
  )
}


#' @exportS3Method base::print
#' @noRd
print.dyadic_counts <- function(x, digits = 0L, ...) {
  print(.format_dyadic_matrix(x, digits = digits), quote = FALSE, ...)
  invisible(x)
}


#' @exportS3Method base::print
#' @noRd
print.dyadic_mle <- function(x, digits = 3L, ...) {
  print(.format_dyadic_matrix(x, digits = digits), quote = FALSE, ...)
  invisible(x)
}


#' @exportS3Method base::summary
#' @noRd
summary.dyadic_counts <- function(object, ...) {
  out <- .summary_dyadic_matrix(
    object = object,
    object_type = "empirical transition counts"
  )
  out$total_count <- sum(object)
  out$row_totals <- rowSums(object)
  class(out) <- c("summary_dyadic_counts", "list")
  out
}


#' @exportS3Method base::summary
#' @noRd
summary.dyadic_mle <- function(object, ...) {
  out <- .summary_dyadic_matrix(
    object = object,
    object_type = "transition probability estimates"
  )
  out$row_totals <- rowSums(object)
  class(out) <- c("summary_dyadic_mle", "list")
  out
}
