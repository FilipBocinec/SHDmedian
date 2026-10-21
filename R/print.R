# Print methods -----------------------------------------------------------

#' Print a scatter halfspace median result
#'
#' @param x An object returned by [SHDmedian()].
#' @param ... Additional arguments retained for S3 compatibility.
#'
#' @return `x`, invisibly.
#'
#' @export
print.SHDmedian <- function(x, ...) {
  cat("Bivariate scatter halfspace median\n")
  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("Status: %s\n", x$status))
  cat(sprintf(
    "Certified: %s\n",
    if (isTRUE(x$certified)) "yes" else "no"
  ))
  cat(sprintf("Sample size: %d\n", x$n))
  cat(sprintf(
    "Search interval: [%d, %d]\n",
    x$k_min,
    x$k_max
  ))

  cat("Center:\n")
  print(x$center)

  if (!is.na(x$best_k)) {
    cat(sprintf("Best k: %d\n", x$best_k))
    cat(sprintf("Depth: %.6f\n", x$depth))
  } else {
    cat("Best k: unresolved\n")
  }

  if (!is.null(x$Sigma)) {
    cat("Scatter median:\n")
    print(x$Sigma)
  } else {
    cat("Scatter median: unavailable\n")
  }

  cat(sprintf(
    "Elapsed time: %.3f seconds\n",
    x$timing$total
  ))

  invisible(x)
}
