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
  cat("Scatter halfspace median\n")
  cat(sprintf("Method: %s (%s)\n", x$method, x$algorithm))
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
    cat(sprintf("%s: %.6f\n", if (x$method == "approx") "Finite-direction depth" else "Depth", x$depth))
  } else {
    cat("Best k: unresolved\n")
  }

  if (!is.null(x$Sigma)) {
    cat("Scatter median:\n")
    if (isTRUE(x$gaussian_consistent))
      cat(sprintf("Gaussian consistency factor: %.6f (depth refers to Sigma_raw)\n",
                  x$consistency_factor))
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
