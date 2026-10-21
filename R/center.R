# Data centering ----------------------------------------------------------

#' Center a bivariate sample
#'
#' The Tukey median is computed using
#' [mrfDepth::hdepthmedian()] with `maxdir = NULL`.
#'
#' @param X An n x 2 numeric matrix.
#' @param method Either `"tukey"` or `"mean"`.
#' @param mu An optional user-supplied center. When supplied, `method` is
#'   ignored.
#'
#' @return A list containing the original observations, the selected center,
#'   the centered observations, the effective centering method, and the
#'   Tukey depth of the center when available.
#'
#' @keywords internal
center_bivariate_data_d2 <- function(
    X,
    method = c("tukey", "mean"),
    mu = NULL) {
  X <- as_bivariate_matrix_d2(X, "X")
  method <- match.arg(method)

  center_depth <- NA_real_

  if (!is.null(mu)) {
    mu_hat <- as_center_vector_d2(mu, "mu")
    effective_method <- "supplied"
  } else if (method == "mean") {
    mu_hat <- colMeans(X)
    effective_method <- "mean"
  } else {
    median_result <- mrfDepth::hdepthmedian(
      x = X,
      maxdir = NULL
    )

    if (is.null(median_result$median)) {
      stop(
        "mrfDepth::hdepthmedian() did not return a median.",
        call. = FALSE
      )
    }

    mu_hat <- as_center_vector_d2(
      median_result$median,
      "computed Tukey median"
    )

    if (!is.null(median_result$depth)) {
      center_depth <- as.numeric(median_result$depth)[1L]
    }

    effective_method <- "tukey"
  }

  Y <- sweep(
    X,
    MARGIN = 2L,
    STATS = mu_hat,
    FUN = "-"
  )

  list(
    X = X,
    mu = mu_hat,
    Y = Y,
    method = effective_method,
    center_depth = center_depth
  )
}
