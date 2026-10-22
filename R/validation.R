# Input validation --------------------------------------------------------

#' Coerce an object to a numeric bivariate matrix
#'
#' @param X An object coercible to a numeric matrix.
#' @param name The object name used in error messages.
#'
#' @return A numeric matrix with exactly two columns and at least two rows.
#'
#' @keywords internal
as_bivariate_matrix_d2 <- function(X, name = "X") {
  X <- as.matrix(X)

  if (!is.numeric(X) || length(dim(X)) != 2L || ncol(X) != 2L) {
    stop(
      sprintf("%s must be a numeric matrix with exactly two columns.", name),
      call. = FALSE
    )
  }

  if (nrow(X) < 2L) {
    stop(
      sprintf("%s must contain at least two observations.", name),
      call. = FALSE
    )
  }

  if (anyNA(X) || any(!is.finite(X))) {
    stop(
      sprintf("%s must contain only finite values.", name),
      call. = FALSE
    )
  }

  storage.mode(X) <- "double"
  X
}


#' Validate a bivariate center
#'
#' @param mu A numeric center vector.
#' @param name The object name used in error messages.
#'
#' @return A finite numeric vector of length two.
#'
#' @keywords internal
as_center_vector_d2 <- function(mu, name = "mu") {
  if (!is.numeric(mu)) {
    stop(
      sprintf("%s must be numeric.", name),
      call. = FALSE
    )
  }

  mu <- as.numeric(mu)

  if (length(mu) != 2L || anyNA(mu) || any(!is.finite(mu))) {
    stop(
      sprintf("%s must be a finite numeric vector of length two.", name),
      call. = FALSE
    )
  }

  mu
}

# Shared validation for dimension-independent public entry points.
as_sample_matrix <- function(X) {
  X <- as.matrix(X)
  if (!is.numeric(X) || is.complex(X) || length(dim(X)) != 2L ||
      any(dim(X) == 0L) || any(!is.finite(X)))
    stop("X must be a nonempty finite numeric matrix.", call. = FALSE)
  storage.mode(X) <- "double"
  X
}
