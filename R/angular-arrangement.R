# Angular arrangement -----------------------------------------------------

#' Construct the angular arrangement for centered bivariate data
#'
#' The critical angles divide the projective circle into cells in which the
#' ordering of the squared projections is constant.
#'
#' @param Y An n x 2 numeric matrix of centered observations.
#' @param tol A nonnegative numerical tolerance.
#'
#' @return An internal angular-arrangement object containing cell boundaries
#'   and projection orderings.
#'
#' @keywords internal
angular_arrangement_d2 <- function(Y, tol = 1e-12) {
  Y <- as_bivariate_matrix_d2(Y, "Y")

  if (
    !is.numeric(tol) ||
    length(tol) != 1L ||
    is.na(tol) ||
    !is.finite(tol) ||
    tol < 0
  ) {
    stop(
      "tol must be a finite nonnegative number.",
      call. = FALSE
    )
  }

  arrangement <- angular_arrangement_d2_cpp(
    Y = Y,
    tol = as.numeric(tol)
  )

  class(arrangement) <- c(
    "shd_angular_arrangement_d2",
    "list"
  )

  arrangement
}
