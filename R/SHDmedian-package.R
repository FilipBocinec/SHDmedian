#' SHDmedian: Scatter Halfspace Medians
#'
#' The package computes the scatter halfspace median for bivariate data using
#' angular sweeps and a point-based second-order cone formulation.
#' It also provides finite-direction SDP approximations in arbitrary dimension
#' through the unified SHDmedian() interface, with a supplied, mean, or Tukey
#' center shared by both computational methods.
#'
#' @keywords internal
#' @useDynLib SHDmedian, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
