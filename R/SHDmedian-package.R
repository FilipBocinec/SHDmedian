#' SHDmedian: Scatter Halfspace Medians
#'
#' The package computes the scatter halfspace median for bivariate data using
#' an angular arrangement and a cutting-plane algorithm.
#' It also provides finite-direction SDP approximations in arbitrary dimension
#' through the unified SHDmedian() interface, with a supplied, mean, or Tukey
#' center shared by both computational methods.
#'
#' @keywords internal
#' @useDynLib SHDmedian, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
