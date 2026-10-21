#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

// Internal kernels receive validated, normalized data from the R wrapper.
// [[Rcpp::export]]
Rcpp::List shd_prepare_cpp(const arma::mat& Y, const arma::mat& U) {
  if (Y.n_cols != U.n_rows || Y.n_rows == 0 || U.n_cols == 0 ||
      !Y.is_finite() || !U.is_finite()) Rcpp::stop("Invalid kernel dimensions or data.");
  arma::mat Q = arma::square(Y * U);
  if (!Q.is_finite()) Rcpp::stop("Squared projections overflowed.");
  for (arma::uword j = 0; j < Q.n_cols; ++j) {
    Rcpp::checkUserInterrupt();
    std::sort(Q.colptr(j), Q.colptr(j) + Q.n_rows);
  }
  const arma::uword d = U.n_rows;
  arma::mat C(U.n_cols, d * (d + 1) / 2);
  arma::uword index = 0;
  for (arma::uword col = 0; col < d; ++col)
    for (arma::uword row = 0; row <= col; ++row, ++index)
      C.col(index) = (row == col ? 1.0 : std::sqrt(2.0)) *
        (U.row(row) % U.row(col)).t();
  return Rcpp::List::create(Rcpp::Named("Q") = Q, Rcpp::Named("C") = C);
}

// [[Rcpp::export]]
Rcpp::List shd_verify_cpp(const arma::mat& Sigma, const arma::mat& U,
                         const arma::mat& Q, const double tolerance = 0) {
  if (Sigma.n_rows != Sigma.n_cols || Sigma.n_rows != U.n_rows ||
      Q.n_cols != U.n_cols || Q.n_rows == 0 || U.n_cols == 0 ||
      !Sigma.is_finite() || !U.is_finite() || !Q.is_finite() ||
      !std::isfinite(tolerance) || tolerance < 0) Rcpp::stop("Invalid verification input.");
  arma::rowvec values = arma::sum(U % (Sigma * U), 0);
  const arma::uword L = U.n_cols;
  Rcpp::NumericVector lower(L), upper(L), expanded_lower(L), expanded_upper(L);
  for (arma::uword j = 0; j < L; ++j) {
    const double* begin = Q.colptr(j);
    const double* end = begin + Q.n_rows;
    lower[j] = std::upper_bound(begin, end, values[j]) - begin;
    upper[j] = end - std::lower_bound(begin, end, values[j]);
    expanded_lower[j] = std::upper_bound(begin, end, values[j] + tolerance) - begin;
    expanded_upper[j] = end - std::lower_bound(begin, end, values[j] - tolerance);
  }
  arma::vec eigenvalues;
  if (!arma::eig_sym(eigenvalues, Sigma)) Rcpp::stop("Eigenvalue calculation failed.");
  return Rcpp::List::create(Rcpp::Named("values") = Rcpp::NumericVector(values.begin(), values.end()),
    Rcpp::Named("lower") = lower, Rcpp::Named("upper") = upper,
    Rcpp::Named("expanded_lower") = expanded_lower,
    Rcpp::Named("expanded_upper") = expanded_upper,
    Rcpp::Named("min_eigenvalue") = eigenvalues.min());
}
