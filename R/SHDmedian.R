#' Compute a scatter halfspace median
#'
#' Unified interface to finite-direction SDP optimization in any dimension and
#' the exact bivariate cutting-plane algorithm. Returns one optimizer, not the
#' barycenter of the median set.
#'
#' @param X Nonempty finite numeric n by d matrix or data frame.
#' @param method Either `"approx"` (default) or `"exact"`. Exact computation
#'   currently requires d = 2 and at least two observations.
#' @param center_method Either `"tukey"` or `"mean"`; ignored when `center` is
#'   supplied. Tukey centering uses [mrfDepth::hdepthmedian()] with `maxdir = NULL`
#'   (the sample median in dimension one). Above dimension two this center is
#'   approximate. If centering fails, supply `center` or use `"mean"`.
#' @param center Optional finite numeric center vector of length d.
#' @param U Approximate method: optional d by L matrix of nonzero directions.
#'   Columns are normalized; no directions are appended. Repeated and antipodal
#'   columns are permitted.
#' @param n_directions Approximate method: total automatic directions, at least
#'   d, including axes. Ignored when `U` is supplied.
#' @param seed Approximate method: optional seed for automatic directions. The
#'   RNG state at direction generation is restored afterwards. This does not
#'   seed Tukey centering. Without a seed, directions consume the ordinary RNG.
#' @param domain Either `"pd"` or `"psd"`. NULL selects `"pd"` for approximation
#'   and `"psd"` for exact computation. Exact currently supports only `"psd"`.
#'   Approximate PD requires a validated strictly positive eigenvalue margin.
#' @param solver Exact method: solver name recognized by [CVXR::psolve()].
#' @param initial_directions Exact method: positive number of initial directions.
#' @param cuts_per_iteration Exact method: maximum cuts added per iteration.
#' @param max_iterations Exact method: maximum iterations per depth level.
#' @param separation_tol Exact method: nonnegative relative oracle tolerance.
#' @param control Approximate method: named list with `abs_tol`, `rel_tol`,
#'   `pd_tol` (defaults 1e-8), `reuse` (TRUE), `kernel` (`"cpp"` or `"R"`),
#'   and `solver` (a list of Clarabel settings). Tolerances apply in internally
#'   scaled scatter units. Explicit solver verbosity overrides `verbose`.
#' @param verbose Whether to report progress.
#'
#' @details The exact backend constructs the complete projective angular
#'   arrangement and searches depth levels with a cutting-plane algorithm.
#'   Certification uses numerical tolerances, not exact arithmetic.
#'   Approximation optimizes inclusive depth on a fixed finite direction set.
#'   Its attained depth is an upper bound on full directional depth of the
#'   returned matrix. A complete finite search does not certify the full optimum.
#'   Unresolved numerical outcomes preserve supported search bounds.
#'   Method-specific settings are used only by the corresponding backend.
#'
#' @return An object of class `"SHDmedian"` with shared fields `method`
#'   (`"approx"` or `"exact"`), `algorithm`, `status`, `certified`, `n`, `d`,
#'   `center`, `center_method`, `Sigma`, `best_k`, `depth`, `k_min`, `k_max`,
#'   `search_bounds`, and `timing` (including `center` and `total` in seconds).
#'   For approximation, `best_k` is the selected count and `depth` is the raw
#'   attained finite-direction depth (`depth_U`); these can differ under boundary
#'   ambiguity. `certified` is always FALSE for approximation. Exact results may
#'   have unavailable `Sigma` or `best_k`. Exact diagnostics are in `search_details`.
#'   Approximation retains `U`, `domain`, `selected_k`, `attained_k`, `depth_U`,
#'   `levels`, raw and tolerance-expanded counts, eigenvalue and constraint
#'   diagnostics, `tolerances`, `scaling`, and `solver` metadata. For exact results,
#'   `search_bounds` contains the initial search interval; for approximation it
#'   contains the final supported optimum bracket.
#'
#' @examples
#' if (requireNamespace("clarabel", quietly = TRUE)) {
#'   X <- matrix(c(-3, -2, -1, 1, 2, 3))
#'   fit <- SHDmedian(X, center = 0, U = matrix(1))
#'   print(fit)
#' }
#' \dontrun{
#' X <- matrix(rnorm(60), ncol = 2)
#' fit <- SHDmedian(X, method = "exact", center_method = "mean")
#' X3 <- matrix(rnorm(90), ncol = 3)
#' fit3 <- SHDmedian(X3, center_method = "mean", n_directions = 100, seed = 42)
#' }
#' @export
SHDmedian <- function(X, method = c("approx", "exact"),
                      center_method = c("tukey", "mean"), center = NULL,
                      U = NULL, n_directions = 2000L, seed = NULL, domain = NULL,
                      solver = "CLARABEL", initial_directions = 20L,
                      cuts_per_iteration = 10L, max_iterations = 200L,
                      separation_tol = 1e-8, control = list(), verbose = FALSE) {
  start <- proc.time()[["elapsed"]]
  method <- match.arg(method)
  X <- as_sample_matrix(X)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  # Add future exact dimensions in the dispatcher, independently of this API.
  backend <- if (method == "approx") shd_approx else shd_exact_backend(ncol(X))
  if (method == "exact" && nrow(X) < 2L)
    stop("X must contain at least two observations.", call. = FALSE)
  if (is.null(domain)) domain <- if (method == "approx") "pd" else "psd"
  domain <- match.arg(domain, c("pd", "psd"))
  if (method == "exact" && domain != "psd")
    stop("The exact method currently supports only domain = 'psd'.", call. = FALSE)
  centered <- shd_center(X, center_method, center)
  center_seconds <- proc.time()[["elapsed"]] - start
  if (method == "approx") {
    if (is.list(control) && (is.null(control$solver) || is.list(control$solver)) &&
        is.null(control$solver$verbose)) control$solver$verbose <- verbose
    result <- backend(X, centered$center, U, n_directions, seed, domain, control)
    result$best_k <- result$selected_k
    result$depth <- result$depth_U
    result$k_min <- unname(result$search_bounds[1L])
    result$k_max <- unname(result$search_bounds[2L])
    result$certified <- FALSE
  } else {
    result <- backend(X, center = centered$center, solver = solver,
      initial_directions = initial_directions, cuts_per_iteration = cuts_per_iteration,
      max_iterations = max_iterations, separation_tol = separation_tol, verbose = verbose)
    result$search_bounds <- c(lower = result$k_min, upper = result$k_max)
  }
  result$algorithm <- result$method
  result$method <- method
  result$d <- ncol(X)
  result$center_method <- centered$method
  result$timing$preprocessing <- result$timing$preprocessing + center_seconds
  result$timing$center <- center_seconds
  result$timing$total <- proc.time()[["elapsed"]] - start
  class(result) <- c("SHDmedian", "list")
  result
}

shd_exact_backend <- function(d) {
  switch(as.character(d), `2` = shd_exact_d2,
    stop("method = 'exact' currently requires exactly two columns.", call. = FALSE))
}
