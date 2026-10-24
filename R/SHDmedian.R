#' Compute a scatter halfspace median
#'
#' Unified interface to finite-direction SDP optimization in any dimension and
#' the exact bivariate point-based SOCP algorithm. Returns one optimizer, not the
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
#'   and for exact computation. Exact currently supports only `"pd"`.
#'   Approximate PD requires a validated strictly positive eigenvalue margin.
#' @param solver Exact method: only `"CLARABEL"` is supported.
#' @param initial_directions,cuts_per_iteration,max_iterations,separation_tol
#'   Legacy cutting-plane settings, validated but unused by SOCP.
#' @param control Approximate method: named list with `abs_tol`, `rel_tol`,
#'   `pd_tol` (defaults 1e-8), `reuse` (TRUE), `kernel` (`"cpp"` or `"R"`),
#'   and `solver` (a list of Clarabel settings). Tolerances apply in internally
#'   scaled scatter units. Explicit solver verbosity overrides `verbose`.
#' @param verbose Whether to report progress.
#' @param gaussian_consistent Logical; if TRUE (default), multiply the raw
#'   scatter median by `1 / stats::qnorm(0.75)^2` for Gaussian consistency.
#'   Set FALSE to return the unscaled scatter median in `Sigma`.
#'
#' @details The exact backend uses per-observation angular sweeps, sparse SOCP
#'   models and binary depth search. Its `control` accepts `cache_geometry`
#'   (TRUE, FALSE or `"auto"`), `cache_limit_mb` (256), `angle_tol` (1e-12),
#'   `solver_tol` (1e-9), `validation_tol` (1e-7), `pd_tol` (1e-9), and
#'   `solver` (a list of Clarabel settings). Singular boundary cases can remain
#'   unresolved. `search_details` includes strict depth, level diagnostics,
#'   numerical bounds, tolerances and component timings.
#'   Certification uses numerical tolerances, not exact arithmetic.
#'   Approximation optimizes inclusive depth on a fixed finite direction set.
#'   Its attained depth is an upper bound on full directional depth of the
#'   raw median matrix. A complete finite search does not certify the full optimum.
#'   Unresolved numerical outcomes preserve supported search bounds.
#'   Method-specific settings are used only by the corresponding backend.
#'   For a Gaussian population with covariance V and the true center, every
#'   squared standardized projection has median `stats::qnorm(0.75)^2`.
#'   Thus the population scatter halfspace median is this constant times V.
#'   Gaussian scaling corrects this population factor; it is not a finite-sample
#'   unbiasedness correction. Depths, certification and all backend diagnostics
#'   refer to `Sigma_raw`, before Gaussian scaling.
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
#'   `search_bounds` contains the final supported optimum bracket; for approximation it
#'   contains the final supported optimum bracket.
#'   `Sigma_raw` stores the unscaled scatter median (NULL if unavailable),
#'   `gaussian_consistent` records the requested correction, and
#'   `consistency_factor` is the multiplier used to obtain `Sigma` from
#'   `Sigma_raw` (one when the correction is disabled).
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
                      separation_tol = 1e-8, control = list(), verbose = FALSE,
                      gaussian_consistent = TRUE) {
  start <- proc.time()[["elapsed"]]
  method <- match.arg(method)
  X <- as_sample_matrix(X)
  if (!is.logical(gaussian_consistent) || length(gaussian_consistent) != 1L ||
      is.na(gaussian_consistent))
    stop("gaussian_consistent must be TRUE or FALSE.", call. = FALSE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  # Add future exact dimensions in the dispatcher, independently of this API.
  backend <- if (method == "approx") shd_approx else shd_exact_backend(ncol(X))
  if (method == "exact" && nrow(X) < 2L)
    stop("X must contain at least two observations.", call. = FALSE)
  if (is.null(domain)) domain <- "pd"
  domain <- match.arg(domain, c("pd", "psd"))
  if (method == "exact" && domain != "pd")
    stop("The exact method currently supports only domain = 'pd'.", call. = FALSE)
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
      max_iterations = max_iterations, separation_tol = separation_tol, verbose = verbose, control = control)
    result$search_bounds <- c(lower = result$k_min, upper = result$k_max)
  }
  result["Sigma_raw"] <- list(result$Sigma)
  result$gaussian_consistent <- gaussian_consistent
  result$consistency_factor <- if (gaussian_consistent) 1 / stats::qnorm(0.75)^2 else 1
  if (!is.null(result$Sigma)) {
    result$Sigma <- result$Sigma * result$consistency_factor
    if (any(!is.finite(result$Sigma)))
      stop("Gaussian-scaled scatter output overflowed; rescale X and center.", call. = FALSE)
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
