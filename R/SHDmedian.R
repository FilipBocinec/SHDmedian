# Main public interface ---------------------------------------------------

#' Compute the bivariate scatter halfspace median
#'
#' `SHDmedian()` computes a scatter halfspace median for bivariate data using
#' a cutting-plane algorithm. The data are first centered, the complete
#' projective angular arrangement is constructed, and the deepest nonempty
#' depth region is found by binary search.
#'
#' The current implementation supports only numeric data with two columns.
#'
#' @param X An n x 2 numeric matrix or data frame.
#' @param center_method Method used to center the observations. Either
#'   `"tukey"` or `"mean"`.
#' @param mu An optional fixed center vector of length two. If supplied,
#'   `center_method` is ignored.
#' @param solver A solver name recognized by [CVXR::psolve()]. The default is
#'   `"CLARABEL"`.
#' @param initial_directions A positive integer giving the requested number
#'   of initial projective directions.
#' @param cuts_per_iteration A positive integer giving the maximum number of
#'   violated cuts added in one cutting-plane iteration.
#' @param max_iterations A positive integer giving the maximum number of
#'   cutting-plane iterations at one depth level.
#' @param separation_tol A nonnegative relative tolerance used by the
#'   separation oracle.
#' @param verbose Whether to report preprocessing, level-search, and
#'   cutting-plane progress.
#'
#' @return An object of class `"SHDmedian"` containing:
#'
#' - `method`: the computational method;
#' - `status`: the final search status;
#' - `certified`: whether the result was certified by the complete angular
#'   oracle;
#' - `n`: the sample size;
#' - `center`: the center used to construct the centered observations;
#' - `Sigma`: the scatter halfspace median;
#' - `best_k`: the optimal integer depth count;
#' - `depth`: `best_k / n`;
#' - `k_min`, `k_max`: the searched depth bounds;
#' - `timing`: preprocessing, optimization, and total elapsed times;
#' - `search_details`: detailed level-search diagnostics.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rnorm(40), ncol = 2)
#'
#' result <- SHDmedian(
#'   X,
#'   center_method = "mean",
#'   verbose = TRUE
#' )
#'
#' print(result)
#' }
#'
#' @export
SHDmedian <- function(
    X,
    center_method = c("tukey", "mean"),
    mu = NULL,
    solver = "CLARABEL",
    initial_directions = 20L,
    cuts_per_iteration = 10L,
    max_iterations = 200L,
    separation_tol = 1e-8,
    verbose = FALSE) {
  total_start <- proc.time()[["elapsed"]]

  X <- as_bivariate_matrix_d2(X, "X")
  center_method <- match.arg(center_method)

  if (
    !is.character(solver) ||
    length(solver) != 1L ||
    is.na(solver) ||
    !nzchar(solver)
  ) {
    stop(
      "solver must be a nonempty character string.",
      call. = FALSE
    )
  }

  positive_integer_arguments <- list(
    initial_directions = initial_directions,
    cuts_per_iteration = cuts_per_iteration,
    max_iterations = max_iterations
  )

  for (argument_name in names(positive_integer_arguments)) {
    argument_value <-
      positive_integer_arguments[[argument_name]]

    valid <-
      is.numeric(argument_value) &&
      length(argument_value) == 1L &&
      !is.na(argument_value) &&
      is.finite(argument_value) &&
      argument_value == round(argument_value) &&
      argument_value >= 1L

    if (!valid) {
      stop(
        sprintf(
          "%s must be a positive integer.",
          argument_name
        ),
        call. = FALSE
      )
    }
  }

  initial_directions <- as.integer(
    initial_directions
  )

  cuts_per_iteration <- as.integer(
    cuts_per_iteration
  )

  max_iterations <- as.integer(
    max_iterations
  )

  if (
    !is.numeric(separation_tol) ||
    length(separation_tol) != 1L ||
    is.na(separation_tol) ||
    !is.finite(separation_tol) ||
    separation_tol < 0
  ) {
    stop(
      "separation_tol must be a finite nonnegative number.",
      call. = FALSE
    )
  }

  if (
    !is.logical(verbose) ||
    length(verbose) != 1L ||
    is.na(verbose)
  ) {
    stop(
      "verbose must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  n <- nrow(X)

  # Keep internal geometric tolerances fixed in the first package version.
  angle_tol <- 1e-12
  duplicate_tol <- 1e-10

  center_start <- proc.time()[["elapsed"]]

  centered <- center_bivariate_data_d2(
    X = X,
    method = center_method,
    mu = mu
  )

  center_seconds <-
    proc.time()[["elapsed"]] -
    center_start

  if (verbose) {
    message(sprintf(
      "Centered %d observations using the %s center.",
      n,
      centered$method
    ))
  }

  arrangement_start <- proc.time()[["elapsed"]]

  arrangement <- angular_arrangement_d2(
    Y = centered$Y,
    tol = angle_tol
  )

  arrangement_seconds <-
    proc.time()[["elapsed"]] -
    arrangement_start

  if (verbose) {
    message(sprintf(
      paste0(
        "Constructed %d angular cells from ",
        "%d unique interior critical angles."
      ),
      arrangement$number_of_cells,
      arrangement$number_of_unique_interior_angles
    ))
  }

  # Conservative certified bounds for the current bivariate implementation.
  k_min <- max(
    1L,
    as.integer(ceiling(n / 4))
  )

  k_max <- as.integer(
    ceiling(n / 2)
  )

  search <- search_levels_d2_cutting_plane(
    arrangement = arrangement,
    Y = centered$Y,
    k_min = k_min,
    k_max = k_max,
    solver = solver,
    initial_directions = initial_directions,
    cuts_per_iteration = cuts_per_iteration,
    max_iterations = max_iterations,
    separation_tol = separation_tol,
    angle_tol = angle_tol,
    duplicate_tol = duplicate_tol,
    solver_verbose = FALSE,
    verbose = verbose
  )

  total_seconds <-
    proc.time()[["elapsed"]] -
    total_start

  timing <- list(
    center = center_seconds,
    arrangement = arrangement_seconds,
    preprocessing =
      center_seconds + arrangement_seconds,
    search = search$elapsed_seconds,
    total = total_seconds
  )

  result <- list(
    method = "cutting-plane",
    status = search$status,
    certified = isTRUE(search$certified),
    n = n,
    center = centered$mu,
    Sigma = search$Sigma,
    best_k = search$best_k,
    depth = search$depth,
    k_min = k_min,
    k_max = k_max,
    timing = timing,
    search_details = search
  )

  class(result) <- c(
    "SHDmedian",
    "list"
  )

  result
}
