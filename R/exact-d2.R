shd_exact_d2 <- function(
    X,
    center,
    solver = "CLARABEL",
    initial_directions = 20L,
    cuts_per_iteration = 10L,
    max_iterations = 200L,
    separation_tol = 1e-8,
    verbose = FALSE) {
  total_start <- proc.time()[["elapsed"]]

  X <- as_bivariate_matrix_d2(X, "X")

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
      argument_value >= 1L &&
      argument_value <= .Machine$integer.max

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

  centered <- list(mu = center, Y = sweep(X, 2L, center, "-"))
  if (any(!is.finite(centered$Y)))
    stop("Centering overflowed; rescale X and center.", call. = FALSE)

  center_seconds <-
    proc.time()[["elapsed"]] -
    center_start

  if (verbose) {
    message(sprintf(
      "Prepared %d centered observations.",
      n
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
