shd_exact_d2 <- function(
    X,
    center,
    solver = "CLARABEL",
    initial_directions = 20L,
    cuts_per_iteration = 10L,
    max_iterations = 200L,
    separation_tol = 1e-8,
    verbose = FALSE, control = list()) {
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

  if (toupper(solver) != "CLARABEL")
    stop("The exact SOCP method requires solver = 'CLARABEL'.", call. = FALSE)
  allowed <- c("cache_geometry", "cache_limit_mb", "angle_tol", "solver_tol",
               "validation_tol", "pd_tol", "solver")
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || any(!names(control) %in% allowed))))
    stop("Unknown exact control setting.", call. = FALSE)
  args <- control
  settings <- args[["solver", exact = TRUE]]
  if (!is.null(settings) && !is.list(settings))
    stop("control$solver must be a list of Clarabel settings.", call. = FALSE)
  args[["solver"]] <- NULL
  args[["control"]] <- if (is.null(settings)) list() else settings
  fit <- do.call(SHDmedian_socp, c(list(X = X, mu = center, verbose = verbose), args))
  list(method = "point-based-socp", status = fit$status, certified = fit$certified,
       n = fit$n, center = fit$mu, Sigma = fit$Sigma, best_k = fit$k,
       depth = fit$depth, k_min = unname(fit$k_interval[1L]),
       k_max = unname(fit$k_interval[2L]),
       timing = list(preprocessing = unname(fit$timing["geometry"]),
                     geometry = unname(fit$timing["geometry"]),
                     search = unname(sum(fit$timing[c("build", "solver", "validation")])),
                     total = proc.time()[["elapsed"]] - total_start),
       search_details = fit, domain = "pd")
}
