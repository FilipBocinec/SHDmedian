# Binary search over depth levels -----------------------------------------

#' Find the deepest feasible level by binary search
#'
#' @noRd
search_levels_d2_cutting_plane <- function(
    arrangement,
    Y,
    k_min,
    k_max,
    solver = "CLARABEL",
    initial_directions = 8L,
    cuts_per_iteration = 10L,
    max_iterations = 100L,
    separation_tol = 1e-8,
    angle_tol = 1e-12,
    duplicate_tol = 1e-10,
    solver_verbose = FALSE,
    verbose = TRUE) {
  Y <- as_bivariate_matrix_d2(Y, "Y")
  n <- nrow(Y)

  if (!inherits(
    arrangement,
    "shd_angular_arrangement_d2"
  )) {
    stop(
      "arrangement must be returned by angular_arrangement_d2().",
      call. = FALSE
    )
  }

  validate_search_bound <- function(value, name) {
    valid <-
      is.numeric(value) &&
      length(value) == 1L &&
      !is.na(value) &&
      is.finite(value) &&
      value == round(value)

    if (!valid) {
      stop(
        sprintf(
          "%s must be a finite integer scalar.",
          name
        ),
        call. = FALSE
      )
    }

    as.integer(value)
  }

  k_min <- validate_search_bound(
    k_min,
    "k_min"
  )

  k_max <- validate_search_bound(
    k_max,
    "k_max"
  )

  if (
    k_min < 1L ||
    k_max > n ||
    k_min > k_max
  ) {
    stop(
      paste0(
        "The search bounds must satisfy ",
        "1 <= k_min <= k_max <= n."
      ),
      call. = FALSE
    )
  }

  if (
    arrangement$n != n ||
    ncol(arrangement$orderings) != n
  ) {
    stop(
      "The arrangement and Y contain different numbers of observations.",
      call. = FALSE
    )
  }

  search_start <- proc.time()[["elapsed"]]

  lower <- k_min
  upper <- k_max

  best_k <- NA_integer_
  best_result <- NULL

  tested_levels <- integer(0L)
  level_results <- list()

  finish_search <- function(
    status,
    unresolved_k = NA_integer_) {
    if (length(level_results) > 0L) {
      names(level_results) <-
        as.character(tested_levels)
    }

    search_certified <- status %in% c(
      "optimal",
      "no_feasible_level"
    )

    list(
      method = "cutting-plane",
      status = status,
      certified = search_certified,
      search_bounds = c(
        k_min = k_min,
        k_max = k_max
      ),
      best_k = best_k,
      depth = if (!is.na(best_k)) {
        best_k / n
      } else {
        NA_real_
      },
      Sigma = if (!is.null(best_result)) {
        best_result$Sigma
      } else {
        NULL
      },
      best_result = best_result,
      unresolved_k = unresolved_k,
      tested_levels = tested_levels,
      number_of_tested_levels =
        length(tested_levels),
      level_results = level_results,
      elapsed_seconds =
        proc.time()[["elapsed"]] -
        search_start
    )
  }

  while (lower <= upper) {
    k <- lower + (upper - lower) %/% 2L

    if (verbose) {
      message(sprintf(
        paste0(
          "Testing level k = %d within ",
          "the current interval [%d, %d]."
        ),
        k,
        lower,
        upper
      ))
    }

    level_cells <- prepare_level_cells_d2(
      arrangement = arrangement,
      Y = Y,
      k = k,
      tol = angle_tol,
      merge_wrap = TRUE,
      verbose = verbose
    )

    level_result <- solve_level_d2_cutting_plane(
      level_cells = level_cells,
      solver = solver,
      initial_directions =
        initial_directions,
      cuts_per_iteration =
        cuts_per_iteration,
      max_iterations = max_iterations,
      separation_tol = separation_tol,
      angle_tol = angle_tol,
      duplicate_tol = duplicate_tol,
      solver_verbose = solver_verbose,
      verbose = verbose
    )

    tested_levels <- c(
      tested_levels,
      k
    )


    level_results[[length(level_results) + 1L]] <- list(level_result)

    if (
      !isTRUE(level_result$certified) ||
      is.na(level_result$feasible)
    ) {
      return(finish_search(
        status = "unresolved_level",
        unresolved_k = k
      ))
    }

    if (isTRUE(level_result$feasible)) {
      best_k <- k
      best_result <- level_result
      lower <- k + 1L
    } else {
      upper <- k - 1L
    }
  }

  if (is.na(best_k)) {
    return(finish_search(
      status = "no_feasible_level"
    ))
  }

  finish_search(status = "optimal")
}
