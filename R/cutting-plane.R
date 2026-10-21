# Cutting-plane algorithm for one level -----------------------------------

#' Generate cuts from violated oracle checks
#'
#' @noRd
generate_violated_cut_batch_d2 <- function(
    oracle_result,
    level_cells,
    existing_cuts,
    iteration,
    cuts_per_iteration = 10L,
    duplicate_tol = 1e-10) {
  if (
    !is.numeric(cuts_per_iteration) ||
    length(cuts_per_iteration) != 1L ||
    is.na(cuts_per_iteration) ||
    !is.finite(cuts_per_iteration) ||
    cuts_per_iteration != round(cuts_per_iteration) ||
    cuts_per_iteration < 1L
  ) {
    stop(
      "cuts_per_iteration must be a positive integer.",
      call. = FALSE
    )
  }

  cuts_per_iteration <- as.integer(
    cuts_per_iteration
  )

  violated_checks <- oracle_result$checks[
    oracle_result$checks$violated,
    ,
    drop = FALSE
  ]

  if (nrow(violated_checks) == 0L) {
    return(list())
  }

  new_cuts <- list()

  # The checks are already sorted from the worst relative margin.
  for (row_index in seq_len(nrow(violated_checks))) {
    check <- violated_checks[
      row_index,
      ,
      drop = FALSE
    ]

    cell_index <- as.integer(
      check$cell_index[[1L]]
    )

    candidate <- make_direction_cut_d2(
      cell = level_cells[[cell_index]],
      cell_index = cell_index,
      theta = check$theta[[1L]],
      side = check$side[[1L]],
      iteration = iteration,
      source = "oracle"
    )

    all_known_cuts <- c(
      existing_cuts,
      new_cuts
    )

    is_duplicate <- direction_cut_is_duplicate_d2(
      candidate = candidate,
      cuts = all_known_cuts,
      tol = duplicate_tol
    )

    if (!is_duplicate) {
      new_cuts[[length(new_cuts) + 1L]] <-
        candidate
    }

    if (length(new_cuts) >= cuts_per_iteration) {
      break
    }
  }

  new_cuts
}


#' Test one depth level using cutting planes
#'
#' @noRd
solve_level_d2_cutting_plane <- function(
    level_cells,
    solver = "CLARABEL",
    initial_directions = 8L,
    cuts_per_iteration = 10L,
    max_iterations = 100L,
    separation_tol = 1e-8,
    angle_tol = 1e-12,
    duplicate_tol = 1e-10,
    trace_bound = NULL,
    solver_verbose = FALSE,
    verbose = TRUE) {
  if (!is.list(level_cells) || length(level_cells) == 0L) {
    stop(
      "level_cells must contain at least one cell.",
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

  tolerance_arguments <- list(
    separation_tol = separation_tol,
    angle_tol = angle_tol,
    duplicate_tol = duplicate_tol
  )

  for (argument_name in names(tolerance_arguments)) {
    argument_value <-
      tolerance_arguments[[argument_name]]

    valid <-
      is.numeric(argument_value) &&
      length(argument_value) == 1L &&
      !is.na(argument_value) &&
      is.finite(argument_value) &&
      argument_value >= 0

    if (!valid) {
      stop(
        sprintf(
          "%s must be a finite nonnegative number.",
          argument_name
        ),
        call. = FALSE
      )
    }
  }

  if (
    !is.logical(solver_verbose) ||
    length(solver_verbose) != 1L ||
    is.na(solver_verbose)
  ) {
    stop(
      "solver_verbose must be TRUE or FALSE.",
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

  if (is.null(trace_bound)) {
    trace_bound <- compute_trace_bound_d2(
      level_cells
    )
  } else if (
    !is.numeric(trace_bound) ||
    length(trace_bound) != 1L ||
    is.na(trace_bound) ||
    !is.finite(trace_bound) ||
    trace_bound < 0
  ) {
    stop(
      "trace_bound must be a finite nonnegative number.",
      call. = FALSE
    )
  }

  cuts <- initialize_direction_cuts_d2(
    level_cells = level_cells,
    number_of_directions =
      initial_directions,
    duplicate_tol = duplicate_tol
  )

  history <- list()
  final_master <- NULL
  final_oracle <- NULL
  Sigma_hat <- NULL

  time_start <- proc.time()[["elapsed"]]

  finish <- function(
    feasible,
    certified,
    status,
    iteration,
    Sigma = Sigma_hat,
    master = final_master,
    oracle = final_oracle) {
    history_frame <- if (length(history) > 0L) {
      result <- do.call(rbind, history)
      rownames(result) <- NULL
      result
    } else {
      data.frame()
    }

    solver_status <- if (!is.null(master)) {
      master$solver_status
    } else {
      NA_character_
    }

    solver_error <- if (!is.null(master)) {
      master$solver_error
    } else {
      NULL
    }

    list(
      method = "cutting-plane",
      feasible = feasible,
      certified = certified,
      status = status,
      solver_status = solver_status,
      solver_error = solver_error,
      Sigma = Sigma,
      eigenvalues = if (!is.null(Sigma)) {
        eigen(
          Sigma,
          symmetric = TRUE,
          only.values = TRUE
        )$values
      } else {
        NULL
      },
      iterations = as.integer(iteration),
      number_of_cells = length(level_cells),
      number_of_cuts = length(cuts),
      cuts = cuts,
      cut_table =
        direction_cuts_to_data_frame_d2(cuts),
      history = history_frame,
      oracle = oracle,
      master = master,
      trace_bound = trace_bound,
      elapsed_seconds =
        proc.time()[["elapsed"]] - time_start
    )
  }

  for (iteration in seq_len(max_iterations)) {
    master <- solve_direction_master_d2(
      cuts = cuts,
      trace_bound = trace_bound,
      solver = solver,
      solver_verbose = solver_verbose
    )

    final_master <- master

    if (identical(master$status, "solver_error")) {
      return(finish(
        feasible = NA,
        certified = FALSE,
        status = "solver_error",
        iteration = iteration,
        Sigma = NULL,
        master = master,
        oracle = NULL
      ))
    }

    if (identical(master$solver_status, "infeasible")) {
      # If a subset of the complete constraints is already infeasible, the
      # complete level is infeasible as well.
      return(finish(
        feasible = FALSE,
        certified = TRUE,
        status = "infeasible",
        iteration = iteration,
        Sigma = NULL,
        master = master,
        oracle = NULL
      ))
    }

    if (
      !(master$solver_status %in%
        c("optimal", "optimal_inaccurate"))
    ) {
      return(finish(
        feasible = NA,
        certified = FALSE,
        status = "unresolved_master",
        iteration = iteration,
        Sigma = NULL,
        master = master,
        oracle = NULL
      ))
    }

    Sigma_hat <- master$Sigma

    oracle <- separate_level_cells_d2(
      Sigma = Sigma_hat,
      level_cells = level_cells,
      separation_tol = separation_tol,
      angle_tol = angle_tol
    )

    final_oracle <- oracle

    history_row <- data.frame(
      iteration = iteration,
      number_of_cuts = length(cuts),
      master_status = master$solver_status,
      master_margin = master$margin,
      worst_margin = oracle$worst_margin,
      worst_relative_margin =
        oracle$worst_relative_margin,
      number_violated =
        oracle$number_violated,
      cuts_added = 0L,
      stringsAsFactors = FALSE
    )

    if (verbose) {
      message(sprintf(
        paste0(
          "Iteration %d: cuts = %d, ",
          "master margin = %.3e, ",
          "worst relative margin = %.3e, ",
          "violations = %d."
        ),
        iteration,
        length(cuts),
        master$margin,
        oracle$worst_relative_margin,
        oracle$number_violated
      ))
    }

    if (isTRUE(oracle$feasible)) {
      history[[length(history) + 1L]] <-
        history_row

      return(finish(
        feasible = TRUE,
        certified = TRUE,
        status = "feasible",
        iteration = iteration,
        Sigma = Sigma_hat,
        master = master,
        oracle = oracle
      ))
    }

    new_cuts <- generate_violated_cut_batch_d2(
      oracle_result = oracle,
      level_cells = level_cells,
      existing_cuts = cuts,
      iteration = iteration,
      cuts_per_iteration =
        cuts_per_iteration,
      duplicate_tol = duplicate_tol
    )

    history_row$cuts_added <- length(new_cuts)

    history[[length(history) + 1L]] <-
      history_row

    if (length(new_cuts) == 0L) {
      return(finish(
        feasible = NA,
        certified = FALSE,
        status = "numerical_stall",
        iteration = iteration,
        Sigma = Sigma_hat,
        master = master,
        oracle = oracle
      ))
    }

    cuts <- c(cuts, new_cuts)
  }

  finish(
    feasible = NA,
    certified = FALSE,
    status = "maximum_iterations",
    iteration = max_iterations,
    Sigma = Sigma_hat,
    master = final_master,
    oracle = final_oracle
  )
}
