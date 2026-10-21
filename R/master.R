# Cutting-plane master problem --------------------------------------------

#' Compute a safe redundant trace bound
#'
#' @noRd
compute_trace_bound_d2 <- function(level_cells) {
  if (!is.list(level_cells) || length(level_cells) == 0L) {
    stop(
      "level_cells must contain at least one cell.",
      call. = FALSE
    )
  }

  upper_norms_squared <- vapply(
    level_cells,
    function(cell) {
      sum(diag(cell$Zhj))
    },
    numeric(1L)
  )

  if (
    anyNA(upper_norms_squared) ||
    any(!is.finite(upper_norms_squared)) ||
    any(upper_norms_squared < 0)
  ) {
    stop(
      "The upper constraint matrices are invalid.",
      call. = FALSE
    )
  }

  maximum_upper_norm_squared <- max(
    upper_norms_squared
  )

  trace_bound <-
    2 * maximum_upper_norm_squared

  trace_bound *
    (1 + 100 * .Machine$double.eps) +
    100 * .Machine$double.eps
}


#' Solve the current linear-cut master problem
#'
#' @noRd
solve_direction_master_d2 <- function(
    cuts,
    trace_bound,
    solver = "CLARABEL",
    solver_verbose = FALSE) {
  if (!is.list(cuts) || length(cuts) == 0L) {
    stop(
      "cuts must contain at least one direction cut.",
      call. = FALSE
    )
  }

  if (
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

  s11 <- CVXR::Variable(name = "s11")
  s12 <- CVXR::Variable(name = "s12")
  s22 <- CVXR::Variable(name = "s22")

  common_margin <- CVXR::Variable(
    name = "common_margin",
    nonneg = TRUE
  )

  constraints <- vector(
    mode = "list",
    length = 2L + length(cuts)
  )

  # A symmetric 2 x 2 matrix is PSD exactly when
  #
  # ||(2 s12, s11 - s22)||_2 <= s11 + s22.
  constraints[[1L]] <- CVXR::SOC(
    s11 + s22,
    CVXR::vstack(
      2 * s12,
      s11 - s22
    )
  )

  constraints[[2L]] <-
    s11 + s22 <= trace_bound

  for (cut_index in seq_along(cuts)) {
    cut <- cuts[[cut_index]]

    projected_scatter <-
      cut$coefficients["s11"] * s11 +
      cut$coefficients["s12"] * s12 +
      cut$coefficients["s22"] * s22

    if (cut$side == "lower") {
      constraints[[2L + cut_index]] <-
        projected_scatter -
        cut$right_hand_side >=
        common_margin
    } else if (cut$side == "upper") {
      constraints[[2L + cut_index]] <-
        cut$right_hand_side -
        projected_scatter >=
        common_margin
    } else {
      stop(
        "Every cut side must be 'lower' or 'upper'.",
        call. = FALSE
      )
    }
  }

  problem <- CVXR::Problem(
    objective = CVXR::Maximize(common_margin),
    constraints = constraints
  )

  solver_error <- NULL

  objective_value <- tryCatch(
    CVXR::psolve(
      problem,
      solver = solver,
      verbose = solver_verbose
    ),
    error = function(error) {
      solver_error <<- conditionMessage(error)
      NA_real_
    }
  )

  if (!is.null(solver_error)) {
    return(list(
      status = "solver_error",
      solver_status = NA_character_,
      solver_error = solver_error,
      objective_value = objective_value,
      margin = NA_real_,
      Sigma = NULL,
      problem = problem
    ))
  }

  solver_status <- CVXR::status(problem)

  optimal <- solver_status %in% c(
    "optimal",
    "optimal_inaccurate"
  )

  Sigma_hat <- NULL
  margin_hat <- NA_real_

  if (optimal) {
    Sigma_hat <- matrix(
      c(
        as.numeric(CVXR::value(s11)),
        as.numeric(CVXR::value(s12)),
        as.numeric(CVXR::value(s12)),
        as.numeric(CVXR::value(s22))
      ),
      nrow = 2L,
      byrow = TRUE
    )

    Sigma_hat <- (
      Sigma_hat + t(Sigma_hat)
    ) / 2

    margin_hat <- as.numeric(
      CVXR::value(common_margin)
    )
  }

  list(
    status = if (optimal) {
      "optimal"
    } else {
      as.character(solver_status)
    },
    solver_status =
      as.character(solver_status),
    solver_error = NULL,
    objective_value = objective_value,
    margin = margin_hat,
    Sigma = Sigma_hat,
    problem = problem
  )
}
