# Analytic separation oracle ----------------------------------------------

#' Evaluate a symmetric quadratic form at angular directions
#'
#' @param Q A numeric 2 x 2 matrix.
#' @param theta A numeric vector of angles.
#'
#' @return Values of u(theta)^T Q u(theta).
#'
#' @noRd
quadratic_form_at_angle_d2 <- function(Q, theta) {
  cosine <- cos(theta)
  sine <- sin(theta)

  Q[1L, 1L] * cosine^2 +
    2 * Q[1L, 2L] * cosine * sine +
    Q[2L, 2L] * sine^2
}


#' Minimize a quadratic form over a closed angular interval
#'
#' @param Q A symmetric numeric 2 x 2 matrix.
#' @param alpha_left The left endpoint of the interval.
#' @param alpha_right The right endpoint of the interval.
#' @param tol A nonnegative numerical tolerance.
#'
#' @return A list containing the minimum, a minimizing angle, and all
#'   evaluated candidate angles.
#'
#' @noRd
minimize_quadratic_on_arc_d2 <- function(
    Q,
    alpha_left,
    alpha_right,
    tol = 1e-12) {
  Q <- as.matrix(Q)

  if (
    !is.numeric(Q) ||
    !all(dim(Q) == c(2L, 2L)) ||
    anyNA(Q) ||
    any(!is.finite(Q))
  ) {
    stop(
      "Q must be a finite numeric 2 x 2 matrix.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(alpha_left) ||
    length(alpha_left) != 1L ||
    is.na(alpha_left) ||
    !is.finite(alpha_left)
  ) {
    stop(
      "alpha_left must be a finite numeric scalar.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(alpha_right) ||
    length(alpha_right) != 1L ||
    is.na(alpha_right) ||
    !is.finite(alpha_right)
  ) {
    stop(
      "alpha_right must be a finite numeric scalar.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(tol) ||
    length(tol) != 1L ||
    is.na(tol) ||
    !is.finite(tol) ||
    tol < 0
  ) {
    stop(
      "tol must be a finite nonnegative number.",
      call. = FALSE
    )
  }

  if (alpha_right < alpha_left - tol) {
    stop(
      "alpha_right must not be smaller than alpha_left.",
      call. = FALSE
    )
  }

  interval_width <- alpha_right - alpha_left

  if (interval_width > pi + tol) {
    stop(
      "The angular interval must have width at most pi.",
      call. = FALSE
    )
  }

  # Remove negligible numerical asymmetry.
  Q <- (Q + t(Q)) / 2

  # Write
  #
  # q(theta) = A + B cos(2 theta) + C sin(2 theta).
  coefficient_A <-
    (Q[1L, 1L] + Q[2L, 2L]) / 2

  coefficient_B <-
    (Q[1L, 1L] - Q[2L, 2L]) / 2

  coefficient_C <- Q[1L, 2L]

  amplitude <- sqrt(
    coefficient_B^2 +
      coefficient_C^2
  )

  candidates <- c(
    alpha_left,
    alpha_right
  )

  coefficient_scale <- max(
    1,
    abs(coefficient_A),
    amplitude
  )

  if (amplitude > tol * coefficient_scale) {
    # The derivative is zero at
    #
    # base_angle + m * pi / 2, m integer.
    base_angle <-
      atan2(
        coefficient_C,
        coefficient_B
      ) / 2

    stationary_step <- pi / 2

    first_index <- ceiling(
      (alpha_left - base_angle) /
        stationary_step -
        10 * tol
    )

    last_index <- floor(
      (alpha_right - base_angle) /
        stationary_step +
        10 * tol
    )

    if (first_index <= last_index) {
      stationary_angles <-
        base_angle +
        seq.int(
          first_index,
          last_index
        ) * stationary_step

      stationary_angles <- stationary_angles[
        stationary_angles >= alpha_left - tol &
          stationary_angles <= alpha_right + tol
      ]

      candidates <- c(
        candidates,
        stationary_angles
      )
    }
  }

  candidate_values <- quadratic_form_at_angle_d2(
    Q,
    candidates
  )

  minimum_index <- which.min(candidate_values)

  list(
    value = as.numeric(
      candidate_values[minimum_index]
    ),
    theta = as.numeric(
      candidates[minimum_index]
    ),
    candidate_angles = candidates,
    candidate_values = candidate_values
  )
}


#' Check a candidate scatter matrix over all level-specific cells
#'
#' @param Sigma A numeric 2 x 2 candidate scatter matrix.
#' @param level_cells A list returned by prepare_level_cells_d2().
#' @param separation_tol A relative feasibility tolerance.
#' @param angle_tol A tolerance passed to the analytic arc minimizer.
#'
#' @return A list containing the feasibility decision, violations, margins,
#'   and positive-semidefiniteness diagnostics.
#'
#' @noRd
separate_level_cells_d2 <- function(
    Sigma,
    level_cells,
    separation_tol = 1e-8,
    angle_tol = 1e-12) {
  Sigma <- as.matrix(Sigma)

  if (
    !is.numeric(Sigma) ||
    !all(dim(Sigma) == c(2L, 2L)) ||
    anyNA(Sigma) ||
    any(!is.finite(Sigma))
  ) {
    stop(
      "Sigma must be a finite numeric 2 x 2 matrix.",
      call. = FALSE
    )
  }

  if (!is.list(level_cells) || length(level_cells) == 0L) {
    stop(
      "level_cells must contain at least one cell.",
      call. = FALSE
    )
  }

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
    !is.numeric(angle_tol) ||
    length(angle_tol) != 1L ||
    is.na(angle_tol) ||
    !is.finite(angle_tol) ||
    angle_tol < 0
  ) {
    stop(
      "angle_tol must be a finite nonnegative number.",
      call. = FALSE
    )
  }

  Sigma <- (Sigma + t(Sigma)) / 2

  checks <- vector(
    "list",
    2L * length(level_cells)
  )

  check_index <- 0L

  for (cell_index in seq_along(level_cells)) {
    cell <- level_cells[[cell_index]]

    alpha_left <- cell$alpha_j
    alpha_right <- cell$alpha_j1

    Q_lower <- Sigma - cell$Zlj
    Q_upper <- cell$Zhj - Sigma

    lower_minimum <- minimize_quadratic_on_arc_d2(
      Q = Q_lower,
      alpha_left = alpha_left,
      alpha_right = alpha_right,
      tol = angle_tol
    )

    upper_minimum <- minimize_quadratic_on_arc_d2(
      Q = Q_upper,
      alpha_left = alpha_left,
      alpha_right = alpha_right,
      tol = angle_tol
    )

    lower_scale <- max(
      1,
      abs(Sigma),
      abs(cell$Zlj)
    )

    upper_scale <- max(
      1,
      abs(Sigma),
      abs(cell$Zhj)
    )

    check_index <- check_index + 1L

    checks[[check_index]] <- data.frame(
      cell_index = cell_index,
      side = "lower",
      theta = lower_minimum$theta,
      theta_projective =
        lower_minimum$theta %% pi,
      margin = lower_minimum$value,
      scale = lower_scale,
      relative_margin =
        lower_minimum$value / lower_scale,
      tolerance =
        separation_tol * lower_scale,
      violated =
        lower_minimum$value <
        -separation_tol * lower_scale,
      stringsAsFactors = FALSE
    )

    check_index <- check_index + 1L

    checks[[check_index]] <- data.frame(
      cell_index = cell_index,
      side = "upper",
      theta = upper_minimum$theta,
      theta_projective =
        upper_minimum$theta %% pi,
      margin = upper_minimum$value,
      scale = upper_scale,
      relative_margin =
        upper_minimum$value / upper_scale,
      tolerance =
        separation_tol * upper_scale,
      violated =
        upper_minimum$value <
        -separation_tol * upper_scale,
      stringsAsFactors = FALSE
    )
  }

  checks <- do.call(
    rbind,
    checks
  )

  checks <- checks[
    order(
      checks$relative_margin,
      checks$margin
    ),
    ,
    drop = FALSE
  ]

  rownames(checks) <- NULL

  Sigma_eigenvalues <- eigen(
    Sigma,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  psd_margin <- min(Sigma_eigenvalues)
  psd_scale <- max(1, abs(Sigma))
  psd_relative_margin <- psd_margin / psd_scale

  psd_feasible <-
    psd_margin >=
    -separation_tol * psd_scale

  number_of_cell_violations <- sum(
    checks$violated
  )

  list(
    feasible =
      !any(checks$violated) &&
      psd_feasible,
    number_violated =
      number_of_cell_violations +
      as.integer(!psd_feasible),
    number_of_cell_violations =
      number_of_cell_violations,
    worst_margin = min(checks$margin),
    worst_relative_margin = min(
      checks$relative_margin,
      psd_relative_margin
    ),
    worst_check =
      checks[1L, , drop = FALSE],
    psd_feasible = psd_feasible,
    psd_margin = psd_margin,
    psd_relative_margin =
      psd_relative_margin,
    checks = checks
  )
}
