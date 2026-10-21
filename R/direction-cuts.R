# Direction cuts ----------------------------------------------------------

#' Construct one linear direction cut
#'
#' @noRd
make_direction_cut_d2 <- function(
    cell,
    cell_index,
    theta,
    side = c("lower", "upper"),
    iteration = 0L,
    source = "oracle") {
  side <- match.arg(side)

  if (
    !is.numeric(theta) ||
    length(theta) != 1L ||
    is.na(theta) ||
    !is.finite(theta)
  ) {
    stop(
      "theta must be a finite numeric scalar.",
      call. = FALSE
    )
  }

  direction <- c(
    cos(theta),
    sin(theta)
  )

  coefficients <- c(
    s11 = direction[1L]^2,
    s12 = 2 * direction[1L] * direction[2L],
    s22 = direction[2L]^2
  )

  active_matrix <- if (side == "lower") {
    cell$Zlj
  } else {
    cell$Zhj
  }

  right_hand_side <- drop(
    crossprod(
      direction,
      active_matrix %*% direction
    )
  )

  list(
    cell_index = as.integer(cell_index),
    side = side,
    theta = as.numeric(theta),
    theta_projective = as.numeric(theta %% pi),
    coefficients = coefficients,
    right_hand_side = as.numeric(right_hand_side),
    iteration = as.integer(iteration),
    source = as.character(source)
  )
}


#' Determine whether a direction cut is already present
#'
#' @noRd
direction_cut_is_duplicate_d2 <- function(
    candidate,
    cuts,
    tol = 1e-10) {
  if (length(cuts) == 0L) {
    return(FALSE)
  }

  for (existing_cut in cuts) {
    if (candidate$side != existing_cut$side) {
      next
    }

    coefficient_scale <- max(
      1,
      abs(candidate$coefficients),
      abs(existing_cut$coefficients)
    )

    right_hand_side_scale <- max(
      1,
      abs(candidate$right_hand_side),
      abs(existing_cut$right_hand_side)
    )

    same_coefficients <-
      max(abs(
        candidate$coefficients -
          existing_cut$coefficients
      )) <= tol * coefficient_scale

    same_right_hand_side <-
      abs(
        candidate$right_hand_side -
          existing_cut$right_hand_side
      ) <=
      tol * right_hand_side_scale

    if (same_coefficients && same_right_hand_side) {
      return(TRUE)
    }
  }

  FALSE
}


#' Convert direction cuts to a data frame
#'
#' @noRd
direction_cuts_to_data_frame_d2 <- function(cuts) {
  if (length(cuts) == 0L) {
    return(data.frame())
  }

  rows <- lapply(
    seq_along(cuts),
    function(cut_index) {
      cut <- cuts[[cut_index]]

      data.frame(
        cut_index = cut_index,
        cell_index = cut$cell_index,
        side = cut$side,
        theta = cut$theta,
        theta_projective =
          cut$theta_projective,
        coefficient_s11 =
          unname(cut$coefficients["s11"]),
        coefficient_s12 =
          unname(cut$coefficients["s12"]),
        coefficient_s22 =
          unname(cut$coefficients["s22"]),
        right_hand_side =
          cut$right_hand_side,
        iteration = cut$iteration,
        source = cut$source,
        stringsAsFactors = FALSE
      )
    }
  )

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}


#' Select well-spread initial angular cells
#'
#' @noRd
select_initial_cells_d2 <- function(
    level_cells,
    number_of_directions = 8L) {
  if (!is.list(level_cells) || length(level_cells) == 0L) {
    stop(
      "level_cells must contain at least one cell.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(number_of_directions) ||
    length(number_of_directions) != 1L ||
    is.na(number_of_directions) ||
    !is.finite(number_of_directions) ||
    number_of_directions != round(number_of_directions) ||
    number_of_directions < 1L
  ) {
    stop(
      "number_of_directions must be a positive integer.",
      call. = FALSE
    )
  }

  number_of_directions <- as.integer(
    number_of_directions
  )

  cell_midpoints <- vapply(
    level_cells,
    function(cell) {
      (
        (cell$alpha_j + cell$alpha_j1) / 2
      ) %% pi
    },
    numeric(1L)
  )

  target_angles <-
    (seq_len(number_of_directions) - 1L) *
    pi / number_of_directions

  selected_indices <- vapply(
    target_angles,
    function(target_angle) {
      raw_distance <- abs(
        cell_midpoints - target_angle
      )

      projective_distance <- pmin(
        raw_distance,
        pi - raw_distance
      )

      which.min(projective_distance)
    },
    integer(1L)
  )

  unique(selected_indices)
}


#' Initialize the cutting-plane master
#'
#' @noRd
initialize_direction_cuts_d2 <- function(
    level_cells,
    number_of_directions = 8L,
    duplicate_tol = 1e-10) {
  selected_cells <- select_initial_cells_d2(
    level_cells = level_cells,
    number_of_directions = number_of_directions
  )

  cuts <- list()

  for (cell_index in selected_cells) {
    cell <- level_cells[[cell_index]]

    theta <-
      (cell$alpha_j + cell$alpha_j1) / 2

    for (side in c("lower", "upper")) {
      candidate <- make_direction_cut_d2(
        cell = cell,
        cell_index = cell_index,
        theta = theta,
        side = side,
        iteration = 0L,
        source = "initial"
      )

      is_duplicate <- direction_cut_is_duplicate_d2(
        candidate = candidate,
        cuts = cuts,
        tol = duplicate_tol
      )

      if (!is_duplicate) {
        cuts[[length(cuts) + 1L]] <- candidate
      }
    }
  }

  cuts
}
