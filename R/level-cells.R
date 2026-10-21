# Level-specific angular cells --------------------------------------------

#' Prepare angular cells for a fixed depth level
#'
#' For a fixed integer level k, the function selects the observations defining
#' the kth and (n - k + 1)th order statistics in every angular cell.
#' Consecutive cells with identical lower and upper matrix constraints are
#' merged.
#'
#' @param arrangement An object returned by angular_arrangement_d2().
#' @param Y An n x 2 matrix of centered observations.
#' @param k An integer depth count.
#' @param tol A nonnegative numerical tolerance.
#' @param merge_wrap Whether matching first and last projective cells should
#'   be merged.
#' @param verbose Whether to report the reduction in the number of cells.
#'
#' @return A list of level-specific angular cells.
#'
#' @noRd
prepare_level_cells_d2 <- function(
    arrangement,
    Y,
    k,
    tol = 1e-12,
    merge_wrap = TRUE,
    verbose = FALSE) {
  Y <- as_bivariate_matrix_d2(Y, "Y")
  n <- nrow(Y)

  if (!inherits(arrangement, "shd_angular_arrangement_d2")) {
    stop(
      "arrangement must be returned by angular_arrangement_d2().",
      call. = FALSE
    )
  }

  if (
    !is.numeric(k) ||
    length(k) != 1L ||
    is.na(k) ||
    !is.finite(k) ||
    k != round(k)
  ) {
    stop(
      "k must be a single integer.",
      call. = FALSE
    )
  }

  k <- as.integer(k)

  if (k < 1L || k > n) {
    stop(
      "k must be between 1 and the number of observations.",
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

  if (
    !is.logical(merge_wrap) ||
    length(merge_wrap) != 1L ||
    is.na(merge_wrap)
  ) {
    stop(
      "merge_wrap must be TRUE or FALSE.",
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

  boundaries <- arrangement$boundaries
  orderings <- arrangement$orderings

  if (!is.numeric(boundaries) || length(boundaries) < 2L) {
    stop(
      "The arrangement contains invalid boundaries.",
      call. = FALSE
    )
  }

  if (
    anyNA(boundaries) ||
    any(!is.finite(boundaries)) ||
    any(diff(boundaries) <= 0)
  ) {
    stop(
      "The arrangement boundaries must be finite and increasing.",
      call. = FALSE
    )
  }

  if (!is.matrix(orderings) || ncol(orderings) != n) {
    stop(
      paste0(
        "The number of columns of arrangement$orderings must equal ",
        "the number of rows of Y."
      ),
      call. = FALSE
    )
  }

  number_before <- nrow(orderings)

  if (
    number_before < 1L ||
    length(boundaries) != number_before + 1L
  ) {
    stop(
      "The arrangement dimensions are inconsistent.",
      call. = FALSE
    )
  }

  upper_rank <- n - k + 1L

  # Only these two columns of the complete ordering matrix are needed at
  # level k.
  lower_indices <- as.integer(
    orderings[, k, drop = TRUE]
  )

  upper_indices <- as.integer(
    orderings[, upper_rank, drop = TRUE]
  )

  alpha_left <- boundaries[-length(boundaries)]
  alpha_right <- boundaries[-1L]

  # Store only the three distinct entries of y_i y_i^T.
  outer_components <- cbind(
    z11 = Y[, 1L]^2,
    z12 = Y[, 1L] * Y[, 2L],
    z22 = Y[, 2L]^2
  )

  outer_products_equal <- function(first_index, second_index) {
    if (first_index == second_index) {
      return(TRUE)
    }

    first <- outer_components[first_index, ]
    second <- outer_components[second_index, ]

    comparison_scale <- max(
      1,
      abs(first),
      abs(second)
    )

    max(abs(first - second)) <= tol * comparison_scale
  }

  same_constraints <- function(first_position, second_position) {
    outer_products_equal(
      lower_indices[first_position],
      lower_indices[second_position]
    ) &&
      outer_products_equal(
        upper_indices[first_position],
        upper_indices[second_position]
      )
  }

  # Find maximal consecutive runs with identical constraints. At this stage
  # only integer indices and numeric vectors are processed.
  group_starts <- integer(number_before)
  group_ends <- integer(number_before)

  group_count <- 1L
  current_start <- 1L

  if (number_before >= 2L) {
    for (cell_index in 2:number_before) {
      contiguous <-
        abs(
          alpha_right[cell_index - 1L] -
            alpha_left[cell_index]
        ) <= tol

      constraints_match <- same_constraints(
        current_start,
        cell_index
      )

      merged_width <-
        alpha_right[cell_index] -
        alpha_left[current_start]

      can_merge <-
        contiguous &&
        constraints_match &&
        merged_width < pi - tol

      if (!can_merge) {
        group_starts[group_count] <- current_start
        group_ends[group_count] <- cell_index - 1L

        group_count <- group_count + 1L
        current_start <- cell_index
      }
    }
  }

  group_starts[group_count] <- current_start
  group_ends[group_count] <- number_before

  group_starts <- group_starts[seq_len(group_count)]
  group_ends <- group_ends[seq_len(group_count)]

  construct_Bj <- function(alpha_left, alpha_right) {
    u_left <- c(
      cos(alpha_left),
      sin(alpha_left)
    )

    u_right <- c(
      cos(alpha_right),
      sin(alpha_right)
    )

    cbind(u_left, u_right)
  }

  construct_outer_product <- function(index) {
    components <- outer_components[index, ]

    matrix(
      c(
        components["z11"],
        components["z12"],
        components["z12"],
        components["z22"]
      ),
      nrow = 2L,
      byrow = TRUE
    )
  }

  construct_level_cell <- function(
    start_position,
    end_position) {
    lower_index <- lower_indices[start_position]
    upper_index <- upper_indices[start_position]

    left_endpoint <- alpha_left[start_position]
    right_endpoint <- alpha_right[end_position]

    list(
      alpha_j = left_endpoint,
      alpha_j1 = right_endpoint,
      lj = lower_index,
      hj = upper_index,
      Zlj = construct_outer_product(lower_index),
      Zhj = construct_outer_product(upper_index),
      Bj = construct_Bj(
        left_endpoint,
        right_endpoint
      )
    )
  }

  # Complete list objects are allocated only for the merged cells.
  merged_cells <- lapply(
    seq_len(group_count),
    function(group_index) {
      construct_level_cell(
        group_starts[group_index],
        group_ends[group_index]
      )
    }
  )

  # The first and last cells are adjacent on the projective circle.
  if (merge_wrap && length(merged_cells) >= 2L) {
    first_cell <- merged_cells[[1L]]
    last_index <- length(merged_cells)
    last_cell <- merged_cells[[last_index]]

    wrap_constraints_match <- same_constraints(
      group_starts[1L],
      group_starts[group_count]
    )

    wrap_width <-
      pi +
      first_cell$alpha_j1 -
      last_cell$alpha_j

    if (
      wrap_constraints_match &&
      wrap_width < pi - tol
    ) {
      # Represent the wrapped projective cell as an unwrapped interval whose
      # right endpoint is larger than pi.
      wrapped_cell <- last_cell
      wrapped_cell$alpha_j1 <-
        pi + first_cell$alpha_j1

      wrapped_cell$Bj <- construct_Bj(
        wrapped_cell$alpha_j,
        wrapped_cell$alpha_j1
      )

      if (length(merged_cells) == 2L) {
        merged_cells <- list(wrapped_cell)
      } else {
        middle_cells <- merged_cells[
          2L:(last_index - 1L)
        ]

        merged_cells <- c(
          middle_cells,
          list(wrapped_cell)
        )
      }
    }
  }

  number_after <- length(merged_cells)

  if (verbose) {
    message(
      sprintf(
        paste0(
          "Reduced the number of level-%d cells from %d to %d ",
          "(%.1f%% reduction)."
        ),
        k,
        number_before,
        number_after,
        100 * (1 - number_after / number_before)
      )
    )
  }

  attr(merged_cells, "number_before") <- number_before
  attr(merged_cells, "number_after") <- number_after
  attr(merged_cells, "n") <- n
  attr(merged_cells, "k") <- k

  merged_cells
}
