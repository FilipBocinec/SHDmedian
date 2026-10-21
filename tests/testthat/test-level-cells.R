test_that("level cells have the expected structure", {
  Y <- rbind(
    c(-1.0, 0.0),
    c(1.0, 0.2),
    c(0.3, -1.0),
    c(0.2, 1.4),
    c(-0.7, 0.8)
  )

  arrangement <- angular_arrangement_d2(Y)

  level_cells <- prepare_level_cells_d2(
    arrangement = arrangement,
    Y = Y,
    k = 2L
  )

  expect_true(length(level_cells) >= 1L)

  expect_equal(
    attr(level_cells, "number_before"),
    arrangement$number_of_cells
  )

  expect_equal(
    attr(level_cells, "number_after"),
    length(level_cells)
  )

  expect_equal(attr(level_cells, "n"), nrow(Y))
  expect_equal(attr(level_cells, "k"), 2L)

  for (cell in level_cells) {
    expect_named(
      cell,
      c(
        "alpha_j",
        "alpha_j1",
        "lj",
        "hj",
        "Zlj",
        "Zhj",
        "Bj"
      )
    )

    expect_true(cell$alpha_j1 > cell$alpha_j)
    expect_lte(cell$alpha_j1 - cell$alpha_j, pi)

    expect_equal(dim(cell$Zlj), c(2L, 2L))
    expect_equal(dim(cell$Zhj), c(2L, 2L))
    expect_equal(dim(cell$Bj), c(2L, 2L))

    expect_equal(
      cell$Zlj,
      tcrossprod(Y[cell$lj, ])
    )

    expect_equal(
      cell$Zhj,
      tcrossprod(Y[cell$hj, ])
    )
  }
})


test_that("merged cells represent every original angular cell", {
  Y <- rbind(
    c(-1.0, 0.0),
    c(1.0, 0.2),
    c(0.3, -1.0),
    c(0.2, 1.4),
    c(-0.7, 0.8),
    c(0.6, 0.4)
  )

  k <- 2L
  arrangement <- angular_arrangement_d2(Y)

  level_cells <- prepare_level_cells_d2(
    arrangement = arrangement,
    Y = Y,
    k = k
  )

  upper_rank <- nrow(Y) - k + 1L

  for (original_cell in seq_len(arrangement$number_of_cells)) {
    alpha_left <- arrangement$boundaries[original_cell]
    alpha_right <- arrangement$boundaries[original_cell + 1L]
    theta_mid <- (alpha_left + alpha_right) / 2

    covered <- vapply(
      level_cells,
      function(cell) {
        candidate_angles <- c(
          theta_mid,
          theta_mid + pi
        )

        any(
          candidate_angles > cell$alpha_j &
            candidate_angles < cell$alpha_j1
        )
      },
      logical(1L)
    )

    expect_equal(sum(covered), 1L)

    merged_cell <- level_cells[[which(covered)]]

    lower_index <-
      arrangement$orderings[original_cell, k]

    upper_index <-
      arrangement$orderings[original_cell, upper_rank]

    expect_equal(
      merged_cell$Zlj,
      tcrossprod(Y[lower_index, ]),
      tolerance = 1e-12
    )

    expect_equal(
      merged_cell$Zhj,
      tcrossprod(Y[upper_index, ]),
      tolerance = 1e-12
    )
  }
})


test_that("cell merging never increases the number of cells", {
  Y <- rbind(
    c(0, 0),
    c(1, 0),
    c(1, 0),
    c(-1, 0),
    c(0, 1),
    c(0, -1)
  )

  arrangement <- angular_arrangement_d2(Y)

  level_cells <- prepare_level_cells_d2(
    arrangement = arrangement,
    Y = Y,
    k = 2L
  )

  expect_lte(
    attr(level_cells, "number_after"),
    attr(level_cells, "number_before")
  )
})


test_that("invalid level-cell arguments are rejected", {
  Y <- rbind(
    c(-1, 0),
    c(1, 0),
    c(0, -1),
    c(0, 1)
  )

  arrangement <- angular_arrangement_d2(Y)

  expect_error(
    prepare_level_cells_d2(arrangement, Y, k = 0),
    "between 1"
  )

  expect_error(
    prepare_level_cells_d2(arrangement, Y, k = 5),
    "between 1"
  )

  expect_error(
    prepare_level_cells_d2(arrangement, Y, k = 1.5),
    "single integer"
  )

  expect_error(
    prepare_level_cells_d2(
      arrangement,
      Y,
      k = 1,
      tol = -1
    ),
    "finite nonnegative"
  )

  expect_error(
    prepare_level_cells_d2(
      list(),
      Y,
      k = 1
    ),
    "angular_arrangement_d2"
  )
})
