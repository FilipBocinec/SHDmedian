test_that("direction cut represents the projected scatter", {
  Y <- rbind(
    c(-1.0, 0.0),
    c(1.0, 0.2),
    c(0.3, -1.0),
    c(0.2, 1.4)
  )

  arrangement <- angular_arrangement_d2(Y)

  level_cells <- prepare_level_cells_d2(
    arrangement,
    Y,
    k = 1L
  )

  cell <- level_cells[[1L]]
  theta <- 0.5 * (cell$alpha_j + cell$alpha_j1)

  cut <- make_direction_cut_d2(
    cell = cell,
    cell_index = 1L,
    theta = theta,
    side = "lower"
  )

  Sigma <- matrix(
    c(
      2, 0.4,
      0.4, 1
    ),
    nrow = 2L,
    byrow = TRUE
  )

  coefficient_vector <- c(
    s11 = Sigma[1L, 1L],
    s12 = Sigma[1L, 2L],
    s22 = Sigma[2L, 2L]
  )

  represented_projection <- sum(
    cut$coefficients *
      coefficient_vector
  )

  direction <- c(cos(theta), sin(theta))

  direct_projection <- drop(
    crossprod(
      direction,
      Sigma %*% direction
    )
  )

  expect_equal(
    represented_projection,
    direct_projection,
    tolerance = 1e-14
  )

  expect_equal(
    cut$right_hand_side,
    drop(crossprod(
      direction,
      cell$Zlj %*% direction
    )),
    tolerance = 1e-14
  )
})


test_that("duplicate direction cuts are detected", {
  Y <- rbind(
    c(-1, 0),
    c(1, 0),
    c(0, -1),
    c(0, 1)
  )

  arrangement <- angular_arrangement_d2(Y)

  level_cells <- prepare_level_cells_d2(
    arrangement,
    Y,
    k = 1L
  )

  cell <- level_cells[[1L]]
  theta <- 0.5 * (cell$alpha_j + cell$alpha_j1)

  lower_cut <- make_direction_cut_d2(
    cell,
    cell_index = 1L,
    theta = theta,
    side = "lower"
  )

  upper_cut <- make_direction_cut_d2(
    cell,
    cell_index = 1L,
    theta = theta,
    side = "upper"
  )

  expect_true(
    direction_cut_is_duplicate_d2(
      lower_cut,
      list(lower_cut)
    )
  )

  expect_false(
    direction_cut_is_duplicate_d2(
      upper_cut,
      list(lower_cut)
    )
  )
})


test_that("initial cuts contain both constraint sides", {
  Y <- rbind(
    c(-1.0, 0.0),
    c(1.0, 0.2),
    c(0.3, -1.0),
    c(0.2, 1.4),
    c(-0.7, 0.8)
  )

  arrangement <- angular_arrangement_d2(Y)

  level_cells <- prepare_level_cells_d2(
    arrangement,
    Y,
    k = 1L
  )

  cuts <- initialize_direction_cuts_d2(
    level_cells,
    number_of_directions = 4L
  )

  expect_true(length(cuts) >= 2L)
  expect_true(all(vapply(
    cuts,
    function(cut) {
      cut$side %in% c("lower", "upper")
    },
    logical(1L)
  )))

  expect_true(all(vapply(
    cuts,
    function(cut) {
      identical(cut$source, "initial")
    },
    logical(1L)
  )))

  cuts_data <- direction_cuts_to_data_frame_d2(cuts)

  expect_equal(nrow(cuts_data), length(cuts))
})
