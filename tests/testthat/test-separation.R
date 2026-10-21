test_that("quadratic form evaluation agrees with matrix multiplication", {
  Q <- matrix(
    c(
      2, 0.5,
      0.5, 3
    ),
    nrow = 2L,
    byrow = TRUE
  )

  theta <- c(0, 0.3, 1.2, pi)

  expected <- vapply(
    theta,
    function(angle) {
      u <- c(cos(angle), sin(angle))
      drop(t(u) %*% Q %*% u)
    },
    numeric(1L)
  )

  expect_equal(
    quadratic_form_at_angle_d2(Q, theta),
    expected,
    tolerance = 1e-14
  )
})


test_that("arc minimizer finds an interior eigenvector", {
  minimum_angle <- 0.7

  rotation <- matrix(
    c(
      cos(minimum_angle),
      sin(minimum_angle),
      -sin(minimum_angle),
      cos(minimum_angle)
    ),
    nrow = 2L
  )

  Q <-
    rotation %*%
    diag(c(1, 4)) %*%
    t(rotation)

  result <- minimize_quadratic_on_arc_d2(
    Q = Q,
    alpha_left = 0.3,
    alpha_right = 1.1
  )

  expect_equal(
    result$value,
    1,
    tolerance = 1e-12
  )

  expect_equal(
    result$theta,
    minimum_angle,
    tolerance = 1e-12
  )
})


test_that("arc minimizer handles constant quadratic forms", {
  result <- minimize_quadratic_on_arc_d2(
    Q = 2 * diag(2L),
    alpha_left = 0.2,
    alpha_right = 1.3
  )

  expect_equal(
    result$value,
    2,
    tolerance = 1e-14
  )
})


test_that("sample second moment is feasible at level one", {
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
    k = 1L
  )

  Sigma <- crossprod(Y) / nrow(Y)

  result <- separate_level_cells_d2(
    Sigma = Sigma,
    level_cells = level_cells,
    separation_tol = 1e-10
  )

  expect_true(result$feasible)
  expect_true(result$psd_feasible)
  expect_equal(result$number_violated, 0)
  expect_false(any(result$checks$violated))
})


test_that("an excessively large scatter matrix violates upper constraints", {
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
    k = 1L
  )

  result <- separate_level_cells_d2(
    Sigma = 100 * diag(2L),
    level_cells = level_cells
  )

  expect_false(result$feasible)
  expect_true(result$psd_feasible)
  expect_true(result$number_of_cell_violations > 0L)

  expect_true(any(
    result$checks$side == "upper" &
      result$checks$violated
  ))
})


test_that("a non-PSD candidate is detected", {
  Y <- rbind(
    c(-1, 0),
    c(1, 0),
    c(0, -1),
    c(0, 1)
  )

  arrangement <- angular_arrangement_d2(Y)

  level_cells <- prepare_level_cells_d2(
    arrangement = arrangement,
    Y = Y,
    k = 1L
  )

  result <- separate_level_cells_d2(
    Sigma = -diag(2L),
    level_cells = level_cells
  )

  expect_false(result$feasible)
  expect_false(result$psd_feasible)
  expect_lt(result$psd_margin, 0)
})


test_that("invalid separation arguments are rejected", {
  expect_error(
    minimize_quadratic_on_arc_d2(
      Q = diag(3L),
      alpha_left = 0,
      alpha_right = 1
    ),
    "2 x 2"
  )

  expect_error(
    minimize_quadratic_on_arc_d2(
      Q = diag(2L),
      alpha_left = 1,
      alpha_right = 0
    ),
    "must not be smaller"
  )

  expect_error(
    minimize_quadratic_on_arc_d2(
      Q = diag(2L),
      alpha_left = 0,
      alpha_right = 2 * pi
    ),
    "width at most pi"
  )

  expect_error(
    separate_level_cells_d2(
      Sigma = diag(3L),
      level_cells = list(list())
    ),
    "2 x 2"
  )

  expect_error(
    separate_level_cells_d2(
      Sigma = diag(2L),
      level_cells = list()
    ),
    "at least one cell"
  )
})
