test_that("angular arrangement has consistent dimensions", {
  Y <- rbind(
    c(-1.0, 0.0),
    c(1.0, 0.2),
    c(0.3, -1.0),
    c(0.2, 1.4)
  )

  arrangement <- angular_arrangement_d2(Y)

  expect_s3_class(
    arrangement,
    "shd_angular_arrangement_d2"
  )

  expect_equal(
    length(arrangement$boundaries),
    arrangement$number_of_cells + 1L
  )

  expect_equal(
    nrow(arrangement$orderings),
    arrangement$number_of_cells
  )

  expect_equal(
    ncol(arrangement$orderings),
    nrow(Y)
  )

  expect_length(
    arrangement$has_ties,
    arrangement$number_of_cells
  )
})


test_that("arrangement boundaries cover the projective circle", {
  Y <- rbind(
    c(-1.0, 0.0),
    c(1.0, 0.2),
    c(0.3, -1.0),
    c(0.2, 1.4)
  )

  arrangement <- angular_arrangement_d2(Y)

  expect_equal(
    arrangement$boundaries[1L],
    0,
    tolerance = 1e-15
  )

  expect_equal(
    tail(arrangement$boundaries, 1L),
    pi,
    tolerance = 1e-15
  )

  expect_true(all(diff(arrangement$boundaries) > 0))
})


test_that("stored orderings agree with midpoint projections", {
  Y <- rbind(
    c(-1.0, 0.0),
    c(1.0, 0.2),
    c(0.3, -1.0),
    c(0.2, 1.4),
    c(-0.7, 0.8)
  )

  arrangement <- angular_arrangement_d2(Y)

  for (cell in seq_len(arrangement$number_of_cells)) {
    alpha_left <- arrangement$boundaries[cell]
    alpha_right <- arrangement$boundaries[cell + 1L]
    theta_mid <- (alpha_left + alpha_right) / 2

    direction <- c(
      cos(theta_mid),
      sin(theta_mid)
    )

    squared_projections <- drop(Y %*% direction)^2

    expected_order <- order(
      squared_projections,
      seq_len(nrow(Y)),
      method = "radix"
    )

    expect_equal(
      arrangement$orderings[cell, ],
      expected_order
    )
  }
})


test_that("ordering remains constant inside every cell", {
  Y <- rbind(
    c(-1.0, 0.0),
    c(1.0, 0.2),
    c(0.3, -1.0),
    c(0.2, 1.4),
    c(-0.7, 0.8)
  )

  arrangement <- angular_arrangement_d2(Y)

  for (cell in seq_len(arrangement$number_of_cells)) {
    alpha_left <- arrangement$boundaries[cell]
    alpha_right <- arrangement$boundaries[cell + 1L]

    theta_first <-
      0.75 * alpha_left +
      0.25 * alpha_right

    theta_second <-
      0.25 * alpha_left +
      0.75 * alpha_right

    direction_first <- c(
      cos(theta_first),
      sin(theta_first)
    )

    direction_second <- c(
      cos(theta_second),
      sin(theta_second)
    )

    order_first <- order(
      drop(Y %*% direction_first)^2,
      seq_len(nrow(Y)),
      method = "radix"
    )

    order_second <- order(
      drop(Y %*% direction_second)^2,
      seq_len(nrow(Y)),
      method = "radix"
    )

    expect_equal(order_first, order_second)
    expect_equal(
      arrangement$orderings[cell, ],
      order_first
    )
  }
})


test_that("degenerate observations are handled", {
  Y <- rbind(
    c(0, 0),
    c(1, 0),
    c(1, 0),
    c(-1, 0),
    c(0, 1)
  )

  arrangement <- angular_arrangement_d2(Y)

  expect_true(arrangement$number_of_cells >= 1L)
  expect_true(all(is.finite(arrangement$boundaries)))
  expect_true(any(arrangement$has_ties))
})


test_that("invalid arrangement arguments are rejected", {
  expect_error(
    angular_arrangement_d2(matrix(1:9, ncol = 3L)),
    "exactly two columns"
  )

  Y <- matrix(1:8, ncol = 2L)

  expect_error(
    angular_arrangement_d2(Y, tol = -1),
    "finite nonnegative"
  )

  expect_error(
    angular_arrangement_d2(Y, tol = c(1e-12, 1e-10)),
    "finite nonnegative"
  )
})
