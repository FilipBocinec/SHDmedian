test_that("trace bound is finite and positive", {
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

  trace_bound <- compute_trace_bound_d2(
    level_cells
  )

  expect_length(trace_bound, 1L)
  expect_true(is.finite(trace_bound))
  expect_gt(trace_bound, 0)
})


test_that("CLARABEL solves a restricted master problem", {
  skip_if_not_installed("CVXR")
  skip_if_not_installed("clarabel")

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

  trace_bound <- compute_trace_bound_d2(
    level_cells
  )

  result <- solve_direction_master_d2(
    cuts = cuts,
    trace_bound = trace_bound,
    solver = "CLARABEL"
  )

  expect_equal(result$status, "optimal")
  expect_null(result$solver_error)
  expect_true(is.finite(result$margin))
  expect_false(is.null(result$Sigma))

  expect_equal(
    result$Sigma,
    t(result$Sigma),
    tolerance = 1e-10
  )

  eigenvalues <- eigen(
    result$Sigma,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  expect_gte(
    min(eigenvalues),
    -1e-6
  )

  expect_lte(
    sum(diag(result$Sigma)),
    trace_bound + 1e-6
  )
})
