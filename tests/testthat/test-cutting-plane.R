test_that("violated oracle checks generate new cuts", {
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

  oracle <- separate_level_cells_d2(
    Sigma = 100 * diag(2L),
    level_cells = level_cells
  )

  expect_false(oracle$feasible)
  expect_true(any(oracle$checks$violated))

  new_cuts <- generate_violated_cut_batch_d2(
    oracle_result = oracle,
    level_cells = level_cells,
    existing_cuts = list(),
    iteration = 1L,
    cuts_per_iteration = 3L
  )

  expect_true(length(new_cuts) >= 1L)
  expect_lte(length(new_cuts), 3L)

  expect_true(all(vapply(
    new_cuts,
    function(cut) {
      identical(cut$source, "oracle")
    },
    logical(1L)
  )))

  expect_true(all(vapply(
    new_cuts,
    function(cut) {
      identical(cut$iteration, 1L)
    },
    logical(1L)
  )))
})


test_that("cutting plane certifies a feasible level", {
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

  result <- solve_level_d2_cutting_plane(
    level_cells = level_cells,
    solver = "CLARABEL",
    initial_directions = 4L,
    cuts_per_iteration = 4L,
    max_iterations = 50L,
    separation_tol = 1e-8,
    verbose = FALSE
  )

  expect_true(result$feasible)
  expect_true(result$certified)
  expect_equal(result$status, "feasible")
  expect_false(is.null(result$Sigma))
  expect_true(result$iterations >= 1L)
  expect_true(result$number_of_cuts >= 2L)

  final_check <- separate_level_cells_d2(
    Sigma = result$Sigma,
    level_cells = level_cells,
    separation_tol = 1e-8
  )

  expect_true(final_check$feasible)
})


test_that("cutting-plane output contains diagnostics", {
  skip_if_not_installed("CVXR")
  skip_if_not_installed("clarabel")

  Y <- rbind(
    c(-1, 0),
    c(1, 0.2),
    c(0, -1),
    c(0.3, 1)
  )

  arrangement <- angular_arrangement_d2(Y)

  level_cells <- prepare_level_cells_d2(
    arrangement,
    Y,
    k = 1L
  )

  result <- solve_level_d2_cutting_plane(
    level_cells = level_cells,
    initial_directions = 3L,
    max_iterations = 30L,
    verbose = FALSE
  )

  expect_named(
    result,
    c(
      "method",
      "feasible",
      "certified",
      "status",
      "solver_status",
      "solver_error",
      "Sigma",
      "eigenvalues",
      "iterations",
      "number_of_cells",
      "number_of_cuts",
      "cuts",
      "cut_table",
      "history",
      "oracle",
      "master",
      "trace_bound",
      "elapsed_seconds"
    )
  )

  expect_equal(result$method, "cutting-plane")
  expect_true(is.data.frame(result$history))
  expect_true(is.data.frame(result$cut_table))
  expect_true(result$elapsed_seconds >= 0)
})


test_that("invalid cutting-plane arguments are rejected", {
  expect_error(
    solve_level_d2_cutting_plane(
      level_cells = list()
    ),
    "at least one cell"
  )

  dummy_cell <- list(
    alpha_j = 0,
    alpha_j1 = pi / 2,
    Zlj = diag(2L),
    Zhj = 2 * diag(2L)
  )

  expect_error(
    solve_level_d2_cutting_plane(
      level_cells = list(dummy_cell),
      max_iterations = 0
    ),
    "positive integer"
  )

  expect_error(
    solve_level_d2_cutting_plane(
      level_cells = list(dummy_cell),
      separation_tol = -1
    ),
    "finite nonnegative"
  )
})
