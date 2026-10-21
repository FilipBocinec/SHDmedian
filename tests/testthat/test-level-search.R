test_that("single-level search finds level one", {
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

  result <- search_levels_d2_cutting_plane(
    arrangement = arrangement,
    Y = Y,
    k_min = 1L,
    k_max = 1L,
    solver = "CLARABEL",
    initial_directions = 4L,
    cuts_per_iteration = 4L,
    max_iterations = 50L,
    verbose = FALSE
  )

  expect_equal(result$status, "optimal")
  expect_true(result$certified)
  expect_equal(result$best_k, 1L)
  expect_equal(result$depth, 1 / nrow(Y))
  expect_false(is.null(result$Sigma))
  expect_equal(result$tested_levels, 1L)
  expect_equal(result$number_of_tested_levels, 1L)
})


test_that("binary search returns a certified deepest level", {
  skip_if_not_installed("CVXR")
  skip_if_not_installed("clarabel")

  Y <- rbind(
    c(-1.1, 0.1),
    c(1.0, 0.3),
    c(0.4, -1.2),
    c(0.1, 1.5),
    c(-0.8, 0.7),
    c(0.7, -0.4)
  )

  arrangement <- angular_arrangement_d2(Y)

  result <- search_levels_d2_cutting_plane(
    arrangement = arrangement,
    Y = Y,
    k_min = 1L,
    k_max = 3L,
    solver = "CLARABEL",
    initial_directions = 6L,
    cuts_per_iteration = 5L,
    max_iterations = 100L,
    separation_tol = 1e-8,
    verbose = FALSE
  )

  expect_equal(result$status, "optimal")
  expect_true(result$certified)

  expect_true(
    result$best_k >= 1L &&
      result$best_k <= 3L
  )

  expect_equal(
    result$depth,
    result$best_k / nrow(Y)
  )

  expect_false(is.null(result$Sigma))
  expect_true(isTRUE(result$best_result$feasible))
  expect_true(isTRUE(result$best_result$certified))

  expect_true(all(
    result$tested_levels >= 1L &
      result$tested_levels <= 3L
  ))

  expect_equal(
    length(unique(result$tested_levels)),
    length(result$tested_levels)
  )
})


test_that("level search rejects invalid bounds", {
  Y <- rbind(
    c(-1, 0),
    c(1, 0),
    c(0, -1),
    c(0, 1)
  )

  arrangement <- angular_arrangement_d2(Y)

  expect_error(
    search_levels_d2_cutting_plane(
      arrangement,
      Y,
      k_min = 0,
      k_max = 2
    ),
    "1 <= k_min"
  )

  expect_error(
    search_levels_d2_cutting_plane(
      arrangement,
      Y,
      k_min = 3,
      k_max = 2
    ),
    "1 <= k_min"
  )

  expect_error(
    search_levels_d2_cutting_plane(
      arrangement,
      Y,
      k_min = 1.5,
      k_max = 2
    ),
    "finite integer scalar"
  )
})


test_that("level search output contains diagnostics", {
  skip_if_not_installed("CVXR")
  skip_if_not_installed("clarabel")

  Y <- rbind(
    c(-1, 0),
    c(1, 0.2),
    c(0, -1),
    c(0.3, 1)
  )

  arrangement <- angular_arrangement_d2(Y)

  result <- search_levels_d2_cutting_plane(
    arrangement,
    Y,
    k_min = 1L,
    k_max = 1L,
    initial_directions = 3L,
    max_iterations = 30L,
    verbose = FALSE
  )

  expect_named(
    result,
    c(
      "method",
      "status",
      "certified",
      "search_bounds",
      "best_k",
      "depth",
      "Sigma",
      "best_result",
      "unresolved_k",
      "tested_levels",
      "number_of_tested_levels",
      "level_results",
      "elapsed_seconds"
    )
  )

  expect_true(is.list(result$level_results))
  expect_true(result$elapsed_seconds >= 0)
})
