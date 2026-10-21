test_that("SHDmedian returns the requested structure", {
  skip_if_not_installed("CVXR")
  skip_if_not_installed("clarabel")

  X <- rbind(
    c(-1.0, 0.0),
    c(1.0, 0.2),
    c(0.3, -1.0),
    c(0.2, 1.4)
  )

  result <- SHDmedian(
    X = X,
    center_method = "mean",
    solver = "CLARABEL",
    initial_directions = 4L,
    cuts_per_iteration = 4L,
    max_iterations = 100L,
    verbose = FALSE
  )

  expect_s3_class(result, "SHDmedian")

  expect_named(
    result,
    c(
      "method",
      "status",
      "certified",
      "n",
      "center",
      "Sigma",
      "best_k",
      "depth",
      "k_min",
      "k_max",
      "timing",
      "search_details"
    )
  )

  expect_equal(result$method, "cutting-plane")
  expect_equal(result$status, "optimal")
  expect_true(result$certified)
  expect_equal(result$n, nrow(X))
  expect_equal(result$center, colMeans(X))

  expect_equal(
    result$depth,
    result$best_k / result$n
  )

  expect_equal(
    result$k_min,
    ceiling(nrow(X) / 4)
  )

  expect_equal(
    result$k_max,
    as.integer(ceiling(nrow(X) / 2))
  )
  expect_false(is.null(result$Sigma))

  expect_equal(
    result$Sigma,
    t(result$Sigma),
    tolerance = 1e-8
  )

  eigenvalues <- eigen(
    result$Sigma,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  expect_gte(min(eigenvalues), -1e-6)

  expect_named(
    result$timing,
    c(
      "center",
      "arrangement",
      "preprocessing",
      "search",
      "total"
    )
  )

  expect_true(result$timing$total >= 0)
  expect_true(is.list(result$search_details))
})


test_that("a supplied center is used by SHDmedian", {
  skip_if_not_installed("CVXR")
  skip_if_not_installed("clarabel")

  X <- rbind(
    c(-1.0, 0.0),
    c(1.0, 0.2),
    c(0.3, -1.0),
    c(0.2, 1.4)
  )

  supplied_center <- c(0.1, -0.2)

  result <- SHDmedian(
    X = X,
    center_method = "tukey",
    mu = supplied_center,
    initial_directions = 4L,
    cuts_per_iteration = 4L,
    max_iterations = 100L
  )

  expect_equal(
    result$center,
    supplied_center
  )
})


test_that("SHDmedian rejects unsupported data", {
  expect_error(
    SHDmedian(matrix(1:9, ncol = 3L)),
    "exactly two columns"
  )

  expect_error(
    SHDmedian(matrix(1:2, nrow = 1L)),
    "at least two observations"
  )

  X <- matrix(1:8, ncol = 2L)

  expect_error(
    SHDmedian(
      X,
      initial_directions = 0
    ),
    "positive integer"
  )

  expect_error(
    SHDmedian(
      X,
      separation_tol = -1
    ),
    "finite nonnegative"
  )
})


test_that("SHDmedian has a print method", {
  skip_if_not_installed("CVXR")
  skip_if_not_installed("clarabel")

  X <- rbind(
    c(-1, 0),
    c(1, 0.2),
    c(0, -1),
    c(0.3, 1)
  )

  result <- SHDmedian(
    X,
    center_method = "mean",
    initial_directions = 3L,
    cuts_per_iteration = 3L,
    max_iterations = 50L
  )

  output <- capture.output(print(result))

  expect_true(any(grepl(
    "Bivariate scatter halfspace median",
    output,
    fixed = TRUE
  )))

  expect_true(any(grepl(
    "Scatter median",
    output,
    fixed = TRUE
  )))
})
