test_that("mean centering is correct", {
  X <- rbind(
    c(0, 0),
    c(2, 4),
    c(4, 2)
  )

  result <- center_bivariate_data_d2(
    X = X,
    method = "mean"
  )

  expect_equal(result$mu, colMeans(X))
  expect_equal(result$Y, sweep(X, 2L, colMeans(X), "-"))
  expect_equal(result$method, "mean")
  expect_true(is.na(result$center_depth))
})


test_that("a supplied center overrides the method", {
  X <- rbind(
    c(0, 0),
    c(2, 4),
    c(4, 2)
  )

  mu <- c(1, 2)

  result <- center_bivariate_data_d2(
    X = X,
    method = "tukey",
    mu = mu
  )

  expect_equal(result$mu, mu)
  expect_equal(result$Y, sweep(X, 2L, mu, "-"))
  expect_equal(result$method, "supplied")
})


test_that("Tukey centering uses mrfDepth", {
  X <- rbind(
    c(-1, 0),
    c(1, 0),
    c(0, -1),
    c(0, 1),
    c(0, 0)
  )

  expected <- mrfDepth::hdepthmedian(
    x = X,
    maxdir = NULL
  )

  result <- center_bivariate_data_d2(
    X = X,
    method = "tukey"
  )

  expect_equal(
    result$mu,
    as.numeric(expected$median),
    tolerance = 1e-12
  )

  expect_equal(result$method, "tukey")
  expect_equal(dim(result$Y), dim(X))
  expect_true(all(is.finite(result$Y)))
})
