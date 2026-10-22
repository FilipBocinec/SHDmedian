test_that("shared centering supports dimensions one, three and five", {
  skip_if_not_installed("clarabel")
  for (d in c(1L, 3L, 5L)) {
    set.seed(100 + d)
    X <- matrix(rnorm(20 * d), 20, d)
    automatic <- SHDmedian(X, center_method = "mean", U = diag(d))
    supplied <- SHDmedian(X, center = colMeans(X), U = diag(d))
    expect_s3_class(automatic, "SHDmedian")
    expect_identical(automatic$method, "approx")
    expect_identical(automatic$d, d)
    expect_equal(automatic$center, colMeans(X))
    expect_equal(automatic$Sigma, supplied$Sigma)
    expect_equal(automatic$search_bounds, supplied$search_bounds)
    expect_equal(automatic$depth, automatic$depth_U)
    expect_false(automatic$certified)
  }
})

test_that("Tukey centering is shared and a supplied center takes precedence", {
  skip_if_not_installed("clarabel")
  expect_equal(SHDmedian(matrix(c(-3, 0, 1)), U = matrix(1))$center, 0)
  set.seed(9)
  X <- matrix(rnorm(60), 20, 3)
  set.seed(8)
  expected <- shd_center(X)$center
  set.seed(8)
  fit <- SHDmedian(X, U = diag(3))
  expect_equal(fit$center, expected)
  supplied <- SHDmedian(X, center = c(0, 1, 2),
                        center_method = "ignored", U = diag(3))
  expect_equal(supplied$center, c(0, 1, 2))
  expect_identical(supplied$center_method, "supplied")
})

test_that("method dispatch rejects unsupported dimensions and domains", {
  expect_error(SHDmedian(matrix(1, 4, 3), method = "exact"), "two columns")
  expect_error(SHDmedian(matrix(1, 4, 2), method = "exact", domain = "pd"),
               "only domain")
  expect_error(SHDmedian(matrix(1), method = "other"), "arg")
  expect_error(SHDmedian(matrix(1), center = Inf), "center")
  expect_error(SHDmedian(matrix(1), verbose = NA), "verbose")
})
