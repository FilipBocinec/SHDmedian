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
  expect_error(SHDmedian(matrix(1, 4, 2), method = "exact", domain = "psd"),
               "only domain")
  expect_error(SHDmedian(matrix(1), method = "other"), "arg")
  expect_error(SHDmedian(matrix(1), center = Inf), "center")
  expect_error(SHDmedian(matrix(1), verbose = NA), "verbose")
})

test_that("exact SOCP dispatch preserves diagnostics and numerical uncertainty", {
  skip_if_not_installed("clarabel")
  X <- rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -3))
  fit <- SHDmedian(X, method = "exact", center = c(0, 0),
    control = list(cache_geometry = FALSE, solver_tol = 1e-11, pd_tol = 1e-8))
  expect_equal(fit$algorithm, "point-based-socp")
  expect_equal(fit$domain, "pd")
  expect_equal(fit$best_k, 1L)
  expect_equal(unname(fit$search_bounds), fit$search_details$k_interval)
  expect_false(fit$search_details$cache_geometry)
  if (!fit$certified) expect_equal(fit$status, "unresolved")
  expect_error(SHDmedian(X, method = "exact", center = c(0, 0),
    control = list(unknown = TRUE)), "control")
  expect_error(SHDmedian(X, method = "exact", center = c(0, 0),
    solver = "SCS"), "CLARABEL")
  expect_error(SHDmedian(X, method = "exact", center = c(0, 0),
    control = list(solver_tol = "bad")), "positive")
})
