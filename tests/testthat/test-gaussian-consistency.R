test_that("Gaussian consistency is enabled by default for both backends", {
  skip_if_not_installed("clarabel")
  X <- rbind(c(-1, 0), c(1, 0.2), c(0.3, -1), c(0.2, 1.4))
  for (method in c("approx", "exact")) {
    raw <- SHDmedian(X, method = method, center = c(0, 0), U = diag(2),
                     gaussian_consistent = FALSE)
    corrected <- SHDmedian(X, method = method, center = c(0, 0), U = diag(2))
    expect_identical(raw$Sigma, raw$Sigma_raw)
    expect_identical(raw$consistency_factor, 1)
    expect_false(raw$gaussian_consistent)
    expect_true(corrected$gaussian_consistent)
    expect_equal(corrected$consistency_factor, 1 / qchisq(0.5, df = 1))
    expect_equal(corrected$Sigma_raw, raw$Sigma)
    expect_equal(corrected$Sigma, raw$Sigma / qnorm(0.75)^2)
    for (field in c("depth", "best_k", "search_bounds", "status", "certified"))
      expect_identical(corrected[[field]], raw[[field]])
  }
})

test_that("Gaussian correction recovers covariance on Gaussian projection medians", {
  skip_if_not_installed("clarabel")
  variance <- 9
  X <- matrix(c(-1, 1) * qnorm(0.75) * sqrt(variance))
  fit <- SHDmedian(X, center = 0, U = matrix(1))
  expect_equal(fit$Sigma, matrix(variance), tolerance = 1e-7)
  expect_equal(fit$Sigma_raw, matrix(variance * qnorm(0.75)^2), tolerance = 1e-7)
})

test_that("Gaussian consistency option requires a single logical value", {
  for (value in list(NA, NULL, logical(), c(TRUE, FALSE), 1, "TRUE"))
    expect_error(SHDmedian(matrix(1), gaussian_consistent = value),
                 "gaussian_consistent must be TRUE or FALSE", fixed = TRUE)
})
