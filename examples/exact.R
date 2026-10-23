library(SHDmedian)

set.seed(42)
X <- matrix(rnorm(100), ncol = 2)
fit <- SHDmedian(X, method = "exact", center_method = "mean",
  control = list(cache_geometry = "auto", solver_tol = 1e-9))
print(fit)
fit$Sigma
fit$search_bounds
fit$search_details$strict_depth
fit$search_details$levels
