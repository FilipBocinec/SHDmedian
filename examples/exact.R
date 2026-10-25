library(SHDmedian)

set.seed(42)

X <- matrix(rnorm(50000), ncol = 2)

fit <- SHDmedian(X, method = "approx", center = c(0,0))

print(fit)
cov(X)


fit$Sigma
fit$search_bounds
fit$search_details$strict_depth
fit$search_details$levels
