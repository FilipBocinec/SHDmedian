library(SHDmedian)
library(DepthDescent)
set.seed(100)
X <- matrix(rnorm(4000), 1000, 4)
fit <- SHDmedian(X, center = rep(0, 4), n_directions = 10000, seed = 42)
print(fit)
print(fit$search_bounds)
print(fit$timing)

fit$Sigma/qnorm(3/4)^2
X <- matrix(rnorm(4000), 1000, 4)

comparison <- compare_shd_depthdescent(
  X,
  center = rep(0, ncol(X)),
  n_directions = 1000L,
  ntry = 200L,
  n_validation = 20000L,
  seed_U = 402L,
  seed_validation = 2026L,
  true_cov = diag(ncol(X))
)

print(comparison$table)
print(comparison$sdp_search_bounds)
print(comparison$fit_sdp)



# Independent extra directions: evaluate on the UNION, keeping this matrix fixed.
V <- matrix(rnorm(3 * 100), 3)
V <- sweep(V, 2, sqrt(colSums(V^2)), "/")
U_union <- cbind(fit$U, V)
Q <- (X %*% U_union)^2
thresholds <- colSums(U_union * (fit$Sigma %*% U_union))
counts <- vapply(seq_len(ncol(Q)), function(j)
  min(sum(Q[,j] <= thresholds[j]), sum(Q[,j] >= thresholds[j])), numeric(1))
print(c(training_depth = fit$depth_U, union_depth = min(counts) / nrow(X)))
stopifnot(min(counts) / nrow(X) <= fit$depth_U)
