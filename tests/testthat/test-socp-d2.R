test_that("SOCP matches independent geometry and depth references", {
skip_if_not_installed("clarabel")
# Independent, deliberately slow full-arrangement reference for small samples.
# It sorts ALL projections in each global cell, unlike the implementation's
# per-observation rank sweep. This is only used in tests.
reference_slacks <- function(Y, Sigma, k) {
  n <- nrow(Y)
  angles <- 0
  for (i in seq_len(n - 1L)) for (j in (i + 1L):n) for (sgn in c(-1, 1)) {
    v <- Y[i, ] + sgn * Y[j, ]
    if (any(v != 0)) angles <- c(angles, atan2(v[1L], -v[2L]) %% pi)
  }
  angles <- sort(unique(angles))
  ends <- c(angles[-1L], angles[1L] + pi)
  answer <- c(lower = Inf, upper = Inf)
  min_on_arc <- function(M, left, right) {
    theta <- ((atan2(2 * M[1L, 2L], M[1L, 1L] - M[2L, 2L]) + pi) / 2) %% pi
    theta <- left + (theta - left) %% pi
    candidates <- c(left, right, if (theta <= right) theta)
    U <- rbind(cos(candidates), sin(candidates))
    min(colSums(U * (M %*% U)))
  }
  for (j in seq_along(angles)) {
    u <- c(cos((angles[j] + ends[j]) / 2), sin((angles[j] + ends[j]) / 2))
    perm <- order(drop(Y %*% u)^2)
    lo <- Y[perm[k], ]
    hi <- Y[perm[n - k + 1L], ]
    answer[1L] <- min(answer[1L], min_on_arc(Sigma - tcrossprod(lo), angles[j], ends[j]))
    answer[2L] <- min(answer[2L], min_on_arc(tcrossprod(hi) - Sigma, angles[j], ends[j]))
  }
  answer
}

set.seed(81427)
samples <- list(
  rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -3)),
  rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1),
        c(3, 0), c(-3, 0), c(0, 3), c(0, -3)),
  rbind(c(0, 0), c(1, 0), c(1, 0), c(0, 2), c(0, -2), c(1, 1)),
  cbind(-2:2, 2 * (-2:2)))
for (n in 3:12) samples[[length(samples) + 1L]] <- matrix(rnorm(2 * n), n, 2)

checks <- 0L
for (Y in samples) {
  p <- SHDprepare_socp(Y, mu = c(0, 0), cache_geometry = TRUE)
  streaming <- SHDprepare_socp(Y, mu = c(0, 0), cache_geometry = FALSE)
  Sigmas <- list(diag(0.6, 2L), crossprod(matrix(rnorm(4), 2L)) + diag(0.1, 2L))
  for (k in seq_len(nrow(Y))) {
    model <- shd_socp_model_cpp(p$ptr, k, p$scale)
    model2 <- shd_socp_model_cpp(streaming$ptr, k, streaming$scale)
    stopifnot(identical(model$structural_empty, model2$structural_empty))
    if (!model$structural_empty)
      stopifnot(identical(model$i, model2$i), identical(model$j, model2$j),
                identical(model$x, model2$x), identical(model$b, model2$b))
    for (S in Sigmas) {
      actual <- shd_socp_validate_cpp(model$level_ptr, S)
      ref <- reference_slacks(Y, S, k)
      stopifnot(max(abs(c(actual$lower_slack, actual$upper_slack) - ref)) < 1e-8)
      oracle <- SHDdepth_d2(Y, S, tolerance = 1e-8)
      if (abs(min(ref) + 1e-8) > 1e-10)
        stopifnot(identical(oracle$k >= k, min(ref) >= -1e-8))
      checks <- checks + 1L
    }
  }
}
cat("Passed", checks, "point-geometry/full-arrangement and depth-oracle comparisons.\n")

# Native CLARABEL integration, ties, and singular-boundary behavior.
two <- rbind(c(1, 0), c(0, 1))
p <- SHDprepare_socp(two, mu = c(0, 0))
one <- SHDlevel_socp(p, 1L)
stopifnot(one$state == "feasible", one$depth_k == 1L,
          max(abs(one$Sigma - diag(0.5, 2L))) < 1e-6)

eight <- samples[[2L]]
p <- SHDprepare_socp(eight, mu = c(0, 0))
four <- SHDlevel_socp(p, 4L, solver_tol = 1e-10)
stopifnot(four$state == "feasible", four$depth_k == 4L,
          max(abs(four$Sigma - diag(0.9, 2L))) < 2e-6)
fit <- SHDmedian_socp(eight, mu = c(0, 0))
stopifnot(fit$k == 4L, fit$certified)

curved <- samples[[1L]]
p <- SHDprepare_socp(curved, mu = c(0, 0))
singular <- SHDlevel_socp(p, 2L, solver_tol = 1e-11, pd_tol = 1e-8)
stopifnot(singular$state %in% c("unknown", "infeasible"))
fit <- SHDmedian_socp(curved, mu = c(0, 0), solver_tol = 1e-11, pd_tol = 1e-8)
stopifnot(fit$k == 1L)
if (singular$state == "unknown") stopifnot(!fit$certified, fit$k_interval[2L] == 2L)

# A curved boundary point: lower equality on a non-endpoint direction.
c0 <- 0.3
S <- matrix(c(c0^2, c0, c0, 2), 2L)
stopifnot(SHDdepth_d2(curved, S, tolerance = 1e-8)$k == 1L)
line <- cbind(-3:3, 2 * (-3:3))
stopifnot(SHDmedian_socp(line, mu = c(0, 0))$k == 0L)

# Median level and validity survive observation permutations and rotations;
# selected matrices need not coincide if a median region has multiple members.
set.seed(3501)
Y <- matrix(rnorm(60), 30L, 2L)
original <- SHDmedian_socp(Y, mu = c(0, 0))
permuted <- SHDmedian_socp(Y[sample.int(nrow(Y)), ], mu = c(0, 0))
theta <- 0.37
Q <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2L)
rotated <- SHDmedian_socp(7 * Y %*% Q, mu = c(0, 0))
stopifnot(original$certified, permuted$certified, rotated$certified,
          original$k == permuted$k, original$k == rotated$k)
cat("Passed native solver, known median, degeneracy, permutation, and rotation/scale tests.\n")
cat("All tests passed.\n")

expect_true(checks > 0L)
})
