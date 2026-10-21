# Column-major, sqrt(2)-scaled symmetric vectorization (Clarabel convention).
shd_svec <- function(S) {
  d <- nrow(S)
  index <- upper.tri(S, diag = TRUE)
  weights <- matrix(sqrt(2), d, d)
  diag(weights) <- 1
  (S * weights)[index]
}

shd_unsvec <- function(x, d) {
  S <- matrix(0, d, d)
  S[upper.tri(S, diag = TRUE)] <- x
  S[upper.tri(S)] <- S[upper.tri(S)] / sqrt(2)
  S + t(S) - diag(diag(S), d)
}

# Transparent reference kernels retained for independent checks.
shd_prepare_R <- function(Y, U) {
  Q <- (Y %*% U)^2
  if (any(!is.finite(Q))) stop("Squared projections overflowed.")
  Q <- vapply(seq_len(ncol(U)), function(j) sort(Q[, j]), numeric(nrow(Y)))
  dim(Q) <- c(nrow(Y), ncol(U))
  p <- nrow(U) * (nrow(U) + 1) / 2
  C <- t(matrix(vapply(seq_len(ncol(U)), function(j) shd_svec(tcrossprod(U[, j])),
                numeric(p)), p, ncol(U)))
  list(Q = Q, C = C)
}

shd_verify_R <- function(Sigma, U, Q, tolerance = 0) {
  values <- colSums(U * (Sigma %*% U))
  lower <- upper <- expanded_lower <- expanded_upper <- numeric(ncol(U))
  for (j in seq_len(ncol(U))) {
    lower[j] <- sum(Q[, j] <= values[j])
    upper[j] <- sum(Q[, j] >= values[j])
    expanded_lower[j] <- sum(Q[, j] <= values[j] + tolerance)
    expanded_upper[j] <- sum(Q[, j] >= values[j] - tolerance)
  }
  list(values = values, lower = lower, upper = upper,
       expanded_lower = expanded_lower, expanded_upper = expanded_upper,
       min_eigenvalue = min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values))
}

shd_bounds <- function(Q, k) {
  if (k == 0L) return(list(lower = rep(0, ncol(Q)), upper = rep(Inf, ncol(Q))))
  list(lower = Q[k, ], upper = Q[nrow(Q) - k + 1L, ])
}

shd_search <- function(n, solve_level, incumbent) {
  low <- 0L
  high <- n + 1L
  levels <- list()
  while (high > low + 1L) {
    k <- floor((low + high) / 2)
    answer <- solve_level(k, FALSE)
    levels[[length(levels) + 1L]] <- c(list(k = k, retry = FALSE), answer)
    if (answer$status == "unresolved") {
      answer <- solve_level(k, TRUE)
      levels[[length(levels) + 1L]] <- c(list(k = k, retry = TRUE), answer)
    }
    if (answer$status == "feasible") {
      low <- k
      incumbent <- answer$Sigma
    } else if (answer$status == "infeasible") {
      high <- k
    } else break
  }
  list(Sigma = incumbent, selected_k = low, lower = low, upper = high - 1L,
       status = if (high == low + 1L) "complete" else "partial_resolution",
       levels = levels)
}
