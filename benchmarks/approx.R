# Run from the package root after installing SHDmedian. Optional first argument:
# output directory. R_LIBS_USER may point at a local dependency library.
library(SHDmedian)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else "benchmarks/results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rows <- list()
fits <- list()
for (seed in c(11L, 29L, 47L)) {
  set.seed(seed)
  X <- matrix(rnorm(500L * 10L), 500L, 10L)
  for (reuse in c(TRUE, FALSE)) {
    fit <- SHDmedian_approx(X, center = rep(0, 10), n_directions = 2000L,
      seed = seed + 100L, control = list(reuse = reuse, solver = list(
        verbose = FALSE, presolve_enable = FALSE, chordal_decomposition_enable = FALSE)))
    key <- paste(seed, reuse, sep = "_")
    fits[[key]] <- fit
    rows[[key]] <- data.frame(seed = seed, n = 500, d = 10, L = 2000,
      reuse = reuse, status = fit$status, selected_k = fit$selected_k,
      attained_k = fit$attained_k, lower = fit$search_bounds[1], upper = fit$search_bounds[2],
      calls = fit$solver$calls, reused = fit$solver$reused,
      preprocessing = fit$timing$preprocessing, model_setup = fit$timing$model_setup,
      solver_setup = fit$timing$solver_setup, feasibility_solves = fit$timing$feasibility_solves,
      verification = fit$timing$verification, total = fit$timing$total)
    print(rows[[key]])
  }
  persistent <- fits[[paste(seed, TRUE, sep = "_")]]
  oneshot <- fits[[paste(seed, FALSE, sep = "_")]]
  level_signature <- function(fit) lapply(fit$levels, function(x)
    x[c("k", "retry", "status")])
  stopifnot(identical(persistent$search_bounds, oneshot$search_bounds),
            identical(level_signature(persistent), level_signature(oneshot)))
}
write.csv(do.call(rbind, rows), file.path(out,"timings.csv"), row.names = FALSE)
saveRDS(fits, file.path(out,"fits.rds"))

# Small 2D comparison uses the existing full angular oracle for every k.
set.seed(71)
X <- matrix(rnorm(16), 8, 2)
approx <- SHDmedian_approx(X, c(0,0), n_directions = 100, seed = 72, domain = "psd")
exact_time <- system.time(exact <- SHDmedian(X, mu = c(0,0), initial_directions = 8,
                                           cuts_per_iteration = 8, verbose = FALSE))
arr <- SHDmedian:::angular_arrangement_d2(X)
full <- vapply(seq_len(nrow(X)), function(k) {
  cells <- SHDmedian:::prepare_level_cells_d2(arr, X, k)
  SHDmedian:::separate_level_cells_d2(approx$Sigma, cells, separation_tol = 0)$feasible
}, logical(1))
comparison <- data.frame(n = nrow(X), L = ncol(approx$U),
  approximate_status = approx$status, approximate_depth_U = approx$depth_U,
  approximate_full_angular_depth = max(c(0, which(full))) / nrow(X),
  exact_status = exact$status, exact_depth = exact$depth,
  approximate_seconds = approx$timing$total, exact_seconds = unname(exact_time["elapsed"]))
print(comparison)
write.csv(comparison, file.path(out,"comparison-2d.csv"), row.names = FALSE)
saveRDS(list(approx = approx, exact = exact), file.path(out,"comparison-2d.rds"))
capture.output({
  print(sessionInfo())
  print(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                     "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")))
  print(extSoftVersion())
  print(vapply(c("SHDmedian", "clarabel", "Matrix", "Rcpp", "RcppArmadillo",
                 "CVXR", "mrfDepth"), function(p) as.character(packageVersion(p)), character(1)))
  print(fits[[1L]]$solver)
}, file = file.path(out,"environment.txt"))
