# Developer benchmarks: use the installed package, including internal diagnostics.
SHDmedian_socp <- getFromNamespace("SHDmedian_socp", "SHDmedian")
SHDprepare_socp <- getFromNamespace("SHDprepare_socp", "SHDmedian")
SHDlevel_socp <- getFromNamespace("SHDlevel_socp", "SHDmedian")
SHDdepth_d2 <- getFromNamespace("SHDdepth_d2", "SHDmedian")
.shd_data <- getFromNamespace(".shd_data", "SHDmedian")
.shd_clock <- getFromNamespace(".shd_clock", "SHDmedian")
# Run the timing experiment from this directory:
# Rscript benchmark_socp.R
# Fixed center mu=0; location estimation and C++ compilation are excluded.
library(SHDmedian)

invisible(SHDlevel_socp(SHDprepare_socp(rbind(c(1, 0), c(0, 1)), mu = c(0, 0)), 1L))
sizes <- c(100L, 250L, 500L, 1000L)
results <- list()
for (n in sizes) for (rep in 1:3) {
  seed <- 31000L + n + rep
  set.seed(seed)
  X <- matrix(rnorm(2 * n), n, 2L)
  invisible(gc())
  fit <- SHDmedian_socp(X, mu = c(0, 0))
  results[[length(results) + 1L]] <- data.frame(n = n, repetition = rep, seed = seed,
    k = fit$k, strict_k = fit$strict_k, bracket_low = fit$k_interval[1L],
    bracket_high = fit$k_interval[2L], status = fit$status,
    geometry_seconds = unname(fit$timing["geometry"]),
    build_seconds = unname(fit$timing["build"]),
    solver_seconds = unname(fit$timing["solver"]),
    validation_seconds = unname(fit$timing["validation"]),
    total_seconds = unname(fit$timing["total"]), cache_mb = fit$cache_mb,
    feasibility_tests = nrow(fit$levels),
    max_lower_soc = if (nrow(fit$levels)) max(fit$levels$lower_soc) else 0L,
    max_multipliers = if (nrow(fit$levels)) max(fit$levels$multipliers) else 0L,
    max_upper_cuts = if (nrow(fit$levels)) max(fit$levels$upper_cuts) else 0L)
  print(fit)
  flush.console()
}
out <- do.call(rbind, results)
write.csv(out, "benchmark_results.csv", row.names = FALSE)
capture.output(sessionInfo(), file = "benchmark_session.txt")
print(aggregate(total_seconds ~ n, out, median))

