# Developer benchmarks: use the installed package, including internal diagnostics.
SHDmedian_socp <- getFromNamespace("SHDmedian_socp", "SHDmedian")
SHDprepare_socp <- getFromNamespace("SHDprepare_socp", "SHDmedian")
SHDlevel_socp <- getFromNamespace("SHDlevel_socp", "SHDmedian")
SHDdepth_d2 <- getFromNamespace("SHDdepth_d2", "SHDmedian")
.shd_data <- getFromNamespace(".shd_data", "SHDmedian")
.shd_clock <- getFromNamespace(".shd_clock", "SHDmedian")
# Source this file to define the comparison helper.

.shd_extract_matrix <- function(fit) {
  if (is.matrix(fit) && identical(dim(fit), c(2L, 2L))) return(fit)
  for (name in c("Sigma", "sigma", "scatter", "median", "Sigma_hat", "Sigma_median")) {
    z <- fit[[name]]
    if (is.matrix(z) && identical(dim(z), c(2L, 2L))) return(z)
  }
  stop("Cannot find the cutting-plane matrix. Supply cutting_matrix = function(fit) fit$YOUR_FIELD.")
}

# X and the SAME precomputed mu go to both methods. No package source is
# modified. Supply an adapter if your installed package uses another API.
compare_shd_methods <- function(X, mu, cutting_fun = NULL,
                                cutting_matrix = .shd_extract_matrix,
                                repetitions = 3L, warmup = TRUE,
                                socp_args = list(), validation_tol = 1e-7,
                                verbose = TRUE) {
  X <- .shd_data(X)
  if (missing(mu) || length(mu) != 2L || any(!is.finite(mu)))
    stop("Supply one precomputed two-dimensional mu shared by both methods.")
  if (repetitions < 1L || repetitions != as.integer(repetitions)) stop("Invalid repetitions.")
  if (is.null(cutting_fun)) {
    if (!requireNamespace("scatterHD", quietly = TRUE) ||
        !"SHDmedian" %in% getNamespaceExports("scatterHD"))
      stop("Supply cutting_fun = function(X, mu) YOUR_EXISTING_SOLVER(X, mu = mu).")
    cutting_fun <- function(X, mu) scatterHD::SHDmedian(X, mu = mu, solver = "CLARABEL")
  }
  
  socp <- function(X, mu) do.call(SHDmedian_socp,
                                  c(list(X = X, mu = mu), socp_args))
  methods <- list(point_socp = socp, cutting_plane = cutting_fun)
  if (warmup) {
    if (verbose) message("One untimed warm-up per method, on the same data.")
    for (method in methods) tryCatch(method(X, mu), error = function(e)
      warning("Warm-up failed: ", conditionMessage(e)))
  }
  Y <- sweep(X, 2L, mu, "-")
  norms <- rowSums(Y * Y)
  scale <- if (any(norms > 0)) median(norms[norms > 0]) else 1
  rows <- list()
  fits <- list()
  for (rep in seq_len(repetitions)) {
    # Alternate order without consuming or changing the user's RNG state.
    method_names <- if (rep %% 2L) names(methods) else rev(names(methods))
    for (name in method_names) {
      invisible(gc())
      started <- .shd_clock()
      fit <- tryCatch(methods[[name]](X, mu), error = identity)
      elapsed <- .shd_clock() - started
      key <- paste(name, rep, sep = "_")
      fits[[key]] <- fit
      row <- data.frame(method = name, repetition = rep, n = nrow(X),
                        elapsed = elapsed, verified_k = NA_integer_, strict_k = NA_integer_,
                        min_eigenvalue = NA_real_, reported_k = NA_real_,
                        numerically_resolved = NA, status = "", message = "")
      if (inherits(fit, "error")) {
        row$status <- "error"
        row$message <- conditionMessage(fit)
      } else {
        S <- tryCatch(if (name == "point_socp") fit$Sigma else cutting_matrix(fit), error = identity)
        if (inherits(S, "error")) {
          row$status <- "extraction_error"
          row$message <- conditionMessage(S)
        } else {
          # Independent validation is outside the timed call for BOTH methods.
          verified <- SHDdepth_d2(X, S, mu, tolerance = validation_tol * scale)
          strict <- SHDdepth_d2(X, S, mu, tolerance = 0)
          row$verified_k <- verified$k
          row$strict_k <- strict$k
          row$min_eigenvalue <- min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
          row$status <- if (name == "point_socp") fit$status else "returned"
          if (is.list(fit) && length(fit$k) == 1L && is.numeric(fit$k)) row$reported_k <- fit$k
          if (is.list(fit) && length(fit$certified) == 1L) row$numerically_resolved <- fit$certified
        }
      }
      rows[[length(rows) + 1L]] <- row
      if (verbose) message(sprintf("%s, rep %d: %.3fs, verified k=%s, %s",
                                    name, rep, elapsed, row$verified_k, row$status))
    }
  }
  runs <- do.call(rbind, rows)
  summaries <- lapply(split(runs, runs$method), function(z) {
    ok <- !is.na(z$verified_k)
    data.frame(method = z$method[1L], successful_runs = sum(ok),
      median_seconds = if (any(ok)) median(z$elapsed[ok]) else NA_real_,
      min_verified_k = if (any(ok)) min(z$verified_k[ok]) else NA_integer_,
      max_verified_k = if (any(ok)) max(z$verified_k[ok]) else NA_integer_)
  })
  list(summary = do.call(rbind, summaries), runs = runs, fits = fits,
       mu = mu, validation_tol = validation_tol, scale = scale,
       note = "Compare times only when verified depth and numerical completion agree; median matrices need not be unique.")
}

# Example, after loading your package and both files:
# mu <- mrfDepth::hdepthmedian(X)$median
# bench <- compare_shd_methods(
#   X, mu,
#   cutting_fun = function(X, mu) scatterHD::SHDmedian(
#     X, mu = mu, solver = "CLARABEL", separation_tol = 1e-8),
#   repetitions = 3L)
# bench$summary
# bench$runs
# write.csv(bench$runs, "comparison.csv", row.names = FALSE)

