.shd_clock <- function() unname(proc.time()["elapsed"])

.shd_data <- function(X) {
  X <- as.matrix(X)
  if (!is.numeric(X) || ncol(X) != 2L || nrow(X) < 1L || any(!is.finite(X)))
    stop("X must be a finite numeric matrix with two columns and at least one row.")
  storage.mode(X) <- "double"
  X
}
.shd_positive <- function(x, name, zero = FALSE) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L ||
      !is.finite(x) || if (zero) x < 0 else x <= 0)
    stop(name, " must be a finite ", if (zero) "nonnegative" else "positive", " number.")
  invisible(x)
}

# Preparation can be shared by several single-level solves. External pointers
# in this object are session-local: do not saveRDS() and reuse in another R session.
SHDprepare_socp <- function(X, mu = NULL, center_method = c("tukey", "mean"),
                            cache_geometry = "auto", cache_limit_mb = 256,
                            angle_tol = 1e-12) {
  started <- .shd_clock()

  X <- .shd_data(X)
  .shd_positive(angle_tol, "angle_tol", zero = TRUE)
  .shd_positive(cache_limit_mb, "cache_limit_mb")
  center_method <- match.arg(center_method)
  ct <- .shd_clock()
  dithered <- FALSE
  if (is.null(mu)) {
    if (center_method == "mean") mu <- colMeans(X) else {
      if (!requireNamespace("mrfDepth", quietly = TRUE))
        stop("For a Tukey center install mrfDepth, or supply the same mu used by your other solver.")
      location <- mrfDepth::hdepthmedian(X)
      mu <- as.numeric(location$median)
      dithered <- isTRUE(location$dithered)
      if (dithered) warning("mrfDepth reports dithering in the location calculation; supply mu explicitly for reproducible comparisons.")
    }
  } else center_method <- "supplied"
  if (length(mu) != 2L || any(!is.finite(mu))) stop("mu must have two finite components.")
  Y <- sweep(X, 2L, as.numeric(mu), "-")
  norms <- rowSums(Y * Y)
  if (any(!is.finite(norms))) stop("Squared coordinates overflow. Rescale the data by a common scalar.")
  scale <- if (any(norms > 0)) stats::median(norms[norms > 0]) else 1
  ct <- .shd_clock() - ct
  if (identical(cache_geometry, "auto")) cache_geometry <- 24 * nrow(Y)^2 <= cache_limit_mb * 1024^2
  if (!is.logical(cache_geometry) || length(cache_geometry) != 1L || is.na(cache_geometry))
    stop("cache_geometry must be TRUE, FALSE, or 'auto'.")
  # Exact collinearity test on the supplied floating-point coordinates.
  v <- Y[which.max(norms), ]
  collinear <- all(Y[, 1L] * v[2L] - Y[, 2L] * v[1L] == 0)
  gt <- .shd_clock()
  geometry <- shd_socp_prepare_cpp(Y, cache_geometry, angle_tol)
  gt <- .shd_clock() - gt
  structure(list(Y = Y, mu = as.numeric(mu), center_method = center_method,
                 location_dithered = dithered, n = nrow(Y), scale = scale,
                 collinear = collinear, ptr = geometry$ptr,
                 cache_geometry = cache_geometry, cache_bytes = geometry$cache_bytes,
                 event_count = geometry$event_count, angle_tol = angle_tol,
                 timing = c(center = ct, geometry = gt, preparation = .shd_clock() - started)),
            class = "shd_socp_prepared")
}

# Independent depth evaluation: roots of q_Sigma(theta)-r_i(theta), rather
# than the point-pair geometry used to construct the optimization problem.
# tolerance is in squared projection units, and is absolute here.
SHDdepth_d2 <- function(X, Sigma, mu = c(0, 0), tolerance = 0, angle_tol = 1e-13) {

  X <- .shd_data(X)
  Sigma <- as.matrix(Sigma)
  if (!identical(dim(Sigma), c(2L, 2L)) || any(!is.finite(Sigma)) ||
      max(abs(Sigma - t(Sigma))) > 1e-12 * max(1, abs(Sigma)))
    stop("Sigma must be a finite symmetric 2 by 2 matrix.")
  if (length(mu) != 2L || any(!is.finite(mu))) stop("Invalid mu.")
  .shd_positive(tolerance, "tolerance", zero = TRUE)
  .shd_positive(angle_tol, "angle_tol", zero = TRUE)
  storage.mode(Sigma) <- "double"
  shd_depth_point_cpp(sweep(X, 2L, mu, "-"), Sigma, tolerance, angle_tol)
}

.shd_status <- function(solution) {
  s <- solution$status
  if (is.character(s)) return(s[1L])
  descriptions <- clarabel::solver_status_descriptions()
  if (length(s) == 1L && is.finite(s) && s >= 1L && s <= length(descriptions))
    return(names(descriptions)[as.integer(s)])
  paste0("Status", s)
}

# The formulation is exact in real arithmetic. All decisions below are
# numerical, and "unknown" is retained instead of treating solver failure or
# a nearly zero positive-definiteness margin as mathematical infeasibility.
SHDlevel_socp <- function(prepared, k, solver_tol = 1e-9,
                          validation_tol = 1e-7, pd_tol = 1e-9,
                          control = list(), keep_model = FALSE, polish = TRUE) {
  if (!inherits(prepared, "shd_socp_prepared")) stop("Use SHDprepare_socp first.")
  if (!requireNamespace("Matrix", quietly = TRUE) || !requireNamespace("clarabel", quietly = TRUE))
    stop("Install Matrix and clarabel first.")
  if (length(k) != 1L || !is.finite(k) || k != as.integer(k) || k < 1L || k > prepared$n)
    stop("k must be an integer between 1 and n.")
  .shd_positive(solver_tol, "solver_tol")
  .shd_positive(validation_tol, "validation_tol", zero = TRUE)
  .shd_positive(pd_tol, "pd_tol", zero = TRUE)
  started <- .shd_clock()
  model <- shd_socp_model_cpp(prepared$ptr, as.integer(k), prepared$scale)
  ans <- list(k = as.integer(k), state = "unknown", Sigma = NULL,
              solver_status = NA_character_, delta = NA_real_, min_eigenvalue = NA_real_,
              depth = NA_real_, strict_depth = NA_real_, depth_k = NA_integer_,
              strict_k = NA_integer_, lower_slack = NA_real_, upper_slack = NA_real_,
              n_lower_soc = model$n_lower_soc, n_multipliers = model$n_multipliers,
              n_upper_cuts = model$n_upper_cuts, n_lower_arcs = model$n_lower_arcs,
              n_upper_arcs = model$n_upper_arcs, solver_iterations = NA_integer_,
              solver_time = 0, validation_time = 0,
              polish_status = "not_needed", interior_margin = NA_real_)
  if (model$structural_empty) {
    ans$state <- "infeasible"
    ans$solver_status <- "UpperArcForcesSingularity"
    ans$reason <- "An upper-active arc contains a direction orthogonal to its observation."
    ans$build_time <- .shd_clock() - started
  } else {
    A <- Matrix::sparseMatrix(i = model$i, j = model$j, x = model$x,
                             dims = c(model$nrow, model$nvar))
    cones <- list(l = as.integer(model$nlin), q = rep.int(3L, model$nsoc))
    settings <- utils::modifyList(list(verbose = FALSE, max_iter = 250L,
      tol_gap_abs = solver_tol, tol_gap_rel = solver_tol, tol_feas = solver_tol,
      tol_infeas_abs = solver_tol, tol_infeas_rel = solver_tol), control)
    ans$build_time <- .shd_clock() - started
    st <- .shd_clock()
    solution <- tryCatch(clarabel::clarabel(A = A, b = model$b, q = model$q,
                         cones = cones, control = settings), error = identity)
    ans$solver_time <- .shd_clock() - st
    if (inherits(solution, "error")) {
      ans$solver_status <- "Error"
      ans$reason <- conditionMessage(solution)
    } else {
      ans$solver_status <- .shd_status(solution)
      ans$solver_iterations <- solution$iterations
      ans$primal_residual <- solution$r_prim
      ans$dual_residual <- solution$r_dual
      if (ans$solver_status == "PrimalInfeasible") {
        ans$state <- "infeasible"
        ans$reason <- "CLARABEL returned its numerical primal-infeasibility certificate."
      } else if (length(solution$x) == model$nvar && all(is.finite(solution$x))) {
        vt <- .shd_clock()
        x <- solution$x
        Sigma <- prepared$scale * matrix(c(x[1L], x[2L], x[2L], x[3L]), 2L)
        ans$Sigma <- Sigma
        ans$delta <- prepared$scale * x[4L]
        ans$min_eigenvalue <- min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)
        valid <- shd_socp_validate_cpp(model$level_ptr, Sigma)
        ans$lower_slack <- valid$lower_slack
        ans$upper_slack <- valid$upper_slack
        strict <- shd_depth_point_cpp(prepared$Y, Sigma, 0, 1e-13)
        relaxed <- shd_depth_point_cpp(prepared$Y, Sigma,
                                                    validation_tol * prepared$scale, 1e-13)
        ans$depth <- relaxed$depth
        ans$strict_depth <- strict$depth
        ans$depth_k <- relaxed$k
        ans$strict_k <- strict$k
        if (ans$min_eigenvalue > pd_tol * prepared$scale && relaxed$k >= k &&
            min(valid$lower_slack, valid$upper_slack) >= -validation_tol * prepared$scale) {
          ans$state <- "feasible"
          ans$reason <- "Positive definite candidate passed both all-direction numerical checks."
        } else if (ans$min_eigenvalue <= pd_tol * prepared$scale) {
          ans$reason <- "Positive-definiteness margin is at or below pd_tol; strict SPD feasibility is unresolved."
        } else ans$reason <- "Candidate failed all-direction validation at the requested tolerance."
        ans$validation_time <- .shd_clock() - vt
      } else ans$reason <- "No finite candidate or accepted infeasibility certificate was returned."
    }
    # If roundoff puts the maximum-eigenvalue-margin candidate just outside
    # a depth boundary, try a relative interior point. Keep half its proven
    # positive eigenvalue margin fixed and maximize a common lower/upper slack.
    # Failure of this OPTIONAL solve never changes a level to infeasible.
    if (polish && ans$state == "feasible" && ans$strict_k < k) {
      bt <- .shd_clock()
      pd_row <- model$nlin + 1L
      low_rows <- model$nlin + 4L + 3L * (seq_len(model$n_lower_soc) - 1L)
      up_rows <- 1L + model$n_multipliers + seq_len(model$n_upper_cuts)
      keep <- !(model$i == pd_row & model$j == 4L)
      rows <- c(low_rows, up_rows)
      Am <- Matrix::sparseMatrix(i = c(model$i[keep], rows),
             j = c(model$j[keep], rep.int(4L, length(rows))),
             x = c(model$x[keep], rep.int(1, length(rows))),
             dims = c(model$nrow, model$nvar))
      bm <- model$b
      bm[pd_row] <- -ans$min_eigenvalue / prepared$scale
      ans$build_time <- ans$build_time + .shd_clock() - bt
      st <- .shd_clock()
      polished <- tryCatch(clarabel::clarabel(A = Am, b = bm, q = model$q,
                           cones = cones, control = settings), error = identity)
      ans$solver_time <- ans$solver_time + .shd_clock() - st
      ans$polish_status <- if (inherits(polished, "error")) "Error" else .shd_status(polished)
      if (!inherits(polished, "error") && ans$polish_status %in% c("Solved", "AlmostSolved") &&
          length(polished$x) == model$nvar && all(is.finite(polished$x))) {
        vt <- .shd_clock()
        x <- polished$x
        S <- prepared$scale * matrix(c(x[1L], x[2L], x[2L], x[3L]), 2L)
        val <- shd_socp_validate_cpp(model$level_ptr, S)
        dep <- shd_depth_point_cpp(prepared$Y, S, 0, 1e-13)
        emin <- min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
        if (emin > pd_tol * prepared$scale && dep$k >= k &&
            min(val$lower_slack, val$upper_slack) >= -validation_tol * prepared$scale) {
          ans$Sigma <- S
          ans$min_eigenvalue <- emin
          ans$lower_slack <- val$lower_slack
          ans$upper_slack <- val$upper_slack
          ans$strict_k <- dep$k
          ans$strict_depth <- dep$depth
          relaxed <- shd_depth_point_cpp(prepared$Y, S,
                                           validation_tol * prepared$scale, 1e-13)
          ans$depth_k <- relaxed$k
          ans$depth <- relaxed$depth
          ans$interior_margin <- prepared$scale * x[4L]
        }
        ans$validation_time <- ans$validation_time + .shd_clock() - vt
      }
    }
    if (keep_model) {
      ans$model <- list(A = A, b = model$b, q = model$q, cones = cones, scale = prepared$scale)
      ans$solution <- solution
    }
  }
  ans$total_time <- .shd_clock() - started
  ans$certificate_type <- "floating-point, tolerance-based; not a formal exact certificate"
  ans
}

SHDmedian_socp <- function(X, mu = NULL, center_method = c("tukey", "mean"),
                           cache_geometry = "auto", cache_limit_mb = 256,
                           angle_tol = 1e-12, solver_tol = 1e-9,
                           validation_tol = 1e-7, pd_tol = 1e-9,
                           control = list(), verbose = FALSE) {
  .shd_positive(solver_tol, "solver_tol")
  .shd_positive(validation_tol, "validation_tol", zero = TRUE)
  .shd_positive(pd_tol, "pd_tol", zero = TRUE)

  started <- .shd_clock()
  prepared <- SHDprepare_socp(X, mu, match.arg(center_method), cache_geometry,
                              cache_limit_mb, angle_tol)
  n <- prepared$n
  if (prepared$collinear) {
    ans <- list(Sigma = diag(prepared$scale, 2L), mu = prepared$mu, n = n,
                k = 0L, depth = 0, strict_k = 0L, strict_depth = 0,
                k_interval = c(0L, 0L), certified = TRUE,
                status = "rank_deficient_depth_zero", levels = data.frame(),
                center_method = prepared$center_method, scale = prepared$scale,
                timing = c(prepared$timing[c("center", "geometry")],
                           build = 0, solver = 0, validation = 0, total = .shd_clock() - started))
  } else {
    Sigma <- crossprod(prepared$Y) / n
    # This initial matrix always satisfies the k=1 inequalities in exact
    # arithmetic. Evaluate its depth; no uncomputed feasible matrix is assumed.
    initial <- shd_depth_point_cpp(prepared$Y, Sigma,
                                                validation_tol * prepared$scale, 1e-13)
    best_k <- min(initial$k, floor(n / 2))
    search_upper <- upper <- floor(n / 2) + 1L
    history <- list()
    while (search_upper - best_k > 1L) {
      k <- as.integer(floor((best_k + search_upper) / 2))
      result <- SHDlevel_socp(prepared, k, solver_tol, validation_tol, pd_tol, control,
                              polish = FALSE)
      if (verbose) cat(sprintf("k=%d: %s (%s), SOC=%d, lambda=%d, upper=%d, %.3fs\n",
                         k, result$state, result$solver_status, result$n_lower_soc,
                         result$n_multipliers, result$n_upper_cuts, result$total_time))
      history[[length(history) + 1L]] <- data.frame(k = k, state = result$state,
        solver_status = result$solver_status, delta = result$delta,
        min_eigenvalue = result$min_eigenvalue, strict_k = result$strict_k,
        validated_k = result$depth_k, lower_slack = result$lower_slack,
        upper_slack = result$upper_slack, lower_soc = result$n_lower_soc,
        multipliers = result$n_multipliers, upper_cuts = result$n_upper_cuts,
        build = result$build_time, solver = result$solver_time,
        validation = result$validation_time, reason = result$reason)
      if (result$state == "feasible") {
        best_k <- k
        Sigma <- result$Sigma
      } else {
        # On UNKNOWN, search lower levels for a useful candidate but do NOT
        # tighten the reported upper bound. This is not an infeasibility claim.
        search_upper <- k
        if (result$state == "infeasible") upper <- min(upper, k)
      }
    }
    tab <- if (length(history)) do.call(rbind, history) else data.frame()
    strict <- shd_depth_point_cpp(prepared$Y, Sigma, 0, 1e-13)
    polish_time <- c(build = 0, solver = 0, validation = 0)
    polish_status <- "not_needed"
    if (best_k > 0L && strict$k < best_k) {
      refined <- SHDlevel_socp(prepared, as.integer(best_k), solver_tol,
                               validation_tol, pd_tol, control, polish = TRUE)
      polish_time <- c(build = refined$build_time, solver = refined$solver_time,
                        validation = refined$validation_time)
      polish_status <- refined$polish_status
      if (refined$state == "feasible" && refined$strict_k > strict$k) {
        Sigma <- refined$Sigma
        strict <- shd_depth_point_cpp(prepared$Y, Sigma, 0, 1e-13)
      }
    }
    relaxed <- shd_depth_point_cpp(prepared$Y, Sigma,
                                                validation_tol * prepared$scale, 1e-13)
    resolved <- best_k == upper - 1L
    ans <- list(Sigma = Sigma, mu = prepared$mu, n = n, k = best_k,
                depth = best_k / n, strict_k = strict$k, strict_depth = strict$depth,
                validated_k = relaxed$k, k_interval = c(best_k, upper - 1L),
                certified = resolved, status = if (resolved) "resolved_numerically" else "unresolved",
                levels = tab, center_method = prepared$center_method, scale = prepared$scale,
                polish_status = polish_status,
                timing = c(prepared$timing[c("center", "geometry")],
                  build = (if (nrow(tab)) sum(tab$build) else 0) + unname(polish_time["build"]),
                  solver = (if (nrow(tab)) sum(tab$solver) else 0) + unname(polish_time["solver"]),
                  validation = (if (nrow(tab)) sum(tab$validation) else 0) + unname(polish_time["validation"]),
                  total = .shd_clock() - started))
  }
  ans$method <- "Point-based exact SOCP formulation in d=2"
  ans$certificate_type <- "floating-point, tolerance-based; not a formal exact certificate"
  ans$cache_geometry <- prepared$cache_geometry
  ans$cache_mb <- prepared$cache_bytes / 1024^2
  ans$tolerances <- c(angle = angle_tol, solver = solver_tol,
                      validation = validation_tol, pd = pd_tol)
  class(ans) <- "shd_socp_median"
  ans
}
