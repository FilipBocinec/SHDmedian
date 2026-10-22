shd_approx <- function(X, center, U = NULL, n_directions = 2000L,
                             seed = NULL, domain = c("pd", "psd"), control = list()) {
  total_start <- proc.time()[["elapsed"]]
  X <- as.matrix(X)
  if (!is.numeric(X) || is.complex(X) || length(dim(X)) != 2L || any(dim(X) == 0L) ||
      any(!is.finite(X))) stop("X must be a nonempty finite numeric matrix.")
  d <- ncol(X)
  n <- nrow(X)
  if (missing(center) || !is.numeric(center) || is.complex(center) || length(center) != d ||
      any(!is.finite(center))) stop("Supply a finite numeric center of length ncol(X).")
  domain <- match.arg(domain)
  defaults <- list(abs_tol = 1e-8, rel_tol = 1e-8, pd_tol = 1e-8,
                   reuse = TRUE, kernel = "cpp", solver = list(verbose = FALSE))
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      any(!names(control) %in% names(defaults)) || anyDuplicated(names(control)))))
    stop("control must be a named list of supported settings.")
  ctrl <- utils::modifyList(defaults, control)
  for (name in c("abs_tol", "rel_tol", "pd_tol")) {
    value <- ctrl[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0)
      stop(name, " must be finite and positive.")
  }
  if (!is.logical(ctrl$reuse) || length(ctrl$reuse) != 1L || is.na(ctrl$reuse))
    stop("reuse must be TRUE or FALSE.")
  ctrl$kernel <- match.arg(ctrl$kernel, c("cpp", "R"))
  if (!is.list(ctrl$solver)) stop("solver settings must be a list.")
  if (is.null(U)) {
    if (!is.numeric(n_directions) || length(n_directions) != 1L ||
        !is.finite(n_directions) || n_directions != floor(n_directions) ||
        n_directions < d || n_directions > .Machine$integer.max)
      stop("n_directions must be an integer at least ncol(X).")
    if (!is.null(seed)) {
      if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
          seed < 0 || seed > .Machine$integer.max || seed != floor(seed))
        stop("seed must be a nonnegative integer.")
      had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
      if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv)
      on.exit(if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv) else
        if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
          rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
      set.seed(seed)
    }
    U <- cbind(diag(d), matrix(stats::rnorm(d * (n_directions - d)), d))
  }
  U <- as.matrix(U)
  if (!is.numeric(U) || is.complex(U) || nrow(U) != d || ncol(U) < 1L || any(!is.finite(U)))
    stop("U must be a finite numeric d by L matrix with L >= 1.")
  maxima <- apply(abs(U), 2L, max)
  if (any(maxima == 0)) stop("U contains a zero direction.")
  U <- sweep(U, 2L, maxima, "/")
  U <- sweep(U, 2L, sqrt(colSums(U^2)), "/")
  Y <- sweep(X, 2L, center, "-")
  if (any(!is.finite(Y))) stop("Centering overflowed; rescale X and center.")
  scale <- max(abs(Y))
  if (scale == 0) scale <- 1
  scatter_scale <- scale^2
  if (!is.finite(scatter_scale) || scatter_scale == 0)
    stop("Scatter scale is outside representable range; rescale X and center.")
  Y <- Y / scale
  if (!requireNamespace("clarabel", quietly = TRUE))
    stop("Install the optional solver with install.packages('clarabel') to use SHDmedian(method = 'approx').")
  solver_defaults <- clarabel::clarabel_control()
  if (length(ctrl$solver) && (is.null(names(ctrl$solver)) ||
      anyDuplicated(names(ctrl$solver)) || any(!names(ctrl$solver) %in% names(solver_defaults))))
    stop("solver must contain named Clarabel control settings.")
  for (name in names(ctrl$solver)) {
    value <- ctrl$solver[[name]]
    default <- solver_defaults[[name]]
    if (length(value) != 1L || anyNA(value) ||
        (is.logical(default) && !is.logical(value)) ||
        (is.numeric(default) && (!is.numeric(value) || value < 0 ||
          (!is.finite(value) && name != "time_limit"))) ||
        (is.integer(default) && value != floor(value)))
      stop("Invalid Clarabel control setting: ", name)
  }
  ctrl$solver <- do.call(clarabel::clarabel_control, ctrl$solver)
  prep <- if (ctrl$kernel == "R") shd_prepare_R(Y, U) else shd_prepare_cpp(Y, U)
  verify <- if (ctrl$kernel == "R") shd_verify_R else shd_verify_cpp
  prep_seconds <- proc.time()[["elapsed"]] - total_start
  setup_start <- proc.time()[["elapsed"]]
  cap <- max(1, max(prep$Q))
  if (!is.finite(cap * scatter_scale))
    stop("Scatter cap is outside representable range; rescale X and center.")
  model <- shd_conic_model(prep$C, d, domain, cap)
  adapter <- shd_solver_adapter(model, ctrl$solver, ctrl$reuse)
  setup_seconds <- proc.time()[["elapsed"]] - setup_start
  verification_seconds <- 0
  tol <- ctrl$abs_tol + ctrl$rel_tol * max(1, max(prep$Q))
  solve_level <- function(k, retry = FALSE) {
    bounds <- shd_bounds(prep$Q, k)
    if (any(bounds$lower > bounds$upper))
      return(list(status = "infeasible", reason = "order_statistic_contradiction"))
    # A zero upper bound in a nonzero direction precludes any PD matrix.
    if (domain == "pd" && any(bounds$upper == 0))
      return(list(status = "infeasible", reason = "zero_upper_bound_precludes_pd"))
    b <- c(bounds$upper, -bounds$lower,
           if (domain == "pd") c(cap, 0), rep(0, model$p))
    answer <- adapter$run(b, retry)
    status_names <- names(clarabel::solver_status_descriptions())
    status <- if (is.numeric(answer$status) && length(answer$status) == 1L &&
                  answer$status >= 1 && answer$status <= length(status_names))
      status_names[answer$status] else as.character(answer$status)
    result <- list(status = "unresolved", solver_status = status, error = answer$error)
    if (identical(status, "PrimalInfeasible")) {
      certificate <- shd_certificate(answer$z, model, b, tol)
      result$certificate <- certificate
      if (isTRUE(certificate$valid)) result$status <- "infeasible"
      return(result)
    }
    if (length(answer$x) < model$p || any(!is.finite(answer$x))) return(result)
    start <- proc.time()[["elapsed"]]
    S <- shd_unsvec(answer$x[seq_len(model$p)], d)
    # Equality rows can be snapped by an independently validated linear correction.
    # This is useful for exact ties; never ridge or clip eigenvalues.
    equal <- which(bounds$lower == bounds$upper)
    if (length(equal)) {
      E <- prep$C[equal, , drop = FALSE]
      sv <- svd(E)
      keep <- which(sv$d > max(sv$d) * 1e-12)
      x <- shd_svec(S)
      if (length(keep)) {
        delta <- sv$v[, keep, drop = FALSE] %*%
          (crossprod(sv$u[, keep, drop = FALSE], bounds$lower[equal] - E %*% x) / sv$d[keep])
        S <- shd_unsvec(x + as.numeric(delta), d)
      }
    }
    check <- verify(S, U, prep$Q, tol)
    violation <- max(0, bounds$lower - check$values, check$values - bounds$upper)
    attained <- min(check$lower, check$upper)
    margin <- if (domain == "pd") utils::tail(answer$x, 1L) else NA_real_
    domain_ok <- if (domain == "pd")
      check$min_eigenvalue > ctrl$pd_tol && margin > ctrl$pd_tol &&
        margin <= cap + tol && check$min_eigenvalue >= margin - tol else
          check$min_eigenvalue >= -tol
    if (violation <= tol && domain_ok && attained >= k) result$status <- "feasible"
    result$reason <- if (attained < k && violation <= tol) "boundary_ambiguity" else
      if (!domain_ok) "domain_margin_unresolved" else "candidate_validation"
    result$Sigma <- S
    result$attained_k <- attained
    result$min_eigenvalue <- check$min_eigenvalue
    result$margin <- margin
    result$max_constraint_violation <- violation
    verification_seconds <<- verification_seconds + proc.time()[["elapsed"]] - start
    result
  }
  search <- shd_search(n, solve_level, diag(cap / 2, d))
  start <- proc.time()[["elapsed"]]
  Sigma <- search$Sigma * scatter_scale
  if (any(!is.finite(Sigma))) stop("Scatter output overflowed; rescale X and center.")
  output_eigenvalue <- min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)
  if (domain == "pd" && output_eigenvalue <= 0)
    stop("Scatter output lost positive definiteness at this scale; rescale X and center.")
  # Multiplication back to original units can round or underflow entries.
  # Verify the actual returned matrix, expressed again in scaled coordinates.
  check <- verify(Sigma / scatter_scale, U, prep$Q, tol)
  if (min(check$lower, check$upper) < search$selected_k) {
    verification_seconds <- verification_seconds + proc.time()[["elapsed"]] - start
    refinement <- solve_level(search$selected_k, TRUE)
    start <- proc.time()[["elapsed"]]
    search$levels[[length(search$levels) + 1L]] <- c(
      list(k = search$selected_k, retry = TRUE, output_refinement = TRUE), refinement)
    if (refinement$status == "feasible") {
      proposed <- refinement$Sigma * scatter_scale
      proposed_check <- verify(proposed / scatter_scale, U, prep$Q, tol)
      proposed_eigenvalue <- min(eigen(proposed, symmetric = TRUE, only.values = TRUE)$values)
      if (min(proposed_check$lower, proposed_check$upper) >= search$selected_k &&
          (domain == "psd" || proposed_eigenvalue > 0)) {
        Sigma <- proposed
        check <- proposed_check
        output_eigenvalue <- proposed_eigenvalue
      }
    }
    if (min(check$lower, check$upper) < search$selected_k)
      search$status <- "boundary_ambiguity"
  }
  bounds <- shd_bounds(prep$Q, search$selected_k)
  violation <- max(0, bounds$lower - check$values, check$values - bounds$upper)
  verification_seconds <- verification_seconds + proc.time()[["elapsed"]] - start
  # Map level diagnostics back to original scatter units as well.
  levels <- lapply(search$levels, function(level) {
    for (name in c("Sigma", "min_eigenvalue", "margin", "max_constraint_violation"))
      if (!is.null(level[[name]])) level[[name]] <- level[[name]] * scatter_scale
    level
  })
  solver <- adapter$diagnostics()
  structure(list(method = "finite-direction SDP", Sigma = Sigma,
    center = as.numeric(center), U = U, domain = domain, n = n,
    requested_k = if (length(levels)) utils::tail(levels, 1L)[[1L]]$k else 0L,
    selected_k = search$selected_k, attained_k = min(check$lower, check$upper),
    depth_U = min(check$lower, check$upper) / n,
    search_bounds = c(lower = search$lower, upper = search$upper), status = search$status,
    min_eigenvalue = output_eigenvalue,
    max_constraint_violation = violation * scatter_scale,
    raw_counts = cbind(lower = check$lower, upper = check$upper),
    tolerance_counts = cbind(lower = check$expanded_lower, upper = check$expanded_upper),
    levels = levels, tolerances = list(absolute = ctrl$abs_tol * scatter_scale,
      relative = ctrl$rel_tol, directional = tol * scatter_scale,
      pd_margin = ctrl$pd_tol * scatter_scale),
    scaling = list(data_divisor = scale, scatter_multiplier = scatter_scale,
                   t_cap = cap * scatter_scale), solver = solver,
    timing = list(preprocessing = prep_seconds, model_setup = setup_seconds,
      solver_setup = solver$setup_seconds, feasibility_solves = solver$solve_seconds,
      verification = verification_seconds, total = proc.time()[["elapsed"]] - total_start)),
    class = "SHDmedian")
}
