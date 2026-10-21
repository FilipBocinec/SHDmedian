# Construct once; only the RHS changes across levels. C is generally dense.
shd_conic_model <- function(C, d, domain, cap) {
  L <- nrow(C)
  p <- ncol(C)
  if (domain == "psd") {
    A <- rbind(C, -C, -diag(p))
    objective <- rep(0, p)
    linear <- 2L * L
  } else {
    A <- rbind(cbind(C, 0), cbind(-C, 0), c(rep(0, p), 1),
               c(rep(0, p), -1), cbind(-diag(p), shd_svec(diag(d))))
    objective <- c(rep(0, p), -1)
    linear <- 2L * L + 2L
  }
  list(A = Matrix::Matrix(A, sparse = TRUE), q = objective,
       cones = list(l = as.integer(linear), s = as.integer(d)),
       linear = linear, d = d, p = p, domain = domain, cap = cap)
}

# Validate a normalized Farkas witness independently of solver status.
shd_certificate <- function(z, model, b, tol) {
  if (length(z) != length(b) || any(!is.finite(z)) || max(abs(z)) == 0)
    return(list(valid = FALSE))
  z <- z / max(abs(z))
  residual <- max(abs(as.numeric(crossprod(model$A, z)))) /
    (1 + max(Matrix::colSums(abs(model$A))))
  cone_min <- min(z[seq_len(model$linear)],
                  eigen(shd_unsvec(z[-seq_len(model$linear)], model$d),
                        symmetric = TRUE, only.values = TRUE)$values)
  separation <- sum(b * z) / (1 + sum(abs(b)))
  list(valid = residual <= tol && cone_min >= -tol && separation < -tol,
       residual = residual, cone_min = cone_min, separation = separation)
}

shd_solver_adapter <- function(model, settings, reuse = TRUE) {
  solver <- NULL
  calls <- 0L
  reused <- 0L
  setup_seconds <- solve_seconds <- 0
  exports <- getNamespaceExports("clarabel")
  persistent <- reuse && all(c("clarabel_solver", "solver_is_update_allowed",
                               "solver_update", "solver_solve") %in% exports)
  run <- function(b, retry = FALSE) {
    ctrl <- settings
    if (retry) {
      ctrl$tol_feas <- min(1e-10, if (is.null(ctrl$tol_feas)) 1e-10 else ctrl$tol_feas)
      ctrl$tol_gap_abs <- ctrl$tol_gap_rel <- 1e-10
    }
    args <- list(A = model$A, b = b, q = model$q, P = NULL,
                 cones = model$cones, control = ctrl)
    # A zero-objective retry seeks an interior witness instead of placing a
    # maximum-margin PD solution on a directional boundary. It is used only
    # for feasibility: no conclusion about a zero optimum is drawn from it.
    if (retry && model$domain == "pd") args$q <- rep(0, length(model$q))
    calls <<- calls + 1L
    start <- proc.time()[["elapsed"]]
    setup_elapsed <- 0
    result <- tryCatch({
      if (persistent && !retry) {
        if (is.null(solver)) {
          solver <<- do.call(clarabel::clarabel_solver, args)
          setup_elapsed <- proc.time()[["elapsed"]] - start
          setup_seconds <<- setup_seconds + setup_elapsed
        } else if (clarabel::solver_is_update_allowed(solver)) {
          clarabel::solver_update(solver, b = b)
          reused <<- reused + 1L
        } else {
          solver <<- NULL
          persistent <<- FALSE
        }
        if (is.null(solver)) do.call(clarabel::clarabel, args) else
          clarabel::solver_solve(solver)
      } else do.call(clarabel::clarabel, args)
    }, error = function(e) list(error = conditionMessage(e)))
    solve_seconds <<- solve_seconds + proc.time()[["elapsed"]] - start - setup_elapsed
    result
  }
  list(run = run, diagnostics = function() list(calls = calls, reused = reused,
       persistent_available = persistent, setup_seconds = setup_seconds,
       solve_seconds = solve_seconds, version = as.character(utils::packageVersion("clarabel")),
       settings = settings))
}
