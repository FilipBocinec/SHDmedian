test_that("column svec agrees with quadratic forms beyond dimension two", {
  set.seed(17)
  for (d in c(3, 5)) {
    S <- matrix(rnorm(d*d), d); S <- S + t(S)
    U <- matrix(rnorm(d*7), d)
    expect_equal(shd_unsvec(shd_svec(S), d), S)
    prep <- shd_prepare_R(matrix(rnorm(4*d), 4), U)
    expect_equal(as.numeric(prep$C %*% shd_svec(S)), colSums(U * (S %*% U)))
  }
  expect_equal(shd_svec(matrix(1:9, 3)), c(1, 4*sqrt(2), 5, 7*sqrt(2), 8*sqrt(2), 9))
})

test_that("bounds and raw inclusive counts preserve ties", {
  Q <- cbind(c(0,1,1,4), c(2,2,2,3))
  expect_equal(shd_bounds(Q, 3), list(lower=c(1,2), upper=c(1,2)))
  expect_equal(shd_bounds(Q, 0)$lower, c(0,0))
  check <- shd_verify_R(matrix(1), matrix(1), matrix(c(1,1)), 0)
  expect_equal(check$lower, 2); expect_equal(check$upper, 2)
  check <- shd_verify_R(matrix(1 + 1e-10), matrix(1), matrix(c(1,1)), 1e-8)
  expect_equal(check$upper, 0); expect_equal(check$expanded_upper, 2)
})

test_that("Rcpp kernels match independent R calculations", {
  set.seed(8)
  for (d in c(1,3,5)) {
    Y <- matrix(rnorm(11*d),11); U <- matrix(rnorm(d*13),d)
    R <- shd_prepare_R(Y,U); cpp <- shd_prepare_cpp(Y,U)
    expect_equal(cpp,R,tolerance=1e-12)
    S <- crossprod(matrix(rnorm(d*d),d))
    expect_equal(shd_verify_cpp(S,U,R$Q,1e-8), shd_verify_R(S,U,R$Q,1e-8), tolerance=1e-12)
  }
  Q <- matrix(c(0,1,1,4))
  expect_equal(shd_verify_cpp(matrix(1),matrix(1),Q), shd_verify_R(matrix(1),matrix(1),Q))
})

test_that("unresolved search outcomes preserve both supported bounds", {
  solve <- function(k,retry) list(status=if(k<=2) "feasible" else
    if(k>=6) "infeasible" else "unresolved", Sigma=matrix(k))
  answer <- shd_search(8,solve,matrix(1))
  expect_equal(c(answer$lower,answer$upper),c(0,8))
  expect_length(answer$levels,2)
  expect_equal(answer$status,"partial_resolution")
  for (optimum in 0:8) {
    solve <- function(k,retry) list(status=if(k<=optimum) "feasible" else "infeasible",Sigma=matrix(k))
    answer <- shd_search(8,solve,matrix(1))
    expect_equal(c(answer$lower,answer$upper), rep(optimum,2))
  }
})

test_that("approximate API validates fixed centers and directions", {
  expect_error(SHDmedian(matrix(1,2,3), center = 0),"center")
  expect_error(SHDmedian(matrix(numeric(),0,2), center = c(0,0)),"nonempty")
  expect_error(SHDmedian(matrix(Inf), center = 0),"finite")
  expect_error(SHDmedian(matrix(1+1i), center = 0),"finite")
  expect_error(SHDmedian(matrix(1), center = 1i),"center")
  expect_error(SHDmedian(matrix(1), center = 0,U=matrix(1i)),"finite")
  expect_error(SHDmedian(matrix(1), center = 0,U=matrix(0)),"zero direction")
  expect_error(SHDmedian(matrix(1,2,3), center = rep(0,3),n_directions=2),"at least")
  expect_error(SHDmedian(matrix(1e200), center = 0),"representable")
  expect_error(SHDmedian(matrix(1), center = 0,control=list(abs_tol=NA_real_)),"positive")
  expect_error(SHDmedian(matrix(1), center = 0,control=list(unknown=TRUE)),"supported")
})

test_that("scalar normalization supports representable extreme magnitudes", {
  skip_if_not_installed("clarabel")
  for (scale in c(1e-150,1e150)) {
    a <- SHDmedian(matrix(c(-1,1)*scale), center = 0,U=matrix(1))
    expect_equal(a$Sigma / scale^2,matrix(1),tolerance=1e-12)
    expect_equal(a$depth_U,1)
    expect_true(is.finite(a$tolerances$directional))
  }
  expect_error(SHDmedian(matrix(c(0,1.6e-162),1), center = c(0,0),U=diag(2)),"rescale")
})

test_that("solver handles tied depth above one half and distinguishes PD", {
  skip_if_not_installed("clarabel")
  for (kernel in c("R","cpp")) {
    a <- SHDmedian(matrix(c(-1,1)), center = 0,U=matrix(c(2,-3),1),control=list(kernel=kernel))
    expect_equal(a$depth_U,1); expect_equal(a$selected_k,2)
    expect_equal(a$search_bounds,c(lower=2,upper=2))
    expect_equal(a$Sigma,matrix(1),tolerance=1e-12)
    pd <- SHDmedian(matrix(0,3,2), center = c(0,0),U=diag(2),control=list(kernel=kernel))
    psd <- SHDmedian(matrix(0,3,2), center = c(0,0),U=diag(2),domain="psd",control=list(kernel=kernel))
    expect_equal(pd$selected_k,0); expect_gt(pd$min_eigenvalue,0)
    expect_equal(pd$depth_U,0); expect_equal(pd$status,"complete")
    expect_equal(psd$depth_U,1); expect_equal(psd$Sigma,matrix(0,2,2),tolerance=1e-12)
  }
})

test_that("fixed directions respect translation, scaling, signs and RNG", {
  skip_if_not_installed("clarabel")
  set.seed(51); X <- matrix(rnorm(60),20,3)
  old <- .Random.seed
  a <- SHDmedian(X, center = rep(0,3),n_directions=25,seed=2)
  expect_identical(.Random.seed,old)
  b <- SHDmedian(sweep(X*4,2,c(1,2,3),"+"), center = c(1,2,3),U=a$U)
  expect_equal(b$Sigma,a$Sigma*16,tolerance=1e-6)
  expect_equal(b$search_bounds,a$search_bounds)
  c <- SHDmedian(X, center = rep(0,3),U=-a$U*7,control=list(reuse=FALSE))
  expect_equal(c$search_bounds,a$search_bounds)
  expect_equal(abs(c$U),abs(a$U),tolerance=1e-14)
  check <- shd_verify_R(a$Sigma,a$U,shd_prepare_R(X,a$U)$Q)
  expect_equal(a$attained_k,min(check$lower,check$upper))
  expect_output(print(a),"Finite-direction")
})

test_that("larger direction sets and independent validation obey monotonicity", {
  skip_if_not_installed("clarabel")
  set.seed(10); X <- matrix(rnorm(80),40,2)
  U <- cbind(diag(2),matrix(rnorm(40),2))
  a <- SHDmedian(X, center = c(0,0),U=U[,1:6])
  b <- SHDmedian(X, center = c(0,0),U=U)
  expect_lte(b$search_bounds[1],a$search_bounds[2])
  verification <- shd_verify_R(a$Sigma,U,shd_prepare_R(X,U)$Q)
  expect_lte(min(verification$lower,verification$upper),a$attained_k)
})

test_that("dual certificates are independently checked", {
  model <- shd_conic_model(matrix(1),1,"psd",1)
  # x <= 0, x >= 1, x >= 0: z=(1,1,0) is a Farkas witness.
  expect_true(shd_certificate(c(1,1,0),model,c(0,-1,0),1e-8)$valid)
  expect_false(shd_certificate(c(-1,-1,0),model,c(0,-1,0),1e-8)$valid)
  expect_false(shd_certificate(c(1,0,0),model,c(0,-1,0),1e-8)$valid)
})

test_that("one-dimensional exhaustive feasibility agrees with search", {
  skip_if_not_installed("clarabel")
  X <- matrix(c(-3,-2,-1,1,2,3))
  a <- SHDmedian(X, center = 0,U=matrix(1),domain="psd")
  Q <- shd_prepare_R(X,matrix(1))$Q
  feasible <- vapply(1:6,function(k) {b<-shd_bounds(Q,k); all(b$lower<=b$upper)},logical(1))
  optimum <- max(which(feasible))
  expect_lte(a$search_bounds[1],optimum); expect_gte(a$search_bounds[2],optimum)
  expect_equal(a$status,"complete")
})

test_that("solver iteration limits do not manufacture infeasibility", {
  skip_if_not_installed("clarabel")
  set.seed(22)
  a <- SHDmedian(matrix(rnorm(40),20,2), center = c(0,0),n_directions=12,
    seed=3,control=list(solver=list(max_iter=0L)))
  expect_equal(a$status,"partial_resolution")
  expect_equal(a$search_bounds,c(lower=0,upper=20))
  expect_true(all(vapply(a$levels,function(x) x$status=="unresolved",logical(1))))
})

test_that("a numerically tiny PD margin stays unresolved", {
  skip_if_not_installed("clarabel")
  X <- rbind(c(-1e-6,-1),c(1e-6,1))
  a <- SHDmedian(X, center = c(0,0),U=diag(2))
  expect_equal(a$status,"partial_resolution")
  expect_equal(a$search_bounds,c(lower=0,upper=2))
  expect_true(all(vapply(a$levels,function(x) x$status=="unresolved",logical(1))))
  expect_gt(a$min_eigenvalue,0)
})

test_that("exhaustive small SDP levels agree with the binary search bracket", {
  skip_if_not_installed("clarabel")
  set.seed(19)
  X <- matrix(rnorm(18),9,2)
  U <- cbind(diag(2),c(1,1)/sqrt(2),c(1,-1)/sqrt(2))
  prep <- shd_prepare_R(X,U)
  model <- shd_conic_model(prep$C,2,"psd",1)
  outcomes <- vapply(1:9,function(k) {
    bounds <- shd_bounds(prep$Q,k)
    if (any(bounds$lower>bounds$upper)) return("infeasible")
    b <- c(bounds$upper,-bounds$lower,rep(0,3))
    r <- clarabel::clarabel(model$A,b,model$q,cones=model$cones,
                            control=list(verbose=FALSE,tol_feas=1e-10))
    S <- shd_unsvec(r$x,2)
    check <- shd_verify_R(S,U,prep$Q)
    if (min(check$lower,check$upper)>=k && check$min_eigenvalue>=-1e-8)
      return("feasible")
    if (shd_certificate(r$z,model,b,1e-8)$valid) return("infeasible")
    "unresolved"
  },character(1))
  low <- max(c(0,which(outcomes=="feasible")))
  high <- min(c(10,which(outcomes=="infeasible")))-1
  a <- SHDmedian(X, center = c(0,0),U=U,domain="psd")
  expect_lte(a$search_bounds[1],high)
  expect_gte(a$search_bounds[2],low)
  if (low==high && a$status=="complete") expect_equal(unname(a$search_bounds),c(low,high))
})

test_that("approximate matrix respects the full 2D angular depth comparison", {
  skip_if_not_installed("clarabel")
  skip_if_not_installed("CVXR")
  set.seed(71)
  X <- matrix(rnorm(16),8,2)
  a <- SHDmedian(X, center = c(0,0),n_directions=100,seed=72,domain="psd")
  exact <- SHDmedian(X,method="exact",center=c(0,0),initial_directions=8,cuts_per_iteration=8)
  arr <- angular_arrangement_d2(X)
  passed <- vapply(1:8,function(k)
    separate_level_cells_d2(a$Sigma,prepare_level_cells_d2(arr,X,k),
                            separation_tol=0)$feasible,logical(1))
  full_k <- max(c(0,which(passed)))
  expect_lte(full_k,a$attained_k)
  expect_lte(full_k,exact$best_k)
  expect_equal(exact$status,"optimal")
  # In the same PSD domain, full-direction optimization cannot beat the U optimum.
  expect_lte(exact$best_k,a$search_bounds[2])
})
