# Finite-direction scatter halfspace median

`SHDmedian_approx()` solves the finite-direction problem for a supplied fixed
center in any positive dimension. It returns one maximizing scatter matrix,
subject to its reported numerical bracket. It does not compute the barycenter
of a median set. The existing bivariate `SHDmedian()` API, centering convention,
cutting-plane optimization, and complete angular oracle are unchanged.

```r
library(SHDmedian)
set.seed(1)
X <- matrix(rnorm(300), 100, 3)
fit <- SHDmedian_approx(X, center = rep(0, 3), n_directions = 200, seed = 42)
fit$Sigma
fit$depth_U
fit$search_bounds
fit$status
```

The fixed center is required, including in dimension two. To use the package's
existing Tukey convention, compute that center first with
`mrfDepth::hdepthmedian(X, maxdir = NULL)$median` for bivariate X and pass it
explicitly. No higher-dimensional location estimator is selected implicitly.
Empty samples and zero-dimensional matrices are rejected; one observation and
dimension one are supported. Inputs must be real finite numeric matrices.

Automatic directions contain the d coordinate axes plus normalized Gaussian
columns, with `n_directions` counting the total. Explicit `U` is a d by L matrix
and is normalized without appending columns. Zero columns are rejected;
repeated directions and sign flips are valid. U stays fixed throughout search
and is returned. `seed` temporarily seeds automatic direction generation and
restores the caller's RNG state; without a seed it consumes the ordinary R RNG
stream. With explicit U there is no random generation and seed is ignored.

For each direction, sorted squared centered projections are stored once.
At integer k, the lower and upper bounds are the kth and (n-k+1)th order
statistics, without interpolation. Both comparisons in depth evaluation are
inclusive. Ties can give depth greater than 1/2, so search covers 0 through n.
The zero level has an identity incumbent and needs no zeroth order statistic.

The default `domain = "pd"` maximizes a capped minimum-eigenvalue margin at
each level. A candidate requires an independently checked positive margin and
smallest eigenvalue exceeding the reported numerical threshold. The threshold
is a numerical acceptance criterion, not an imposed eigenvalue floor in the
optimization problem. Near-zero margins remain unresolved, except for direct
contradictions: a zero directional upper bound precludes every PD matrix.
Thus all-zero centered data has PD depth zero but PSD depth one.
`domain = "psd"` uses the PSD feasibility problem directly.

The conic adapter uses Clarabel's column-major upper-triangular svec, scaled
by sqrt(2) off the diagonal. The matrix cone uses matrix size d, not packed
size d(d+1)/2. Its coefficients C are generally dense even though the final
Clarabel A uses sparse storage. C and A construction/copy overhead is included
in the model-setup timing. The reference R kernels remain available internally
for checks; production defaults use RcppArmadillo matrix products, sorted
columns, binary-search counts, and symmetric eigendecomposition.

The center is subtracted first; centered observations are divided by their
maximum absolute entry (one for all-zero data). Scatter is multiplied by the
square of that divisor on return. Sigma, eigenvalues, margins, residuals, and
absolute tolerances are returned in original scatter units. Dimensionless
certificate diagnostics remain in scaled conic coordinates. Centering overflow,
unrepresentable squared scale, cap overflow, and projection overflow produce
actionable rescaling errors. Representable extreme magnitudes benefit from the
scalar normalization; no whitening or covariance calibration is applied.

`selected_k` is the last supported feasible search level; `requested_k` is the
last attempted level, which may be infeasible or unresolved. `attained_k` and
`depth_U` come from independent raw inclusive counts in the scaled coordinates,
without treating near-ties as exact ties. `raw_counts` and `tolerance_counts`
are separate diagnostics. Default absolute, relative and PD tolerances are
1e-8 in scaled scatter units. `tolerances` records their original-unit values.
Floating-point scalar normalization and matrix arithmetic may affect exact
boundary ties, so interpret `levels` and the bracket alongside attained depth.
The mapped output is checked again after conversion to original scatter units;
loss of positive definiteness from underflow requests rescaling. If mapping
changes a raw boundary count below the selected level, the completion status
becomes `"boundary_ambiguity"` rather than claiming attainment of that level.

Numerical infeasibility requires a direct order-statistic contradiction or a
supported `PrimalInfeasible` status with an independently validated dual witness.
The witness is normalized by its largest absolute entry. Stationarity is scaled
by one plus the maximum absolute column sum of A; separation by one plus the
absolute sum of b. Cone membership checks nonnegative entries and the PSD
slack eigenvalues. The validation tolerance is the absolute tolerance plus the
relative tolerance times max(1, max(Q)). These are numerical checks, not proofs
in exact arithmetic. Weak/inaccurate solver statuses do not by themselves
establish either feasibility or infeasibility.

Each unresolved level gets one retry with tighter solver tolerances. In PD
mode the retry uses a zero objective to seek an interior witness rather than
a boundary maximum-margin optimizer. Only validated positive witnesses are
accepted; this retry makes no zero-margin infeasibility claim. Equality rows
may receive an SVD linear correction; every condition is checked afterward.
No eigenvalue clipping or ridge is applied. Raw depth below the attempted
level is recorded as boundary ambiguity. If retry cannot resolve it, neither
search bound moves and `status` is `"partial_resolution"`. `search_bounds`
is the conditional numerical bracket for the optimum count, including when
only the k=0 incumbent is available.

The adapter checks the installed Clarabel exports and whether the actual solver
instance permits RHS updates. It retains a one-shot path (`control =
list(reuse = FALSE)`) and uses that path for retries. Solver version, complete
settings, call and update counts are returned. Persistent solver setup is timed
separately; one-shot setup and retry setup are included in feasibility-solves
time because the interface combines setup and solve. Timings use elapsed wall
time and may be zero at the platform clock's resolution. No speedup is assumed.

Increasing a nested direction set cannot increase the mathematical optimum.
For a fixed candidate, evaluate on the union of training U and independently
generated V for an empirical check: union depth cannot exceed training depth.
It remains an upper bound on full depth. Evaluating on V alone does not give
this monotonic relation. `examples/approx.R` demonstrates union validation.

Run `Rscript benchmarks/approx.R` from the repository root after installing the
package. It runs three seeds at n=500, d=10, L=2000, compares persistent and
one-shot paths, and performs a small PSD-domain 2D comparison using the full
angular oracle and existing exact optimizer. CSV results, full diagnostic RDS
objects, session/package versions, BLAS information and thread environment
variables go into `benchmarks/results/`. See `doc/approximate-validation.md`
for the recorded run and package-check evidence.

The adapter follows the official Clarabel R documentation for the
[one-shot interface](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel.html),
[persistent solver](https://oxfordcontrol.github.io/clarabel-r/reference/clarabel_solver.html),
and [RHS updates](https://oxfordcontrol.github.io/clarabel-r/reference/solver_update.html).
Installed version 0.11.3 was inspected and exercised during implementation.
