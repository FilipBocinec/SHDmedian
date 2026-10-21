# Implementation task: approximate scatter halfspace median using finite directions

This is a specification for the coding agent working in the user's existing R package.
The package already supports exact calculations in dimension two. Inspect the actual
repository before choosing names, defaults, or integration points.

## 1. Goal and scope

Add an approximate scatter halfspace median for arbitrary dimension. For a fixed set
of directions, maximize empirical directional depth by searching nested convex
feasibility problems. Use Clarabel for the semidefinite problems, R for orchestration,
and Rcpp/RcppArmadillo for preprocessing and validation.

Keep the existing exact two-dimensional API and behavior. Add an explicit approximate
method or a new function consistent with the package. Do not assume the existing code
computes both depth evaluation and median optimization: determine this in the audit.

The deliverable is a working implementation, documentation, meaningful tests,
an executable example, and a benchmark script. Follow the milestones below through
implementation and validation; an audit alone is not the requested outcome.

The result is one maximizing matrix, subject to numerical qualifications. It is not
automatically the barycenter of the median set. Do not describe the selected matrix as
the barycenter unless that separate computation is actually implemented.

## 2. Mathematical contract

### Data, center, and directions

Let X be a finite real n-by-d data matrix, with observations in rows, and let mu be a
fixed d-vector. Define Y_i = X_i - mu. The center is fixed throughout optimization.

The manuscript uses a location halfspace median as its center. Accept an explicit
center, and reuse the package's established centering convention where applicable.
If a suitable higher-dimensional centering method is absent, require the center
argument. Do not silently replace it by the sample mean or componentwise median.
Joint optimization of center and scatter is outside this task.

Use a d-by-L matrix U whose columns are nonzero directions normalized to unit length.
Validate dimensions and finite values. Reject zero directions and nonfinite inputs.
Document whether zero-length dimensions and empty samples are supported or rejected.

For automatically generated U, a practical initial rule is to include the d coordinate
directions and generate the remaining L-d normalized Gaussian directions, with L >= d.
L denotes the total number of directions. Explicit user directions must be honored;
do not silently append directions to a user-supplied U. Respect the package's R RNG
conventions, document the seed behavior, and return the actual U.

Keep U fixed for every depth level in one optimization. Antipodal directions impose
identical constraints. Do not resample directions during binary search.

### Directional depth and order statistics

For each i and j define

$$
q_{ij}=(u_j^\top Y_i)^2,\qquad
q_{(1)j}\le\cdots\le q_{(n)j}.
$$

For a positive semidefinite Sigma, define

$$
D_U(\Sigma)=\frac1n\min_{1\le j\le L}
\min\left\{
\#\{i:q_{ij}\le u_j^\top\Sigma u_j\},
\#\{i:q_{ij}\ge u_j^\top\Sigma u_j\}
\right\}.
$$

Both comparisons are inclusive. An observation equal to the threshold contributes
to both counts. For integer 1 <= k <= n,

$$
D_U(\Sigma)\ge k/n
\quad\Longleftrightarrow\quad
q_{(k)j}\le u_j^\top\Sigma u_j\le q_{(n-k+1)j}
\quad\text{for every }j.
$$

Use these order statistics directly. Do not use interpolated sample quantiles.
If a lower bound exceeds its upper bound, the level is infeasible.
Handle k=0 separately; there is no zeroth order statistic.

Ties can make finite-direction depth exceed 1/2. Search the safe integer range
0,...,n; do not impose a universal floor(n/2) upper bound.

### Symmetric vectorization

Let p=d(d+1)/2 and use x=svec_col(Sigma) internally, with upper-triangular entries
stacked column by column and off-diagonal entries multiplied by sqrt(2).
For d=3 the order is

$$
(\sigma_{11},\sqrt2\sigma_{12},\sigma_{22},
 \sqrt2\sigma_{13},\sqrt2\sigma_{23},\sigma_{33}).
$$

Set C[j, ] = svec_col(u_j u_j^T)^T. Then (C x)_j = u_j^T Sigma u_j.
Implement tested pack and unpack helpers.

The manuscript's svec uses the upper triangle ROW BY ROW. If the existing package
uses that convention, retain its external contract and use an explicit permutation
at the solver boundary. The two conventions coincide for d=2 but differ for d>=3.

### PSD feasibility problem

For a fixed k, set ell_j=q_(k)j and h_j=q_(n-k+1)j and solve

$$
\ell\le Cx\le h,\qquad \Sigma\succeq0.
$$

This is an SDP. Directional inequalities alone do not guarantee that Sigma is PSD.
Use the matrix cone directly; a Cholesky parametrization would destroy linearity
in the decision variables.

Clarabel uses A x + s = b with s in a product cone. The PSD feasibility form is

$$
A=\begin{bmatrix}C\\-C\\-I_p\end{bmatrix},\quad
b=\begin{bmatrix}h\\-\ell\\0_p\end{bmatrix},\quad
q_{\rm objective}=0_p,\quad
K=\mathbb R_+^{2L}\times\mathcal S_+^d.
$$

The last slack block is svec_col(Sigma). Use a zero quadratic objective.
Here the PSD cone size is d, although its packed slack has p entries.

### Positive definiteness: preserve the intended estimator

The manuscript defines a scatter median over POSITIVE DEFINITE matrices.
Implement the PSD problem as the core and as an explicitly selectable domain,
but expose a PD-aware mode for the intended median. Do not silently return a
singular PSD optimizer as a PD median.

For a fixed k, PD feasibility can be tested through

$$
\max_{\Sigma,t}\ t
\quad\text{subject to}\quad
\ell\le C\,\operatorname{svec}_{col}(\Sigma)\le h,\quad
\Sigma-tI_d\succeq0,\quad 0\le t\le t_{\rm cap},
$$

where t_cap is any finite positive number in scatter units. Choose and record a
reasonable data-dependent scale, with a positive fallback for all-zero data.
The upper cap prevents an unbounded objective for incomplete direction sets.
In exact arithmetic, a PD feasible matrix exists if and only if the optimum is
positive. The cap does not remove any possible PD-feasibility conclusion.

For w=(x,t), put v=svec_col(I_d). An explicit conic encoding is

$$
A_{\rm PD}=
\begin{bmatrix}
C&0\\-C&0\\0^\top&1\\0^\top&-1\\-I_p&v
\end{bmatrix},\quad
b_{\rm PD}=
\begin{bmatrix}h\\-\ell\\t_{\rm cap}\\0\\0_p\end{bmatrix},
\quad q_{\rm PD}=(0_p,-1).
$$

The nonnegative cone has 2L+2 entries and the PSD cone has matrix size d.
The final slack is svec_col(Sigma-t I_d).

A strictly positive, validated t is a PD witness. An optimum numerically near zero
is not, by itself, proof that PD feasibility is impossible. Return an unresolved
status when the numerical evidence cannot distinguish zero from a small positive
margin. A PSD infeasibility certificate also rules out PD feasibility.

An optional fixed eigenvalue floor Sigma >= epsilon I can be implemented later
as an explicitly named restricted-domain mode. Positive epsilon changes the
optimization problem; it must not silently stand in for the strict PD definition.
Adding a diagonal ridge or clipping eigenvalues after solving can violate the
directional upper bounds. Revalidate every condition after any modification.

## 3. Numerical decisions and search

Use three outcomes for a level: feasible, infeasible, unresolved. Report the
numerical tolerances used. Ordinary floating-point solver output is not an
exact-arithmetic proof.

- Feasible requires a candidate that passes an independent reconstruction,
  directional-constraint check and domain check. For PD, inspect the actual
  smallest eigenvalue and the margin; merely obtaining a solver success code
  is insufficient.
- Infeasible requires a supported infeasibility result with a valid numerical
  certificate, or a direct mathematical contradiction such as ell_j > h_j.
- Time limits, insufficient progress, numerical errors and weak/inaccurate
  termination statuses are unresolved unless independently justified.

For the PSD conic feasibility problem, a dual infeasibility witness satisfies
A^T z = 0, b^T z < 0, and z in the dual cone. Check scaled residuals and cone
membership when relying on this witness. The cones here are self-dual.
Read status meanings and result fields from the installed R interface rather
than guessing status integers or copying Python field names.

Start binary search with low=0 and high=n+1, where high is an impossible sentinel.
Store a valid incumbent, such as a positive multiple of the identity at k=0.
At each iteration, try k=floor((low+high)/2):

- If feasible, store the candidate and set low=k.
- If infeasible, set high=k.
- If unresolved, perform a bounded retry with justified scaling or tolerances.
  If it remains unresolved, preserve the bracket and report partial resolution.
  Do not silently move either bound.

Stop when high=low+1 or when numerical uncertainty prevents further progress.
The bracket for the optimum is [low, high-1], conditional on the documented
numerical validations. In PD mode, unresolved zero-margin levels may leave a gap.

Independently calculate the attained D_U of the returned matrix with inclusive
counts. Report raw floating-point counts separately from tolerance-expanded
diagnostics. Tolerances must not silently turn near-ties into exact ties.
If the attained raw depth is below the attempted level, try a bounded numerical
refinement and otherwise report boundary ambiguity rather than claiming that
the returned matrix attains that level.

Use scale-aware absolute and relative residual tolerances. Handle all-zero
data and extreme finite magnitudes explicitly; check for overflow in squared
projections. If scalar rescaling is introduced, record its transformation and
map Sigma, eigenvalue margins and tolerances back consistently.

## 4. Software architecture

Adapt names and object classes to the package; the following division of work
is the intended design:

1. Public R wrapper: validate input, resolve the fixed center, generate or accept
   U, choose domain, orchestrate search, format the result.
2. Preparation kernel: form Q=(Y U)^2 elementwise, sort each column once, construct
   C, and return dimensions and reusable data.
3. Solver adapter: construct the conic matrix once for each selected domain and
   change the RHS as k changes. Keep solver conventions inside this adapter.
4. Search routine: implement the monotone integer search and unresolved-state
   handling independently of preprocessing.
5. Verification kernel: recompute directional quadratic forms, counts, residuals
   and eigenvalue information from the returned matrix.

Begin with small transparent R reference implementations for the preparation
and count calculations. They are independent correctness checks for the Rcpp
versions. Move the projection, per-column sorting and repeated verification
work into Rcpp/RcppArmadillo after the reference implementation works.

Use matrix multiplication/BLAS for Y U. Avoid allocating an n-by-d-by-L object
or lifting all observations into p dimensions unnecessarily. At first, retain
one sorted n-by-L projection matrix and reuse it. Counts can use lower_bound
and upper_bound on sorted columns; preserve inclusive tie semantics.

For large dense C, do not claim that the coefficients become sparse merely
because a sparse storage class is used. Convert the final solver data to its
required format and account for construction/copy overhead in the benchmark.

### R solver integration

The documented one-shot interface is:

~~~r
clarabel::clarabel(A, b, q, P = NULL, cones, control = list())
~~~

Use integer cone sizes, for example a nonnegative block of size 2L followed by
a PSD block of size d. In R, write as.integer(2 * L), not 2L: the latter is the
integer literal two.

Current documentation also exposes persistent solver construction and RHS updates.
Check the installed package exports and help before choosing this path:

~~~r
solver <- clarabel::clarabel_solver(
  A = A, b = b, q = objective, P = NULL,
  cones = cones, control = control
)
clarabel::solver_is_update_allowed(solver)
clarabel::solver_update(solver, b = b_next)
result <- clarabel::solver_solve(solver)
~~~

Use this only when those functions exist and updates are allowed for the actual
instance. Structural transformations can prevent updates. Keep a correct one-shot
fallback and compare the two on the same levels. Reusing an object does not justify
promising a particular speedup. Record the installed version and settings.

Treat Clarabel as an optional dependency if consistent with package policy:
existing exact routines should remain usable when the optional solver is absent.
Provide an actionable error from the approximate routine in that case.
Use established Rcpp export/documentation generation and package tooling.

Suggested public interface, to adapt after inspection:

~~~r
shd_median_approx(
  X, center, U = NULL, n_directions = 2000L,
  seed = NULL, domain = c("pd", "psd"), control = list()
)
~~~

Return at least Sigma, center, U, domain, requested/selected k, attained k,
depth_U, search lower/upper bounds, completion status, smallest eigenvalue,
maximum constraint violation, per-level statuses, tolerances and elapsed times.
Use explicit names to distinguish attained depth from the solver's target.
Do not label a floating-point result simply "exact" or "certified".

## 5. Tests that protect the mathematics

Use the project's existing test framework and add focused tests for:

1. Packing/unpacking and the quadratic-form identity in d=3 and d=5, including
   nonzero off-diagonal entries. Dimension two alone misses the ordering bug.
2. Small hand-computed order-statistic bounds, including repeated values and
   levels above floor(n/2). For example, in d=1, mu=0, X=(-1,1), Sigma=1 has
   depth 1, not 1/2.
3. PSD versus PD behavior. For all-zero centered data, Sigma=0 has PSD depth 1,
   whereas every PD matrix has finite-direction depth 0 for any nonempty U.
   The implementation must not report a singular answer as a PD median.
4. Comparison of small R reference kernels with Rcpp outputs.
5. Exhaustive enumeration of k=0,...,n for small examples versus binary search.
   Compare outcomes and brackets, not identical optimizer matrices.
6. Simulated numerical failure/unresolved statuses: they must not become
   infeasibility or advance an unsupported search bound.
7. Repeated directions, sign flips and normalization of nonunit directions.
8. Translation when X and center shift together, and scalar rescaling with the
   corresponding Sigma transformation.
9. Increasing a nested direction set cannot increase the mathematical optimum
   when center and matrix domain are unchanged, allowing for reported numerical
   uncertainty.
10. Existing exact 2D routines, insofar as actually available. Evaluate the full
    exact depth of the approximate matrix; if exact optimization is available,
    compare attained depths with its optimum. Nonunique medians need not produce
    identical matrices. Preserve existing regression tests.

For every candidate in the same domain and with the same center,
D_full(Sigma) <= D_U(Sigma). Maximizing over U does not remove approximation
error from omitted directions.

For additional empirical validation, generate independent directions V and
evaluate on U union V. This depth can only decrease relative to the original U
for a fixed candidate. It remains an upper bound on full depth, not a certificate
that all directions have been checked. Comparing on V alone has no such
monotonic relationship with the value on U.

## 6. Milestones and completion evidence

Execute the following milestones in order, using small reviewable changes.
Follow existing AGENTS.md instructions and preserve unrelated working changes.

1. Audit package interfaces, exact 2D capabilities, dependencies, tests and build
   commands. Run the baseline checks that the environment supports.
2. Add and validate reference preparation/vectorization/depth routines and a
   working solver for one fixed k, including the intended domain.
3. Add the search, numerical status handling and public result object.
4. Add the Rcpp kernels and persistent-solver reuse where supported. Retain
   reference checks and the one-shot fallback.
5. Add user documentation, an executable example, regression checks and benchmark.

Benchmark at n=500, d=10, L=2000 with several fixed seeds, plus a small 2D
comparison. Record preprocessing, solver setup, all feasibility solves,
verification and total elapsed time separately; record solver-call counts,
statuses, package versions and relevant BLAS/thread settings.
Avoid making a hard runtime promise from one machine or from a Python benchmark.

Run the existing package checks and targeted new tests. Where practical,
run R CMD check. Clearly report anything not executed and why.
Produce a concise final report of changed files, example usage, tests,
benchmark results and remaining numerical limitations.

Do not introduce calibration factors for covariance consistency, trace or
determinant constraints, whitening, adaptive directions or a new centering
estimator unless separately specified. Those alter the target or add another
approximation and are outside this initial implementation.

## 7. References for the coding agent

The mathematical finite-direction formulation above is the task specification.
Inspect the installed implementation as well as current official documentation:

- R solver interface: https://oxfordcontrol.github.io/clarabel-r/reference/clarabel.html
- Cone vectorization and examples: https://oxfordcontrol.github.io/clarabel-r/articles/clarabel.html
- Persistent solver: https://oxfordcontrol.github.io/clarabel-r/reference/clarabel_solver.html
- RHS updates: https://oxfordcontrol.github.io/clarabel-r/reference/solver_update.html
- Update availability: https://oxfordcontrol.github.io/clarabel-r/reference/solver_is_update_allowed.html

Prepared on 21 October 2026. Documentation availability is not a guarantee that
the user's installed package version exposes every function listed above.
