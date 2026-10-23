# Exact bivariate scatter median

Use `SHDmedian(X, method = "exact", center = mu)`. C++ is compiled during
package installation; no runtime compilation is needed. Both methods share
centering. Configure geometry and tolerances through `control`.

## Numerical results and status

The formulation below is exact in real arithmetic. This implementation uses
floating-point geometry and a numerical conic solver, so its certificates are
**tolerance-based, not formal exact certificates**.

- `fit$best_k`, `fit$depth`: the accepted depth level, using the stated tolerances.
- `fit$search_details$strict_k`, `fit$search_details$strict_depth`: an independent angular root-sweep check
  with zero projection tolerance, still computed in floating-point arithmetic.
- `fit$search_bounds`: numerical lower and upper bounds on the deepest level.
- `fit$certified = TRUE`: that numerical interval has collapsed. Read this
  together with `fit$search_details$certificate_type`, not as a symbolic certificate.
- `fit$status = "unresolved"`: a solver/boundary decision remains open. The
  returned matrix is a candidate at the lower endpoint; it is not claimed to be
  a proven median. The reason is recorded in `fit$search_details$levels`.

Default tolerances are `solver_tol = 1e-9`, `validation_tol = 1e-7`,
`pd_tol = 1e-9`, and `angle_tol = 1e-12` radians. Projection validation and the
positive-definiteness threshold are relative to the common squared-data scale.
Nearby angular events are grouped using `angle_tol`; extremely close events
and ill-conditioned data therefore require care. `angle_tol = 0` disables this
grouping tolerance, but does not provide exact arithmetic.

The auxiliary optimization maximizes δ subject to Σ − δI being PSD and δ ≥ 0.
A validated positive definite candidate establishes numerical feasibility.
CLARABEL's `PrimalInfeasible` status establishes numerical infeasibility of the
closed conic problem. A zero or very small δ is **not** interpreted as a proof
that the region contains no strictly positive definite matrix.

When roundoff puts the final candidate just outside a depth boundary, an
optional extra SOCP seeks an interior point by maximizing a common lower/upper
slack while retaining a fixed positive eigenvalue margin. This often makes the
zero-tolerance depth agree with the accepted level. Lower-dimensional median
regions may not admit such a strictly interior point.

The chosen representative is not the maximum-determinant member of the median
region. The eigenvalue-margin objective and numerical selection do not define
an affine-equivariant selection rule. The exact depth region itself is the
original one.

## Mathematical formulation

Write Y_i = X_i − μ, S_i = Y_i Y_iᵀ, u(θ) = (cos θ, sin θ),
and r_i(θ) = (u(θ)ᵀY_i)², with projective angles 0 ≤ θ < π.
The level-k region is

\[
\Sigma\succ0,\qquad
r_{(k)}(\theta)\le u(\theta)^\top\Sigma u(\theta)
\le r_{(n-k+1)}(\theta)\quad\text{for every }\theta.
\]

For each observation i separately, the C++ sweep sorts only the angles
perpendicular to Y_i − Y_j and Y_i + Y_j. It updates the number of projections
strictly smaller than r_i, handling identical outer products and simultaneous
events. It does not construct global arrangement cells or their permutations.

**Lower constraints.** Let A_i be the union of positive-length closed angular
arcs on which at most k−1 observations have projection square strictly below
r_i. Requiring uᵀ(Σ−S_i)u ≥ 0 on every A_i is equivalent to the lower system.
Isolated active directions can be omitted from this per-point description:
the complete system recovers their order-statistic constraint by continuity
from adjacent nonempty open intervals.

Put M_i = Σ−S_i and

\[
\alpha_i=(M_{i,11}+M_{i,22})/2,\qquad
\beta_i=((M_{i,11}-M_{i,22})/2,M_{i,12})^\top.
\]

Then uᵀM_i u = α_i + β_iᵀw with w = (cos 2θ, sin 2θ).
The convex hull of the doubled-angle active arcs is the unit disk clipped
by one halfspace for each inactive angular gap (a,b), unwrapped so b>a:

\[
G_i=\{w:\|w\|_2\le1,\ N_iw\le d_i\},\qquad
N_{ig}=(\cos(a+b),\sin(a+b)),\quad d_{ig}=\cos(b-a).
\]

Conic duality gives the exact finite condition

\[
\boxed{\quad
\exists\lambda_i\ge0:\quad
\|\beta_i+N_i^\top\lambda_i\|_2
\le\alpha_i-d_i^\top\lambda_i.
\quad}
\]

Thus there is at most **one three-dimensional SOC per nonzero observation**,
plus nonnegative multipliers for its inactive gaps. If all directions are
active, there are no multipliers and this is simply Σ−S_i PSD. There is one
additional three-dimensional SOC for the positive-definiteness margin of Σ.

**Upper constraints.** The same local sweep identifies arcs where at most
k−1 projection squares are strictly larger than r_i. For Σ PSD, on an arc
where uᵀY_i has constant sign, it suffices to check its two endpoint rays:

\[
u_a^\top\Sigma u_a\le (u_a^\top Y_i)^2,\qquad
u_b^\top\Sigma u_b\le (u_b^\top Y_i)^2.
\]

Indeed, every direction in the arc is a positive combination of its endpoint
rays; the triangle inequality for sqrt(uᵀΣu) and linearity of |uᵀY_i| on that
cone extend the endpoint bounds to the entire arc. These bounds are linear
in the three entries of Σ. If the arc contains a direction perpendicular
to Y_i, its upper bound is zero there and no SPD matrix can satisfy it.
Consequently, an external computation of a symmetrized location-depth region
is not needed for this implementation.

Binary search uses nested depth regions and the upper bound floor(n/2).
It starts with the independently evaluated depth of crossprod(Y)/n,
rather than assuming that a stored matrix already exists at a theoretical
lower bound. This also accommodates explicitly supplied centers.

## Cost and practical limits

For each i there are at most 2(n−1) angular events. Geometric preprocessing
therefore costs O(n² log n). Cached event angles and ranks use O(n²) memory;
for a generic n=1000 sample this cache is approximately 22.9 MiB. Each
subsequent level scans these cached ranks in O(n²), then builds and solves
its sparse SOCP. Binary search makes O(log n) level decisions.

The number of lower SOC blocks is at most n, but the number of gap multipliers
and upper endpoint constraints can still be O(n²). There is no claim of a
worst-case linear-size problem or a universal speed advantage over cutting
planes. With `cache_geometry = FALSE`, local geometry is recomputed each level:
this saves the O(n²) event cache but not the memory needed by the conic model.
The default `"auto"` disables caching above a 256 MiB estimated event budget;
this is not a bound on total process memory.


