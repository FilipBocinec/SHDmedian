# Implementation and validation record

Run on 21 October 2026, Windows 11 x64, R 4.5.3, GCC 14.2.0, default C++17,
release optimization, Clarabel 0.11.3, Rcpp 1.1.2, RcppArmadillo 15.6.2-1,
Matrix 1.7-4, CVXR 1.9.3. Full versions, LAPACK/BLAS information, solver settings,
and thread environment variables are in `benchmarks/results/environment.txt`.
The thread environment variables were unset; no thread-count guarantee is made.

## Repository audit and milestones

The specification lives in `doc/shd_sdp_implementation_spec.md` (singular doc).
No AGENTS.md files were found in the repository or its ancestor directories.
The initial checkout contained substantial existing uncommitted package work;
that work was preserved. Existing exact algorithm sources and tests were not
edited. Existing public `SHDmedian()` performs median optimization, not just
depth evaluation. It uses CVXR, the angular arrangement Rcpp kernel, bivariate
Tukey/mean/supplied centering, linear direction cuts and a complete angular
separation oracle. Its public inputs and results remain unchanged.

1. Inspected the public API, exact search/master/oracle, dependencies, build
   configuration and regression tests. The initial baseline attempt could not
   load the test tooling because dependencies were missing. After installing
   workspace-local dependencies, the existing regression suite passed unchanged
   alongside the new implementation.
2. Added R reference preparation, column svec, unpacking, order-statistic bounds
   and inclusive counts. Exercised Clarabel directly and inspected its installed
   exports, controls, status descriptions and actual `x`, `z`, `s` result fields.
   Fixed-level tied examples worked through the R reference path before adding
   compiled kernels.
3. Added integer search, PD/PSD domains, independent candidate and dual-witness
   validation, bounded retry and conservative unresolved handling, public
   results and print method.
4. Added RcppArmadillo preparation/verification kernels, generated Rcpp exports
   with `Rcpp::compileAttributes()`, tested against the independent R kernels,
   and exercised both persistent RHS updates and the one-shot fallback.
5. Added generated help, mathematical/numerical documentation, executable
   examples, regression tests and reproducible benchmarks with saved diagnostics.

Changed implementation files: `R/SHDmedian_approx.R`, `R/approx-kernels.R`,
`R/approx-solver.R`, `src/approx_kernels.cpp`, generated Rcpp exports,
`src/Makevars`, `src/Makevars.win`, DESCRIPTION and generated NAMESPACE.
Package documentation metadata and README gained approximation information;
the exact algorithm implementation was preserved. Clarabel is a direct optional
dependency for the new wrapper. CVXR and mrfDepth retain their existing package
dependency roles. Build exclusions preserve an old `src/scatterHD.dll` on disk
while preventing that binary from entering source distributions.

## Tests and examples

The full source suite passed 558 assertions, with zero failures, errors,
warnings or skipped tests. It covers d=3/d=5 vectorization, R versus Rcpp
kernels, inclusive ties above depth one half, all-zero PSD versus PD behavior,
tiny PD margins, exhaustive small-level comparisons, numerical solver iteration
limits, mocked unresolved search outcomes, dual certificates, normalization,
translation, scaling, RNG restoration, nested directions and the full 2D oracle.
Existing exact regression tests passed unchanged.

`examples/approx.R` ran against the final installed release build: n=100, d=3,
L=100 produced a complete bracket [44,44], attained depth 0.44, and depth 0.44
on the union with 100 independent validation directions. Union validation is
empirical and does not establish full directional depth.

## Benchmark

Reproduce from the package root after installation:

```sh
Rscript benchmarks/approx.R
Rscript examples/approx.R
```

The final benchmark used n=500, d=10, L=2000 and default PD validation
tolerances. Presolve and chordal decomposition were disabled for this comparison
to permit RHS updates; all controls are recorded. Both solver paths used
identical fixed directions and produced identical per-level k/retry/status
sequences and final brackets, checked by the benchmark script.

| Seed | Attained k / depth | Optimum count bracket | Status | Persistent seconds | One-shot seconds | Calls per path |
| --- | --- | --- | --- | --- | --- | --- |
| 11 | 230 / 0.460 | [230,230] | complete | 3.28 | 2.72 | 13 |
| 29 | 226 / 0.452 | [226,233] | partial_resolution | 2.03 | 2.20 | 10 |
| 47 | 230 / 0.460 | [230,231] | partial_resolution | 3.44 | 2.68 | 13 |

Persistent runs performed 8, 6 and 8 RHS updates respectively. The unresolved
levels had solver primal-infeasibility statuses but their scaled dual separation
was too weak for the chosen independent certificate tolerance, including after
retry. The implementation therefore preserved its bracket instead of treating
the status alone as sufficient evidence. These are qualified numerical results.
Persistent reuse did not consistently reduce elapsed time on this machine.
The first model setup includes package/method initialization overhead.

`benchmarks/results/timings.csv` separates preprocessing, conic matrix setup,
persistent solver setup, feasibility solves (including one-shot setup/retries),
verification and total time. Timings can be zero at the clock's resolution.
`fits.rds` contains all level statuses, certificates, tolerances, U and matrices;
it is generated locally and excluded from version control. Runs on other machines
can differ with solver versions, floating-point arithmetic and system load.

The small PSD-domain 2D comparison (n=8, L=100) returned a complete finite-direction
depth of 0.375 in 0.02 seconds. Its full angular depth was 0.25; the existing
exact optimizer reached depth 0.375 in 2.47 seconds. This demonstrates omitted
direction error rather than establishing agreement of optimizer matrices.
Raw results are in `benchmarks/results/comparison-2d.csv`.

## Package tooling

Dependencies were installed into a workspace-local validation library because
the default R libraries lacked them. Rcpp attributes and Rd/NAMESPACE files were
generated using the package's established tooling. The release package compiled
and installed successfully with Rtools45 and default C++17.

Windows repeatedly denied the final staged-install directory rename, including
outside the sandbox. Non-staged installation succeeded in a fresh local library.
The package check therefore uses the explicit install option:

```sh
R CMD build --no-build-vignettes --no-manual .
R CMD check --no-build-vignettes --install-args=--no-staged-install SHDmedian_0.0.0.9000.tar.gz
```

The final check completed with **0 errors, 1 warning and 3 notes**. All 558
assertions passed in the installed-package check; examples, code/documentation
consistency, compiled-code checks and the PDF manual check passed. The remaining
findings are:

- The existing vignette source has no VignetteBuilder configuration and no
  installed built output (one warning and one note). Vignette rebuilding was
  not executed; no new vignette build system was introduced.
- The existing LICENSE file is not referenced by DESCRIPTION (one note).
- Non-staged installation was used to work around Windows directory-rename
  failures (one note). Standard staged installation could not be completed.

The final check log is
`.validation/final-check/SHDmedian.Rcheck/00check.log`, and installed test output
is `.validation/final-check/SHDmedian.Rcheck/tests/testthat.Rout`.
`git diff --check` passed. Validation was performed on Windows with R 4.5.3;
other platforms and R versions were not tested.
