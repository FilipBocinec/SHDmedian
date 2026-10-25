# SOCP benchmarks

Run `benchmark_socp.R` in this directory after installing SHDmedian.
It measures the same fixed-center SOCP implementation used by
`SHDmedian(X, method = "exact", center = c(0, 0))`, including geometry,
level searches, validation and refinement. Internal helpers expose component
timings for benchmarking; they are not public package functions.

`compare_shd_methods.R` defines a comparison helper accepting a supplied
cutting-plane implementation and a common center.

Running the benchmark creates `benchmark_results.csv` and
`benchmark_session.txt` in this directory. These generated results are ignored
by Git.
