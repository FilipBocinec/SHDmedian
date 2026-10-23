# SOCP benchmarks

Run `benchmark_socp.R` in this directory after installing SHDmedian.
It measures the same fixed-center SOCP implementation used by
`SHDmedian(X, method = "exact", center = c(0, 0))`, including geometry,
level searches, validation and refinement. Internal helpers expose component
timings for benchmarking; they are not public package functions.

`compare_shd_methods.R` defines a comparison helper accepting a supplied
cutting-plane implementation and a common center.

The CSV files and `benchmark_session.txt` are historical results supplied
with the original standalone implementation, not timings from this integration.
Running the benchmark overwrites the timing CSV and session file.
