# Benchmarks

The 2D and 3D Poisson matrices come from finite-difference grids with five
and seven nonzeros per interior row. They are deterministic test problems,
not samples from a Poisson distribution. These sparse, symmetric positive
definite systems suit AMG; results do not predict performance on every
sparse matrix.

## Comparison with Matrix

Measured on 27 October 2026: Intel Core Ultra 7 265H, 16 logical processors,
about 64 GiB RAM, Windows 11, R 4.5.3, Matrix 1.7-4, and Rtools GCC 14.3.
The run used Windows R's supplied Rblas.dll. BLAS thread environment
variables requested one thread; vendor thread counts were not queried.

Median wall-clock seconds for the largest cases: **100,489 unknowns in 2D**
and **103,823 in 3D**. Each configuration has three measurements.

| Matrix | Method | Threads | Setup (s) | First solve (s) | Total, 1 RHS (s) | Total, 10 RHS (s) | Peak RSS (MiB) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2D | AMGCL CG | 1 | 0.0907 | 0.1466 | 0.2373 | 1.3974 | 260.0 |
| 2D | AMGCL CG | 2 | 0.0641 | 0.1055 | 0.1781 | 1.0093 | 261.4 |
| 2D | AMGCL CG | 4 | 0.0582 | 0.0777 | 0.1359 | 0.9511 | 262.1 |
| 2D | AMGCL CG | 8 | 0.0452 | 0.0805 | 0.1264 | 0.6844 | 261.8 |
| 2D | Matrix Cholesky | 1 | 0.3076 | 0.0133 | 0.3209 | 0.4586 | 313.3 |
| 3D | AMGCL CG | 1 | 0.1767 | 0.1868 | 0.3478 | 2.1090 | 290.1 |
| 3D | AMGCL CG | 2 | 0.0902 | 0.0868 | 0.1722 | 1.0408 | 289.5 |
| 3D | AMGCL CG | 4 | 0.0873 | 0.0780 | 0.1654 | 1.0845 | 291.6 |
| 3D | AMGCL CG | 8 | 0.0530 | 0.0337 | 0.0871 | 0.4618 | 292.1 |
| 3D | Matrix Cholesky | 1 | 19.7793 | 0.0778 | 19.8571 | 20.5156 | 802.8 |

AMGCL has lower setup costs in these cases. Cholesky's cheap repeated solves
give it the lower ten-RHS total in 2D; AMGCL wins both totals in 3D. BLAS
implementation, matrix structure, size, accuracy and thread count can change
this comparison.

### Method

- The same symmetric `dsCMatrix` and ten distinct known solutions generate
  the inputs for both methods. Tested sizes are 10,000, 40,000 and 100,489
  in 2D, and 10,648, 42,875 and 103,823 in 3D.
- AMGCL uses CG with `tol = 1e-8` and thread caps 1, 2, 4 and 8.
  Matrix uses `Cholesky(A, perm = TRUE, LDL = FALSE, super = NA)`, letting
  CHOLMOD choose its factor type. Every measured case chose `dCHMsuper`.
- Every setup receives a fresh matrix object with an empty factor cache.
  Both methods retain their solver/factor for the ten solves. One setup and
  solve warm up each configuration before measurement.
- Setup and solve are timed separately with a monotonic, high-resolution
  clock. Matrix/RHS construction, explicit garbage collection and external
  validation are outside the timers. AMGCL's own residual check remains
  inside its public solve time.
- `total_seconds` is setup plus the first solve; `batch_total_seconds` is
  setup plus all ten solves. Totals are calculated within each repetition
  **before** taking medians, so displayed column medians need not add up.
- A fresh R process runs each configuration, with a 300-second timeout.
  Peak RSS is its lifetime high-water mark, including R, input generation,
  warmup, all repetitions and validation; it is not incremental factor memory.

All **900 measured solves passed**: relative residual at most `1e-8` and
relative solution error at most `1e-5`. Maximum observed residuals were
`9.94e-9` for AMGCL and `2.28e-15` for Cholesky; maximum solution errors
were `5.41e-7` and `4.15e-14`. Cholesky thus delivered much tighter accuracy.

[Raw measurements](windows-comparison.csv), [session and hashes](windows-comparison.csv.session.txt),
and [run report](windows-comparison.csv.report.txt). The harness is pinned at
`91588a0`; the installed solver code is from `771186f`, unchanged by the
benchmark and documentation updates.

### Reproduce

Install an OpenMP build of RcppAMGCL. Copy `tools/benchmark-memory.cpp` into
a temporary directory and run `R CMD SHLIB benchmark-memory.cpp` there.
From the repository root, supply the resulting library's path:

```sh
Rscript tools/compare-solvers.R benchmarks/results-local.csv 10000,40000,100000 3 10 /path/to/benchmark-memory.dll
```

Use the resulting `.so` on Linux/macOS. The helper was executed on Windows
for this dataset. Per-case logs and metadata are saved under
`.local/compare-solvers/`. The harness respects a lower inherited
`OMP_THREAD_LIMIT`; otherwise it caps the run at eight threads.

## Larger OpenMP scaling run

These earlier measurements use 1,000,000 unknowns in each dimension
(4,996,000 nonzeros in 2D; 6,940,000 in 3D), on the same host and solver
source `771186f`. Three warmed repetitions, CG, `tol = 1e-8`:

| Matrix | Threads | Setup (s) | Solve (s) | Total (s) |
| --- | ---: | ---: | ---: | ---: |
| 2D | 1 | 0.50 | 0.68 | 1.18 |
| 2D | 2 | 0.31 | 0.41 | 0.72 |
| 2D | 4 | 0.22 | 0.27 | 0.49 |
| 2D | 8 | 0.17 | 0.20 | 0.36 |
| 3D | 1 | 1.11 | 0.96 | 2.07 |
| 3D | 2 | 0.69 | 0.61 | 1.28 |
| 3D | 4 | 0.46 | 0.42 | 0.89 |
| 3D | 8 | 0.32 | 0.27 | 0.59 |

All 24 solves converged. [Raw data](windows-openmp-review.csv) and
[session](windows-openmp-review.csv.session.txt). Reproduce with:

```sh
Rscript tools/benchmark.R benchmarks/scaling-local.csv 3 1000000
```

Historical [smaller](windows-openmp.csv) and [million-row](windows-openmp-large.csv)
runs remain available, along with [native memory/residual diagnostics](native-memory-residual.md).
Total columns in existing CSVs were derived from their original measurements;
those measurements were not rerun.

Benchmarks are excluded from the source package. CRAN tests and examples
use at most two threads; these performance experiments run separately.
