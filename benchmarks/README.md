# Setup and solve benchmarks

Measured on 27 October 2026 using the checked OpenMP package from source
commit `ec8bd801c6583a26601aee9140547be80111d680`. Host: Intel Core Ultra 7
265H, Windows 11, R 4.5.3, Rtools GCC 14.3.0. `amgcl_info()` reported
16 available threads; `OMP_THREAD_LIMIT` and `OMP_NUM_THREADS` were unset.

The main run used two sparse SPD Poisson matrices, each with 1,000,000
unknowns: 2D with 4,996,000 nonzeros and 3D with 6,940,000 nonzeros.
CG used the same manufactured RHS and `tol = 1e-8` at every thread count.
Each configuration was warmed up once, followed by three measurements.
Matrix/RHS construction and garbage collection were outside timed sections.

Setup includes conversion from R's CSC matrix to owned CSR storage and
construction of the AMG hierarchy. Solve reuses that hierarchy and includes
RHS conversion, iteration, the independent residual calculation, and
returning the result to R. Times below are median wall-clock seconds.

| Matrix | Threads | Setup (s) | Solve (s) |
| --- | ---: | ---: | ---: |
| 2D Poisson | 1 | 0.50 | 0.75 |
| 2D Poisson | 2 | 0.31 | 0.50 |
| 2D Poisson | 4 | 0.23 | 0.36 |
| 2D Poisson | 8 | 0.18 | 0.33 |
| 3D Poisson | 1 | 1.05 | 0.92 |
| 3D Poisson | 2 | 0.64 | 0.60 |
| 3D Poisson | 4 | 0.42 | 0.44 |
| 3D Poisson | 8 | 0.30 | 0.36 |

All 24 solves converged. CG took 12 iterations for 2D and 13 for 3D;
maximum independently measured relative residuals were approximately
`4.74e-9` and `6.23e-9`, respectively. Effective requested thread caps
matched 1, 2, 4, and 8. Timings vary with hardware, workload, and system load.
These measurements cover CG; BiCGStab correctness is covered by the test suite.

Raw data and environment records:

- [Million-unknown measurements](windows-openmp-large.csv)
- [Million-unknown session metadata](windows-openmp-large.csv.session.txt)
- [160,000/166,375-unknown measurements](windows-openmp.csv)
- [Smaller-workload session metadata](windows-openmp.csv.session.txt)

The smaller run shows more overhead and timing variability, including
slower 2D solves above two threads. It is retained alongside the larger run.

To reproduce with an installed OpenMP build:

```sh
Rscript tools/benchmark.R benchmarks/results-local.csv 3 1000000
```

The script defaults to a minimum of 160,000 unknowns when the third argument
is omitted. Benchmark scripts and records are excluded from the CRAN source
package and are never run as CRAN tests or examples.
