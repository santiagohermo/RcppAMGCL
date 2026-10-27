# Validation

Latest recorded checks: 27 October 2026, source
[`9de2012`](https://github.com/santiagohermo/RcppAMGCL/commit/9de20125f1d34e4988499888868b484459337288).

## Local check

`R CMD check --as-cran` on Windows 11, R 4.5.3 and Rtools GCC 14.3.0:
**0 errors, 0 warnings, 1 NOTE: New submission**. Tests, examples and both
PDF and HTML manuals passed, using an OpenMP build with a two-thread cap.
See the [check log](windows-standalone-as-cran.log).

Tested source archive SHA256:
`c3fe5634ef570e63e38bf216de8dab3af87eed7cab4d3ad54b63dc4846d1e829`.

## GitHub Actions

All six configurations in [this run](https://github.com/santiagohermo/RcppAMGCL/actions/runs/40488244633)
passed with `Status: OK` using R 4.6.1.

| Platform | Compiler | OpenMP on | OpenMP off |
| --- | --- | --- | --- |
| Ubuntu 24.04 | GCC 13.3.0 | Passed | Passed |
| Windows | Rtools GCC 14.3.0 | Passed | Passed |
| macOS | Apple clang 21.0.0 | Passed | Passed |

Each job installs the package, verifies the selected backend, runs R and
standalone C++14 tests, and runs `R CMD check --as-cran --no-manual` with
remote incoming checks disabled. The local check includes the manual and
incoming checks. Linux jobs needed retries for slow/stalled R installation;
the final combined result is in attempt 3, retaining completed jobs.

Tests cover both solvers, symmetric and nonsymmetric matrices, reusable
objects, input validation, ownership, extreme RHS scales, independent
residuals, and thread limits/restoration. Tests and examples use at most
two threads. The separate [performance comparison](../benchmarks/README.md)
completed 900 measured solves with no accuracy failures on Windows.
