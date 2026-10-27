# Validation

## Release preparation: 28 October 2026

The final attribution-update source archive was built and checked on
Windows 11 with R-devel 4.7.0 (2026-10-26 r90766), Rtools GCC 14.3.0,
Rcpp 1.1.2, Matrix 1.7-6 and V8 8.2.0. Both builds used R's default
C++20 standard and a two-thread cap. The complete
`R CMD check --as-cran` checks included tests, examples, PDF manuals,
and HTML manuals with math rendering enabled.

| Build | Remote incoming checks | Errors | Warnings | Notes |
| --- | --- | ---: | ---: | --- |
| OpenMP on | Enabled | 0 | 0 | 1: New submission |
| OpenMP off | Disabled | 0 | 0 | 0 |

Check logs: [OpenMP on](windows-r-devel-openmp-as-cran.log) and
[OpenMP off](windows-r-devel-serial-as-cran.log). The absolute repository
path in these logs is replaced with `<workspace>`.

Both checks used the same final `RcppAMGCL_0.1.0.tar.gz`, SHA256:
`97c5a9f05b6a93db4a2bca84fc7e63a5b99af2b4b29d1105a4f20bcb6df2e2bd`.
The archive contains the copyright inventory and third-party notices.
The inventory was checked against all 15 bundled headers with additional
copyright statements; the reproduced license text was checked against
the five headers containing the Kratos, LAPACK, SciPy and PETSc notices.

The workflow now configures 12 jobs: R-release on Linux, Windows and
macOS; R-devel on Linux and Windows; and R 4.3.3 on Linux. Each runs with
OpenMP on and off. R-devel jobs install V8 to enable HTML math-rendering
checks. The newly added jobs have not yet run on GitHub; the results above
are local R-devel checks. The earlier six-job CI evidence is retained below.

## Earlier cross-platform validation

Recorded checks: 27 October 2026, source
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
