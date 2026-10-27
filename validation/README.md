# Validation

Latest recorded checks: 27 October 2026, source
[`771186f`](https://github.com/santiagohermo/RcppAMGCL/commit/771186f44af6b215d5844e8b6e729a3e94a6dcc4).

## Local check

`R CMD check --as-cran` on Windows 11, R 4.5.3 and Rtools GCC 14.3.0:
**0 errors, 0 warnings, 1 NOTE: New submission**. Tests, examples and both
PDF and HTML manuals passed. The serial build also compiled and passed
the R tests. See the [check log](windows-review-as-cran.log).

Tested source archive SHA256:
`ba0b0a9a92efc0c31bcc5998041f1c5640c65c6f4d9fb1ebc108c6d3b75afc4b`.

## GitHub Actions

All six configurations in [this run](https://github.com/santiagohermo/RcppAMGCL/actions/runs/40431758912)
passed with `Status: OK` using R 4.6.1.

| Platform | Compiler | OpenMP on | OpenMP off |
| --- | --- | --- | --- |
| Ubuntu 24.04 | GCC 13.3.0 | Passed | Passed |
| Windows | Rtools GCC 14.3.0 | Passed | Passed |
| macOS | Apple clang 21.0.0 | Passed | Passed |

Each job installs the package, verifies the selected backend, runs R and
standalone C++14 tests, and runs `R CMD check --as-cran --no-manual` with
remote incoming checks disabled. The local check includes the manual and
incoming checks. Linux installation was retried after slow system-package
downloads; the successful results are in attempt 2.

Tests cover both solvers, symmetric and nonsymmetric matrices, reusable
objects, input validation, ownership, extreme RHS scales, independent
residuals, and thread limits/restoration. Tests and examples use at most
two threads. [Performance benchmarks](../benchmarks/README.md) run separately.
