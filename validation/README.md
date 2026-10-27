# Executed validation

Validated on 27 October 2026. The tested package source is commit
[`ec8bd80`](https://github.com/santiagohermo/RcppAMGCL/commit/ec8bd801c6583a26601aee9140547be80111d680).
Later documentation and benchmark-record changes do not alter package code.

## Full local CRAN check

Executed `R CMD check --as-cran RcppAMGCL_0.1.0.tar.gz` on Windows 11,
R 4.5.3, Rtools GCC 14.3.0, with OpenMP enabled and test threads limited
to two. Result: **0 errors, 0 warnings, 1 NOTE**. The only note is
**New submission** from CRAN incoming feasibility. Both PDF and HTML
manual checks passed. See the [complete check log](windows-as-cran.log).

The tested source archive has SHA256:

```text
d84f28e339aef9603bc4768cdb92b42cc506833cba1fb928d749a5f7612e3f38
```

Windows `R CMD build` reports that it corrected Unix execute permissions
on `configure` and `cleanup`; the resulting archive passed its checks.
Both files also have executable bits recorded in Git.

## Cross-platform CI

All six jobs in [the final CI run](https://github.com/santiagohermo/RcppAMGCL/actions/runs/40425990481)
completed successfully. Each check log ended with `Status: OK`.

| Platform | R | Compiler | OpenMP on | OpenMP off |
| --- | --- | --- | --- | --- |
| Ubuntu 24.04 | 4.6.1 | GCC 13.3.0 | Passed | Passed |
| Windows | 4.6.1 | Rtools GCC 14.3.0 | Passed | Passed |
| macOS | 4.6.1 | Apple clang 21.0.0 | Passed | Passed |

Each job installed the actual package, verified `amgcl_info()$openmp`,
ran the R tests, compiled and ran the standalone C++ test, and ran
`R CMD check --as-cran --no-manual`. CI disables remote incoming checks;
the full local check above also covers incoming checks and the PDF manual.
The macOS OpenMP job used Homebrew libomp and explicit compiler/linker
configuration. CRAN tests, examples, and standalone CI tests use at most
two threads.

## Coverage

R tests cover both solvers on a 1,024-unknown sparse Poisson system,
nonsymmetric BiCGStab, one-off and reusable solves, independent residuals,
zero RHS, iteration exhaustion, matrix ownership after mutation/GC,
malformed input, foreign/null/serialized pointers, and a fresh process
with `OMP_THREAD_LIMIT=1`.

Standalone C++ tests compile without R headers or the R runtime. They
cover C++14 compatibility, actual OpenMP team sizes, numerical correctness,
tiny and large RHS values, nonsymmetric CSC conversion, and restoration
of thread-count and dynamic settings after success and setup exceptions.
On every platform the OpenMP executable observed teams of two, then one
under `OMP_THREAD_LIMIT=1`; serial executables observed one in both runs.

The [benchmark records](../benchmarks/README.md) measure setup and solving
separately at 1, 2, 4, and 8 threads. These are developer benchmarks,
excluded from CRAN tests and the source package.

There are no remaining test failures or check warnings. This validation
does not add singular-system identification or AKM estimation to the
minimal solver package.
