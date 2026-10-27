# Executed validation

Validated on 27 October 2026. The tested package source is commit
[`771186f`](https://github.com/santiagohermo/RcppAMGCL/commit/771186f44af6b215d5844e8b6e729a3e94a6dcc4),
including the review fixes. Later changes to excluded benchmark tools and
validation records do not alter the tested package source.

## Full local CRAN check

Executed `R CMD check --as-cran RcppAMGCL_0.1.0.tar.gz` on Windows 11,
R 4.5.3, Rtools GCC 14.3.0, with OpenMP enabled and test threads limited
to two. Result: **0 errors, 0 warnings, 1 NOTE**. The only note is
**New submission** from CRAN incoming feasibility. Both PDF and HTML
manual checks passed. See the [complete check log](windows-review-as-cran.log).
The [initial implementation's log](windows-as-cran.log) is retained separately.

The tested source archive has SHA256:

```text
ba0b0a9a92efc0c31bcc5998041f1c5640c65c6f4d9fb1ebc108c6d3b75afc4b
```

Windows `R CMD build` reports that it corrected Unix execute permissions
on `configure` and `cleanup`; the resulting archive passed its checks.
Both files also have executable bits recorded in Git.

The final source archive also compiled with `RCPPAMGCL_OPENMP=off` locally.
Its complete R test script passed and `amgcl_info()$openmp` was false.
Regenerating the Rcpp bindings produced no changes.

## Cross-platform CI

All six configurations in [this CI run](https://github.com/santiagohermo/RcppAMGCL/actions/runs/40431758912)
completed successfully. Each check log ended with `Status: OK`.

The initial Linux jobs were canceled during slow Ubuntu system-package
downloads and retried on the same package source. The Windows and macOS
successes were retained. Attempt 2 contains the completed retry results.

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

The review adds both triangle layouts of symmetric `dsCMatrix` inputs,
natural `crossprod()` output, dense/sparse one-column RHS inputs, integer
RHS values, malformed symmetric storage, missing/zero diagonals, fixed
solver parameter enforcement, print behavior, and contextual native setup
errors. A 32,769-row system exercises the parallel residual and its final
partial block at RHS scales `1e-200`, `1`, and `1e200`, plus zero RHS.

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
