## Resubmission

This is a resubmission. In this version I have:

* Replaced relative links in README.md to files excluded from the package
  build (benchmarks/README.md, INTERFACE.md, LICENSE.md,
  validation/README.md) with absolute GitHub URLs.

## Submission

This is the first submission of RcppAMGCL (0.1.0).

## Test environments

* Windows 11: R-devel 4.7.0 (2026-10-26 r90766), OpenMP enabled.
* GitHub Actions: Ubuntu (R 4.3.3, release, devel), Windows (release,
  devel), and macOS (release), each with OpenMP enabled and disabled.

Tests and examples use at most two threads.

## R CMD check results

0 errors, 0 warnings, 1 NOTE:

* New submission.

The full `R CMD check --as-cran` run included online incoming checks,
tests, examples, and PDF/HTML manuals.
All 12 GitHub Actions configurations passed.

There are no known CRAN reverse dependencies.
