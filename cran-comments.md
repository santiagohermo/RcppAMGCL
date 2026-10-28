## Submission status

Prepared for the first CRAN submission of RcppAMGCL 0.1.0.
The repository is currently private: its URL and BugReports links return
HTTP 404 to anonymous visitors. Make the repository public and rerun the
online incoming checks on the final source archive before submitting.
This file is excluded from the source package.

## Package purpose and attribution

RcppAMGCL provides sparse CG and BiCGStab solvers with algebraic multigrid
preconditioning, reusable solver objects, optional OpenMP and an installed
C++ interface. It bundles AMGCL 1.4.5; provenance and local patches are
documented in inst/AMGCL-VERSION. The copyright inventory is in
inst/COPYRIGHTS, and the original notices are retained in the headers and
inst/licenses. This product includes Kratos Multi-Physics technology.

DESCRIPTION references Demidov (2020), DOI 10.1016/j.simpa.2020.100037.
The installed CITATION file provides references for the R package and AMGCL.

## Check environments and results

Checked on 28 October 2026 using R-devel 4.7.0 (2026-10-26 r90766) on
Windows 11, Rtools GCC 14.3.0, Rcpp 1.1.2, Matrix 1.7-6 and V8 8.2.0.
OpenMP was enabled with a two-thread cap. The complete
`R CMD check --as-cran` run included remote incoming checks, tests,
examples, and PDF/HTML manuals with math rendering enabled.

Result: 0 errors, 0 warnings, 1 NOTE. The incoming-feasibility NOTE
reports "New submission" and HTTP 404 for the private repository and
issue-tracker URLs. The DOI resolved successfully. The installed
`citation("RcppAMGCL")` returned both expected references.

Checked archive: RcppAMGCL_0.1.0.tar.gz, SHA256
`15440d6f16cbdf425fb0e5dc99172d0cc494d2dfc3c6ecf2afe8c46c602c3d80`.
The check log is retained in
validation/windows-r-devel-metadata-as-cran.log.

Earlier checks are documented in validation/README.md: both R-devel
backends passed on the preceding attribution-update archive, and six
R-release CI configurations passed on the earlier source revision.
Those historical results do not establish a check of this updated archive.

## Expected and outstanding notes

* "New submission" is expected for the first submission.
* The private GitHub URLs must become publicly accessible before submission.

## Reverse dependencies and resource use

This is a new package, with no known CRAN reverse dependencies.
Tests and examples use at most two threads. Larger benchmarks are excluded
from the source package. Package installation uses the bundled AMGCL headers
and does not download external solver code.
