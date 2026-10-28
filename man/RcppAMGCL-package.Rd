\name{RcppAMGCL-package}
\alias{RcppAMGCL-package}
\alias{RcppAMGCL}
\docType{package}
\title{Sparse Iterative Solvers with Algebraic Multigrid Preconditioning}
\description{
An R interface to the header-only AMGCL library using Rcpp, with a built-in
CPU backend and optional OpenMP acceleration. Use
\code{\link{amgcl_solve}} for a single sparse system, or
\code{\link{amgcl_setup}} to reuse a preconditioner across right-hand sides.
}
\details{
The package supports general and symmetric double-precision compressed
sparse column matrices (\code{dgCMatrix} and \code{dsCMatrix}) from
\pkg{Matrix}, with vector or one-column matrix right-hand sides.
See \code{\link{amgcl_info}} for build capabilities.
The installed \file{RcppAMGCL.h} header provides an R-independent C++
interface for downstream code; the bundled AMGCL headers retain their
upstream license notices.
This product includes Kratos Multi-Physics technology. The installed
\file{COPYRIGHTS} file inventories bundled-code ownership and licenses;
additional notices are installed under \file{licenses}.
Use \code{citation("RcppAMGCL")} for references to this package and AMGCL.
}
\references{
Demidov, D. (2020). AMGCL -- A C++ library for efficient solution of large
sparse linear systems. \emph{Software Impacts}, 6, 100037.
\doi{10.1016/j.simpa.2020.100037}.
}
