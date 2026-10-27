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
The package supports double-precision compressed sparse column matrices
from \pkg{Matrix}. See \code{\link{amgcl_info}} for build capabilities.
The installed \file{RcppAMGCL.h} header provides an R-independent C++
interface for downstream code; the bundled AMGCL headers retain their
upstream license notices.
}
