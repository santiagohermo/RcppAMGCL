\name{amgcl_solve}
\alias{amgcl_solve}
\title{Solve a Sparse Linear System}
\usage{
amgcl_solve(A, b, solver = c("cg", "bicgstab"), tol = 1e-8,
            maxiter = 1000L, threads = 1L)
}
\arguments{
  \item{A}{A nonempty, square \code{Matrix::dgCMatrix} or symmetric
    \code{Matrix::dsCMatrix} with finite entries and a nonzero diagonal
    entry in every row, or a reusable
    \code{amgcl_solver} created by \code{\link{amgcl_setup}}.}
  \item{b}{A finite real numeric right-hand-side vector, including an
    integer vector, of length equal to the matrix dimension. A one-column
    base numeric matrix or one-column double-precision \pkg{Matrix} object
    with the same number of rows is also accepted and converted to a vector.
    Multiple columns, arrays with more than two dimensions, and logical or
    complex inputs are rejected.}
  \item{solver}{\code{"cg"} or \code{"bicgstab"}. Must be omitted when
    \code{A} is a reusable solver.}
  \item{tol}{Finite positive tolerance for AMGCL's relative residual
    estimate. Must be omitted when \code{A} is a reusable solver.}
  \item{maxiter}{Positive integer iteration limit. Must be omitted when
    \code{A} is a reusable solver.}
  \item{threads}{Positive integer requested thread count for this call,
    defaulting to one. Capped by the capacity reported by
    \code{\link{amgcl_info}}; one thread is used without OpenMP.}
}
\description{
Solve \eqn{Ax = b} with AMG-preconditioned conjugate gradients or BiCGStab,
either by constructing a solver for one call or reusing an existing solver.
}
\details{
Every solve starts from a zero initial guess. Inspect \code{converged} and
\code{residual} before using a result. CG requires a symmetric positive
definite matrix; convergence does not validate that assumption. See
\code{\link{amgcl_setup}} for matrix requirements, the direct coarse LU
preconditioner used for systems with at most 64 rows, and setup errors.

When \code{A} is a matrix, \code{threads} controls both setup and solving.
When it is a reusable solver, \code{threads} controls this solve only; the
setup thread count is not retained. OpenMP settings are restored after each
operation. See \code{\link{amgcl_setup}} for thread and object-lifetime
details.
}
\value{
A list with components:
\describe{
  \item{x}{Numeric solution vector.}
  \item{iterations}{Number of iterations reported by AMGCL.}
  \item{error}{AMGCL's final relative residual estimate.}
  \item{residual}{Independently computed Euclidean residual norm divided
    by the Euclidean norm of \code{b}. For a zero right-hand side, the
    absolute residual norm is returned.}
  \item{converged}{Logical; whether the independently computed
    \code{residual} is at most \code{tol}.}
  \item{threads}{Effective requested thread cap after applying build and
    capacity limits. This is not a measurement of the worker count used
    inside each parallel region.}
}
}
\seealso{\code{\link{amgcl_setup}}, \code{\link{amgcl_info}}}
\examples{
A <- methods::as(Matrix::bandSparse(10, k = c(-1, 0, 1),
    diagonals = list(rep(-1, 9), rep(3, 10), rep(-1, 9))), "dgCMatrix")
result <- amgcl_solve(A, rep(1, 10))
result$x
result$residual

# Symmetric sparse storage and a one-column Matrix right-hand side.
S <- Matrix::bandSparse(10, k = c(0, 1),
    diagonals = list(rep(3, 10), rep(-1, 9)), symmetric = TRUE)
rhs <- S \%*\% rep(1, 10)
amgcl_solve(S, rhs)$x
}
