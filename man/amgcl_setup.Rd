\name{amgcl_setup}
\alias{amgcl_setup}
\alias{amgcl_solver}
\alias{print.amgcl_solver}
\title{Create a Reusable Algebraic Multigrid Solver}
\usage{
amgcl_setup(A, solver = c("cg", "bicgstab"), tol = 1e-8,
            maxiter = 1000L, threads = 1L)

\method{print}{amgcl_solver}(x, ...)
}
\arguments{
  \item{A}{A nonempty, square \code{Matrix::dgCMatrix} containing finite
    double-precision values and a nonzero diagonal entry in every row.}
  \item{solver}{Iterative method: \code{"cg"} for conjugate gradients or
    \code{"bicgstab"} for stabilized bi-conjugate gradients.}
  \item{tol}{Finite positive convergence tolerance for AMGCL's relative
    residual estimate.}
  \item{maxiter}{Positive integer maximum number of iterations.}
  \item{threads}{Positive integer requested thread count, defaulting to one.
    Counts above the available capacity are capped. A build without OpenMP
    always uses one thread.}
  \item{x}{An \code{amgcl_solver} object.}
  \item{...}{Additional arguments, currently unused.}
}
\description{
Build an algebraic multigrid (AMG) preconditioner and retain it for repeated
solutions with the same matrix. The matrix is copied into storage owned by
the solver, so subsequent changes to the R matrix do not change the solver.
}
\details{
The package uses AMGCL's built-in CPU backend, smoothed aggregation and
SPAI0 relaxation. A nonzero diagonal is required by the selected AMG
hierarchy and its unpivoted coarse solver; matrices with a missing or zero
diagonal entry are rejected. Conjugate gradients require a symmetric positive definite
matrix. The package does not verify positive definiteness. BiCGStab accepts
nonsymmetric systems, but convergence is not guaranteed for every matrix.
Singular systems and nullspace constraints are not handled automatically;
for example, AKM models must be normalized before solving.

The solver owns its method, tolerance and iteration limit. These cannot be
overridden in subsequent calls to \code{amgcl_solve()}; create a new solver
to change them. The \code{threads} argument applies only to the current
operation: set it separately when solving with this object.

OpenMP, when available, accelerates both setup and solving. Thread counts
are capped by processor capacity and \env{OMP_THREAD_LIMIT}. OpenMP may
provide fewer worker threads than requested. The calling task's thread
count and dynamic-thread setting are restored after each operation. No R
API functions are called from worker threads.

Solver objects contain an R external pointer and are valid only in the
process where they were created. Do not save and restore them with
\code{saveRDS()}, send them to another R process, or access the same object
concurrently. Recreate a solver from its matrix in each process. Native
storage is released automatically when the object is garbage collected.
}
\value{
An object of class \code{amgcl_solver}, containing an external pointer and
metadata \code{n}, \code{solver}, \code{tol} and \code{maxiter}. Treat its
components as read-only. The print method returns the object invisibly.
}
\seealso{\code{\link{amgcl_solve}}, \code{\link{amgcl_info}}}
\examples{
A <- methods::as(Matrix::bandSparse(10, k = c(-1, 0, 1),
    diagonals = list(rep(-1, 9), rep(3, 10), rep(-1, 9))), "dgCMatrix")
s <- amgcl_setup(A)
print(s)
amgcl_solve(s, rep(1, 10))$x
amgcl_solve(s, seq_len(10))$residual
}
