# API and build contract

## R interface

```r
amgcl_setup(A, solver = c("cg", "bicgstab"), tol = 1e-8,
            maxiter = 1000L, threads = 1L)
amgcl_solve(A, b, solver = c("cg", "bicgstab"), tol = 1e-8,
            maxiter = 1000L, threads = 1L)
amgcl_info()
```

`A` is a nonempty square `dgCMatrix` or `dsCMatrix` with finite values and
nonzero diagonal entries. Symmetric storage is validated before expansion.
`b` is a finite real numeric vector, a one-column base numeric matrix, or
a one-column double-precision `Matrix` object of matching size. Logical,
complex, multiple-column, and higher-dimensional inputs are rejected.

Setup returns an `amgcl_solver` containing a tagged external pointer and
metadata `n`, `solver`, `tol`, and `maxiter`. Solve accepts either a matrix
or this reusable object. Its method, tolerance and iteration limit are
fixed at setup; `threads` applies separately to each operation. Objects own
their matrix storage and are process-local; serialization and concurrent
use of one object are unsupported.

Solve returns `x`, `iterations`, AMGCL's estimate `error`, independently
computed relative Euclidean `residual` (absolute for zero `b`),
`converged = residual <= tol`, and the effective requested cap `threads`.
CG requires a symmetric positive definite matrix; convergence does not
validate that assumption. Systems with at most 64 rows use direct coarse
LU preconditioning and often converge in one iteration. Singular systems
and nullspace constraints are not handled automatically.

`amgcl_info()` returns `openmp`, `max_threads`, `thread_limit`,
`current_threads`, and `amgcl_version`. Thread requests are capped by
processor capacity and `OMP_THREAD_LIMIT`. Calling-task OpenMP thread and
dynamic settings are restored on success or error. Worker threads never
call the R API. CRAN tests and examples use at most two threads.

## C++ interface and build

The R-independent `RcppAMGCL.h` exposes `rcppamgcl::Matrix` (owned CSR),
`rcppamgcl::Solver`, `Solver::solve()`, and `thread_info()`. The built-in CPU
backend uses smoothed aggregation and SPAI0 relaxation. The hierarchy owns
the original matrix used for independent residuals; no second solver-level
CSR copy is retained. C++ indices retain 64-bit capacity. Setup exceptions
include context and preserve the underlying AMGCL diagnostic.

Native Rcpp exports are `amgcl_setup_cpp(A, solver, tol, maxiter, threads)`,
`amgcl_solve_cpp(pointer, b, threads)`, and `amgcl_info_cpp()`. The R layer
passes validated general sparse storage and a double vector to these
functions. Native code validates matrix storage and pointer tags/addresses.

Requires R >= 4.3.0 with R's default C++ compiler standard; the installed
header also supports C++14. Rcpp and Matrix are the only non-base R
dependencies. AMGCL 1.4.5 is pinned in `inst/AMGCL-VERSION`, with license
notices retained and `AMGCL_NO_BOOST` enabled.

`RCPPAMGCL_OPENMP=auto|on|off` controls the compile/link probe using R's
`SHLIB_OPENMP_CXXFLAGS`. `auto` falls back to serial; `on` requires OpenMP.
Downstream builds select their own compiler and OpenMP flags. CI covers
Linux, Windows and macOS with OpenMP enabled and disabled; executed results
are recorded under `validation/`. Benchmarks are separate from CRAN tests.
