# RcppAMGCL

RcppAMGCL solves sparse linear systems in R with AMGCL's built-in CPU
backend. Version 0.1.0 provides AMG-preconditioned conjugate gradients (CG)
and BiCGStab, reusable solver objects, and optional OpenMP acceleration for
both setup and solving. It also installs an R-independent C++ interface for
use in downstream projects.

The initial scope is a foundation for larger sparse estimation tasks such
as AKM. Singular systems, identification constraints, fixed-effect
construction, and inference are outside this version's scope.

## Install from source

Install R 4.3.0 or later, the `Rcpp` and `Matrix` packages, and the compiler toolchain for
your platform. Rtools is required on Windows. From the repository directory:

```sh
R CMD build .
R CMD INSTALL RcppAMGCL_0.1.0.tar.gz
```

No Boost dependency is required. AMGCL 1.4.5 is bundled at commit
`dde264f0f8d5d0a292a3bb3709952cc4e004604e`; its MIT license is retained in
[`inst/licenses/AMGCL-LICENSE.md`](inst/licenses/AMGCL-LICENSE.md).

OpenMP is detected with R's compiler flags during installation. The
`RCPPAMGCL_OPENMP` environment variable accepts `auto` (the default), `on`
(require OpenMP support), or `off` (compile a serial build). For example,
in a POSIX shell:

```sh
RCPPAMGCL_OPENMP=off R CMD INSTALL RcppAMGCL_0.1.0.tar.gz
RCPPAMGCL_OPENMP=on R CMD INSTALL RcppAMGCL_0.1.0.tar.gz
```

In PowerShell, set `$env:RCPPAMGCL_OPENMP = "off"` before running the install
command. A serial build exposes the same R interface.

## Solve in R

```r
library(RcppAMGCL)

A <- methods::as(Matrix::bandSparse(
  1000, k = c(-1, 0, 1),
  diagonals = list(rep(-1, 999), rep(3, 1000), rep(-1, 999))
), "dgCMatrix")
b <- rep(1, nrow(A))

fit <- amgcl_solve(A, b)
fit$converged
fit$residual
head(fit$x)

# Reuse the matrix and AMG hierarchy for another right-hand side.
s <- amgcl_setup(A, solver = "cg", tol = 1e-8, maxiter = 1000)
fit1 <- amgcl_solve(s, b)
fit2 <- amgcl_solve(s, seq_len(nrow(A)))

# For a nonsymmetric system, select BiCGStab.
B <- A + Matrix::sparseMatrix(i = 1, j = 2, x = 0.2, dims = dim(A))
fit3 <- amgcl_solve(B, b, solver = "bicgstab")
```

Inputs must be nonempty, square `Matrix::dgCMatrix` matrices with finite
entries and a nonzero diagonal entry in every row, and finite numeric
right-hand-side vectors. The selected AMG hierarchy and unpivoted coarse
solver require the nonzero diagonal. CG requires a
symmetric positive definite matrix; positive definiteness is not checked
automatically. Normalize or otherwise remove the nullspace of an AKM
system before solving. Convergence is not guaranteed for arbitrary
matrices.

Every solve returns `x`, `iterations`, AMGCL's residual estimate `error`, an
independently computed relative Euclidean `residual`, `converged`, and the
effective requested thread cap `threads`. For a zero right-hand side,
`residual` is the absolute residual norm. `converged` compares this
independent residual with the configured tolerance. Inspect these fields
before relying on a solution.

Reusable solvers own a copy of the matrix and retain their method,
tolerance and iteration limit. Create a new solver to change those
settings. Their external pointers are process-local: do not save/restore
them, send them to another R process, or share one object concurrently.
Recreate the solver from its matrix in each process. R releases the native
storage when the object is garbage collected.

## Thread control

```r
amgcl_info()
s <- amgcl_setup(A, threads = 2)
fit <- amgcl_solve(s, b, threads = 2)
```

Both entry points default to `threads = 1`. Set `threads` explicitly for
each operation that should use parallelism; a reusable solver does not
retain the setup thread count. The package caps requests by processor
capacity and `OMP_THREAD_LIMIT`. Set that environment variable before
starting R. A serial build always uses one thread.

OpenMP thread-count and dynamic-thread settings are scoped to each
operation and restored afterwards, including on errors. Worker threads
never call the R API. The result's `threads` field is a requested cap, not
a measured team size; runtime constraints may cause OpenMP to use fewer
workers.

`amgcl_info()` reports `openmp`, `max_threads`, `thread_limit`,
`current_threads`, and `amgcl_version`. Tests and examples use at most two
threads. The separate benchmark script under `tools/` measures setup and
solve times for 1, 2, 4, and 8 threads on larger sparse matrices.

```sh
Rscript tools/benchmark.R benchmarks/results-local.csv 3 160000
```

Recorded measurements and reproduction details are in
[`benchmarks/README.md`](benchmarks/README.md).

## C++ use

The installed `RcppAMGCL.h` exposes the `rcppamgcl` namespace with owned CSR
matrix storage, reusable solvers, and scoped thread controls. The header
does not depend on R or Rcpp. Downstream R packages can add `RcppAMGCL` to
`LinkingTo` and include:

```cpp
#include <RcppAMGCL.h>

rcppamgcl::Matrix A;
A.n = 2;
A.ptr = {0, 2, 4};       // CSR row pointers
A.col = {0, 1, 0, 1};    // Zero-based columns, sorted within each row
A.val = {2, -1, -1, 2};

rcppamgcl::Solver solver(std::move(A), "cg", 1e-8, 1000, 1);
auto result = solver.solve(std::vector<double>{1, 0}, 1);
// result.x is {2/3, 1/3}; inspect result.converged and result.residual.
```

The package uses R's default C++ compiler standard. The installed header
and bundled AMGCL code are also compatible with C++14. Downstream code must
enable and link OpenMP in its own build to use it there; the installed R
package's build mode does not automatically set compiler flags for another
package. See the installed header for the C++ API. A standalone project
can include the installed package's `include` directory without linking
to the R package shared library.

## Development

```sh
R CMD build .
R CMD check --as-cran RcppAMGCL_0.1.0.tar.gz
```

GitHub Actions covers Linux, Windows, and macOS with OpenMP enabled and
disabled. These jobs explicitly check the installed build's OpenMP
capability. The repository is initially private; CI run results, rather
than the presence of the workflow, establish whether each platform passes.
The `URL` and `BugReports` fields are omitted from `DESCRIPTION` while the
repository is private, because unauthenticated CRAN checks cannot access
it. Add the public repository and issue tracker URLs when publishing.

See [`validation/README.md`](validation/README.md) for executed checks,
the tested source revision, and the remaining CRAN submission note.

The R wrapper is MIT licensed. Bundled AMGCL retains its upstream license
and copyright notices.
