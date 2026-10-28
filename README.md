

<!--
README.md is generated from README.qmd.
With Quarto and the R packages knitr and rmarkdown installed, run:
    quarto render README.qmd
Rendering reads the recorded CSV; it does not rerun the benchmarks.
The R usage examples are displayed without execution.
-->

# RcppAMGCL

A minimal Rcpp wrapper for AMGCL: sparse CG and BiCGStab solvers with
algebraic multigrid preconditioning, reusable solver objects, optional
OpenMP, and an R-independent C++ interface.

## Install

Requires R 4.3.0 or later, `Rcpp`, `Matrix`, and an R compiler toolchain
(Rtools on Windows). From this repository:

``` sh
R CMD build .
R CMD INSTALL RcppAMGCL_0.1.0.tar.gz
```

OpenMP is detected automatically. Set `RCPPAMGCL_OPENMP=on` to require
it or `RCPPAMGCL_OPENMP=off` to build without it. Both builds expose the
same API.

## Use from R

``` r
library(RcppAMGCL)

A <- Matrix::bandSparse(1000, k = c(0, 1),
  diagonals = list(rep(3, 1000), rep(-1, 999)), symmetric = TRUE)
b <- rep(1, nrow(A))

fit <- amgcl_solve(A, b)
fit$x
fit$converged
fit$residual

# Reuse setup for another right-hand side.
s <- amgcl_setup(A, solver = "cg", tol = 1e-8)
amgcl_solve(s, 2 * b)$x
```

`A` accepts `dgCMatrix` or `dsCMatrix`; `b` accepts a real numeric
vector or one-column numeric matrix. Use `solver = "bicgstab"` for
nonsymmetric systems. CG requires a symmetric positive definite matrix,
which is not verified automatically. Matrices must be square, finite,
and have nonzero diagonal entries. Check `converged` and the
independently computed `residual` before using a result. Reusable
solvers are process-local and cannot be saved and restored. See
`?amgcl_setup` and `?amgcl_solve` for the complete input and lifetime
contracts.

## Threads

``` r
amgcl_info()
s <- amgcl_setup(A, threads = 2)
fit <- amgcl_solve(s, b, threads = 2)
```

Each operation defaults to one thread. Requests are capped by processor
capacity and `OMP_THREAD_LIMIT`; serial builds use one thread. OpenMP
settings are restored after each operation. Set `OMP_THREAD_LIMIT`
before starting R.

## Performance

Median seconds on an Intel Core Ultra 7 265H, Windows 11, R 4.5.3 and
Matrix 1.7-4 with R’s supplied BLAS. Both methods reuse setup for ten
distinct right-hand sides; three repetitions per configuration.

| Poisson system | Method | Threads | Setup + 1 RHS (s) | Setup + 10 RHS (s) |
|:---|:---|---:|---:|---:|
| 2D, 100,489 unknowns | AMGCL CG | 1 | 0.237 | 1.397 |
| 2D, 100,489 unknowns | AMGCL CG | 8 | 0.126 | 0.684 |
| 2D, 100,489 unknowns | Matrix Cholesky | 1 | 0.321 | 0.459 |
| 3D, 103,823 unknowns | AMGCL CG | 1 | 0.348 | 2.109 |
| 3D, 103,823 unknowns | AMGCL CG | 8 | 0.087 | 0.462 |
| 3D, 103,823 unknowns | Matrix Cholesky | 1 | 19.857 | 20.516 |

All solves met a relative residual target of `1e-8`; Cholesky achieved
much tighter accuracy. AMGCL had cheaper setup in these cases, while
Cholesky was faster for ten RHS in 2D. These deterministic grid problems
suit AMG; results depend on the matrix, BLAS, hardware and accuracy
target.

See [benchmark methods and results](benchmarks/README.md) for matrix
sizes, setup and solve timings, totals, and reproduction commands.

## Use from C++

Add `RcppAMGCL` to a downstream R package’s `LinkingTo` field and
include `RcppAMGCL.h`. The header also works in standalone C++ without R
or Rcpp:

``` cpp
#include <RcppAMGCL.h>

rcppamgcl::Matrix A;
A.n = 2;
A.ptr = {0, 2, 4};       // CSR row pointers
A.col = {0, 1, 0, 1};    // Sorted, zero-based column indices
A.val = {2, -1, -1, 2};
rcppamgcl::Solver solver(std::move(A), "cg", 1e-8, 1000, 1);
auto result = solver.solve(std::vector<double>{1, 0}, 1);
```

The header is compatible with C++14. Downstream builds supply their own
OpenMP compiler and linker flags; no link to the R package library is
required. See the [API and build contract](INTERFACE.md).

AMGCL [1.4.5 is pinned](inst/AMGCL-VERSION) and bundled with its [MIT
license](inst/licenses/AMGCL-LICENSE.md); Boost is not required. The R
wrapper is also [MIT licensed](LICENSE.md). See [validation](validation/README.md) for
executed checks and platform results.

This product includes Kratos Multi-Physics technology. See the
[bundled-code copyright inventory](inst/COPYRIGHTS) and
[third-party notices](inst/licenses/AMGCL-THIRD-PARTY-NOTICES.txt).

For publication references to this package and AMGCL, run
`citation("RcppAMGCL")` in R. The AMGCL reference is Demidov (2020),
[AMGCL -- A C++ library for efficient solution of large sparse linear systems](https://doi.org/10.1016/j.simpa.2020.100037).

## Disclaimer

The original AMGCL solver code was developed by Denis Demidov and
upstream contributors in
[ddemidov/amgcl](https://github.com/ddemidov/amgcl), not in this
repository. RcppAMGCL brings that library to R through Rcpp, bundling
AMGCL 1.4.5 with its original license notices and the compatibility
patches documented in [inst/AMGCL-VERSION](inst/AMGCL-VERSION).

GPT-6 Astra, working with delegated coding agents, implemented the
RcppAMGCL wrapper and its R and C++ interfaces, build configuration,
tests, documentation, and benchmark tooling. Fable 5.1 provided review
and feedback on the package.
