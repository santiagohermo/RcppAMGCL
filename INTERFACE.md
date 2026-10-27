# Implementation contract

Package RcppAMGCL 0.1.0; R >= 4.3.0 with R's default C++ compiler standard;
installed C++ header also compatible with C++14. Rcpp and Matrix are the only
non-base R dependencies (methods is a base R package).
Vendor AMGCL 1.4.5 pinned to its resolved tag commit, retain licence; define AMGCL_NO_BOOST.

Public R API:
- `amgcl_setup(A, solver = c("cg", "bicgstab"), tol = 1e-8, maxiter = 1000L, threads = 1L)` accepts dgCMatrix or dsCMatrix and returns class `amgcl_solver` wrapping a tagged external pointer and metadata. Validate symmetric storage before expansion with `methods::as(A, "generalMatrix")`; native setup receives dgCMatrix.
- `amgcl_solve(A, b, solver = c("cg", "bicgstab"), tol = 1e-8, maxiter = 1000L, threads = 1L)` accepts dgCMatrix, dsCMatrix or amgcl_solver as A. Reusable object owns solver/tol/maxiter; threads applies to each operation. Returns list x, iterations, error (AMGCL estimate), residual (independently computed relative 2-norm, absolute for zero b), converged, threads (effective requested cap).
- `b` accepts a finite real numeric vector (including integer), a one-column base numeric matrix, or a one-column real double Matrix object. Validate shape and type before flattening; reject multiple columns, arrays with more than two dimensions, and logical/complex inputs. Its length or row count must equal the coefficient matrix dimension.
- `amgcl_info()` returns openmp logical, max_threads capacity (min num processors, thread limit; >=1), thread_limit, current_threads (calling-task setting), amgcl_version.
- print.amgcl_solver is concise; objects process-local, cannot be saved/restored or shared concurrently.

Native Rcpp exports (root generates RcppExports):
- `amgcl_setup_cpp(Rcpp::S4 A, std::string solver, double tol, int maxiter, int threads)` returns tagged external pointer (R layer wraps as list pointer, n, solver, tol, maxiter).
- `amgcl_solve_cpp(SEXP pointer, Rcpp::NumericVector b, int threads)` returns result list above.
- `amgcl_info_cpp()` returns info list above.

C++ integration owns src/solver.cpp and inst/include/RcppAMGCL.h (R-independent reusable interface). Prefer smoothed aggregation with SPAI0 for symmetric CG-compatible preconditioner. Correct CSC->CSR conversion, owned storage, strict validation, null/tag pointer checks, finite result checks. Thread guard saves/restores calling task's OpenMP num_threads and dynamic settings even on exceptions; cap threads by OMP_THREAD_LIMIT and processor capacity, serial fallback. No R API in parallel regions.

Systems with at most 64 rows use direct coarse LU preconditioning and often
converge in one iteration; convergence never establishes symmetry or positive
definiteness. Hierarchy setup errors add context while preserving AMGCL's
diagnostic, without classifying every failure as singularity.

R interface agent owns DESCRIPTION, NAMESPACE, R/ (except generated RcppExports.R), man/, README.md, LICENSE, NEWS.md. Author from git: Santiago Hermo <santiago.hermo@monash.edu>. MIT + file LICENSE; acknowledge vendored author in Authors@R. Documentation must state CG requires SPD; no singular-system handling yet (AKM must normalize/drop nullspace).

Infrastructure agent owns configure, configure.win, cleanup, cleanup.win, src/Makevars.in, tests/, .github/, .Rbuildignore, .gitignore, tools/ benchmark/check helpers. Use R's SHLIB_OPENMP_CXXFLAGS, compile/link probe; RCPPAMGCL_OPENMP=auto|on|off, optional fallback. Six CI combinations Linux/Windows/macOS x on/off; explicitly verify capability; CRAN tests/examples <=2 threads. Base R test script preferred to avoid dependency. Benchmark setup/solve separately 1,2,4,8 threads, default >=100k unknowns, outside tests.

Root owns vendoring, local environment/build, generated exports, Git repository/remote, integration review and fixes. Coordinate before changing this interface.
