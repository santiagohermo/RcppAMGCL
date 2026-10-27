# RcppAMGCL 0.1.0

* Add AMG-preconditioned CG and BiCGStab for general `Matrix::dgCMatrix`
  and symmetric `Matrix::dsCMatrix` inputs.
* Accept finite real numeric vector and one-column matrix right-hand sides.
* Add one-off solves and reusable solver objects with owned matrix storage.
* Add optional OpenMP setup and solve acceleration, with one thread by
  default and scoped thread settings.
* Add build capability reporting and an R-independent C++ interface.
* Explain direct coarse LU preconditioning for systems with at most 64 rows;
  convergence does not establish positive definiteness.
* Include setup context in AMGCL hierarchy errors while retaining the
  underlying diagnostic.
* Bundle pinned AMGCL 1.4.5 with its license notices.
