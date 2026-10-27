# RcppAMGCL 0.1.0

* Add AMG-preconditioned CG and BiCGStab for general `Matrix::dgCMatrix`
  and symmetric `Matrix::dsCMatrix` inputs.
* Accept finite real numeric vector and one-column matrix right-hand sides.
* Add one-off solves and reusable solver objects with owned matrix storage.
* Compute independent residuals from the hierarchy's owned matrix without
  retaining a second CSR copy; preserve 64-bit C++ indices.
* Add optional OpenMP setup and solve acceleration, with one thread by
  default and scoped thread settings.
* Parallelize independent residuals on large systems using fixed blocks and
  overflow-resistant norms, with the same serial fallback and thread cap.
* Add build capability reporting and an R-independent C++ interface.
* Include setup context in AMGCL hierarchy errors while retaining the
  underlying diagnostic.
* Bundle pinned AMGCL 1.4.5 with its license notices.
