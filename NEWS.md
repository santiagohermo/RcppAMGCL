# RcppAMGCL 0.1.0

* Add AMG-preconditioned CG and BiCGStab for `Matrix::dgCMatrix` inputs.
* Add one-off solves and reusable solver objects with owned matrix storage.
* Add optional OpenMP setup and solve acceleration, with one thread by
  default and scoped thread settings.
* Add build capability reporting and an R-independent C++ interface.
* Bundle pinned AMGCL 1.4.5 with its license notices.
