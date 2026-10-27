#include <Rcpp.h>
#include <RcppAMGCL.h>

namespace {

SEXP solver_tag() {
    return Rf_install("RcppAMGCL::Solver");
}

rcppamgcl::Matrix matrix_from_R(const Rcpp::S4& A) {
    if (!A.is("dgCMatrix"))
        Rcpp::stop("A must be a Matrix::dgCMatrix");
    SEXP dimensions = A.slot("Dim");
    SEXP pointers = A.slot("p");
    SEXP indices = A.slot("i");
    SEXP values = A.slot("x");
    if (TYPEOF(dimensions) != INTSXP || Rf_xlength(dimensions) != 2 ||
            TYPEOF(pointers) != INTSXP || TYPEOF(indices) != INTSXP ||
            TYPEOF(values) != REALSXP)
        Rcpp::stop("Invalid dgCMatrix slot types");
    const int *dim = INTEGER(dimensions);
    if (dim[0] < 1 || dim[0] != dim[1])
        Rcpp::stop("A must be a nonempty square matrix");
    if (Rf_xlength(indices) != Rf_xlength(values))
        Rcpp::stop("Invalid dgCMatrix slot lengths");
    // All R objects are inspected here, on the calling thread. Conversion and
    // the solver then use only owned C++ data, including all OpenMP workers.
    return rcppamgcl::Matrix::from_csc(static_cast<std::size_t>(dim[0]),
        INTEGER(pointers), static_cast<std::size_t>(Rf_xlength(pointers)),
        INTEGER(indices), REAL(values), static_cast<std::size_t>(Rf_xlength(values)));
}

rcppamgcl::Solver& solver_from_pointer(SEXP pointer) {
    if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrTag(pointer) != solver_tag())
        Rcpp::stop("Invalid AMGCL solver pointer");
    if (R_ExternalPtrAddr(pointer) == nullptr)
        Rcpp::stop("AMGCL solver pointer is null; recreate objects after serialization");
    return *static_cast<rcppamgcl::Solver *>(R_ExternalPtrAddr(pointer));
}

} // namespace

// [[Rcpp::export(rng = false)]]
SEXP amgcl_setup_cpp(Rcpp::S4 A, std::string solver, double tol,
                     int maxiter, int threads) {
    rcppamgcl::Matrix matrix = matrix_from_R(A);
    Rcpp::checkUserInterrupt();
    std::unique_ptr<rcppamgcl::Solver> owner(new rcppamgcl::Solver(
        std::move(matrix), solver, tol, maxiter, threads));
    Rcpp::XPtr<rcppamgcl::Solver> pointer(owner.get(), true, solver_tag());
    owner.release();
    return pointer;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List amgcl_solve_cpp(SEXP pointer, Rcpp::NumericVector b, int threads) {
    rcppamgcl::Solver& solver = solver_from_pointer(pointer);
    std::vector<double> rhs(b.begin(), b.end());
    Rcpp::checkUserInterrupt();
    const rcppamgcl::SolveResult result = solver.solve(rhs, threads);
    Rcpp::checkUserInterrupt();
    return Rcpp::List::create(
        Rcpp::Named("x") = Rcpp::wrap(result.x),
        Rcpp::Named("iterations") = static_cast<int>(result.iterations),
        Rcpp::Named("error") = result.error,
        Rcpp::Named("residual") = result.residual,
        Rcpp::Named("converged") = result.converged,
        Rcpp::Named("threads") = result.threads);
}

// [[Rcpp::export(rng = false)]]
Rcpp::List amgcl_info_cpp() {
    const rcppamgcl::ThreadInfo info = rcppamgcl::thread_info();
    return Rcpp::List::create(
        Rcpp::Named("openmp") = info.openmp,
        Rcpp::Named("max_threads") = info.max_threads,
        Rcpp::Named("thread_limit") = info.thread_limit,
        Rcpp::Named("current_threads") = info.current_threads,
        Rcpp::Named("amgcl_version") = rcppamgcl::amgcl_version());
}
