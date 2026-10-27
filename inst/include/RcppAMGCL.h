#ifndef RCPPAMGCL_H
#define RCPPAMGCL_H

// R-independent, header-only C++ interface. Include this before AMGCL headers.
// Downstream builds choose OpenMP using their own compiler and linker flags.
#ifndef AMGCL_NO_BOOST
#define AMGCL_NO_BOOST
#endif
#ifndef AMGCL_RCPP_NO_OUTPUT
#define AMGCL_RCPP_NO_OUTPUT
#endif
#if defined(RCPPAMGCL_DISABLE_OPENMP) && defined(_OPENMP)
#error "RCPPAMGCL_DISABLE_OPENMP requires compiling without OpenMP support"
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/cg.hpp>

namespace rcppamgcl {

struct ThreadInfo {
    bool openmp;
    int max_threads;
    int thread_limit;
    int current_threads;
};

inline ThreadInfo thread_info() {
#ifdef _OPENMP
    const int limit = std::max(1, omp_get_thread_limit());
    return ThreadInfo{true, std::max(1, std::min(omp_get_num_procs(), limit)),
                      limit, std::max(1, omp_get_max_threads())};
#else
    return ThreadInfo{false, 1, 1, 1};
#endif
}

inline const char *amgcl_version() { return "1.4.5"; }

// Changes only the calling OpenMP task's settings and restores them on exit,
// including C++ exceptions. No R code is called while this guard is active.
class ThreadGuard {
public:
    explicit ThreadGuard(int requested) : effective_(1) {
        if (requested < 1)
            throw std::invalid_argument("threads must be a positive integer");
        effective_ = std::min(requested, thread_info().max_threads);
#ifdef _OPENMP
        previous_threads_ = omp_get_max_threads();
        previous_dynamic_ = omp_get_dynamic();
        omp_set_dynamic(0);
        omp_set_num_threads(effective_);
#endif
    }

    ~ThreadGuard() noexcept {
#ifdef _OPENMP
        omp_set_num_threads(previous_threads_);
        omp_set_dynamic(previous_dynamic_);
#endif
    }

    int threads() const { return effective_; }
    ThreadGuard(const ThreadGuard&) = delete;
    ThreadGuard& operator=(const ThreadGuard&) = delete;

private:
    int effective_;
#ifdef _OPENMP
    int previous_threads_;
    int previous_dynamic_;
#endif
};

// Owned, zero-based CSR storage for a nonempty square matrix. Within each row,
// column indices must be strictly increasing. A nonzero diagonal is required
// by the selected AMG hierarchy and its unpivoted coarse solver.
struct Matrix {
    std::size_t n = 0;
    std::vector<std::ptrdiff_t> ptr;
    std::vector<std::ptrdiff_t> col;
    std::vector<double> val;

    void validate() const {
        if (n == 0 || n > static_cast<std::size_t>(
                std::numeric_limits<std::ptrdiff_t>::max()))
            throw std::invalid_argument("A must be a nonempty square matrix");
        if (ptr.size() != n + 1 || ptr.front() != 0 || col.size() != val.size() ||
                val.size() > static_cast<std::size_t>(
                    std::numeric_limits<std::ptrdiff_t>::max()) ||
                ptr.back() != static_cast<std::ptrdiff_t>(val.size()))
            throw std::invalid_argument("Invalid sparse matrix storage");
        for (std::size_t row = 0; row < n; ++row) {
            if (ptr[row] < 0 || ptr[row] > ptr[row + 1] ||
                    ptr[row + 1] > static_cast<std::ptrdiff_t>(val.size()))
                throw std::invalid_argument("Invalid sparse matrix pointers");
            std::ptrdiff_t previous = -1;
            bool diagonal = false;
            for (std::ptrdiff_t k = ptr[row]; k < ptr[row + 1]; ++k) {
                if (col[k] <= previous || col[k] < 0 ||
                        static_cast<std::size_t>(col[k]) >= n)
                    throw std::invalid_argument("Invalid or unsorted sparse matrix indices");
                if (!std::isfinite(val[k]))
                    throw std::invalid_argument("A must contain only finite values");
                if (col[k] == static_cast<std::ptrdiff_t>(row) && val[k] != 0)
                    diagonal = true;
                previous = col[k];
            }
            if (!diagonal)
                throw std::invalid_argument("AMG requires a nonzero diagonal in every row");
        }
    }

    // The pointer lengths are explicit so malformed R slot storage can be
    // rejected before dereferencing it. The input arrays are never retained.
    static Matrix from_csc(std::size_t size, const int *p, std::size_t p_size,
                           const int *i, const double *x, std::size_t nnz) {
        if (size == 0 || size > static_cast<std::size_t>(
                std::numeric_limits<int>::max()) || p_size != size + 1 ||
                p == nullptr || (nnz && (i == nullptr || x == nullptr)) ||
                nnz > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::invalid_argument("Invalid dgCMatrix dimensions or slot lengths");
        if (p[0] != 0 || p[size] != static_cast<int>(nnz))
            throw std::invalid_argument("Invalid dgCMatrix column pointers");

        Matrix A;
        A.n = size;
        A.ptr.assign(size + 1, 0);
        for (std::size_t column = 0; column < size; ++column) {
            if (p[column] < 0 || p[column] > p[column + 1] ||
                    p[column + 1] > static_cast<int>(nnz))
                throw std::invalid_argument("Invalid dgCMatrix column pointers");
            int previous = -1;
            for (int k = p[column]; k < p[column + 1]; ++k) {
                if (i[k] <= previous || i[k] < 0 ||
                        static_cast<std::size_t>(i[k]) >= size)
                    throw std::invalid_argument("Invalid or unsorted dgCMatrix row indices");
                if (!std::isfinite(x[k]))
                    throw std::invalid_argument("A must contain only finite values");
                ++A.ptr[static_cast<std::size_t>(i[k]) + 1];
                previous = i[k];
            }
        }
        for (std::size_t row = 0; row < size; ++row)
            A.ptr[row + 1] += A.ptr[row];
        A.col.resize(nnz);
        A.val.resize(nnz);
        std::vector<std::ptrdiff_t> next(A.ptr);
        for (std::size_t column = 0; column < size; ++column) {
            for (int k = p[column]; k < p[column + 1]; ++k) {
                const std::ptrdiff_t dest = next[i[k]]++;
                A.col[dest] = static_cast<std::ptrdiff_t>(column);
                A.val[dest] = x[k];
            }
        }
        A.validate();
        return A;
    }
};

struct SolveResult {
    std::vector<double> x;
    std::size_t iterations = 0;
    double error = 0;
    double residual = 0;
    bool converged = false;
    int threads = 1;
};

class Solver {
    using Backend = amgcl::backend::builtin<double>;
    using Preconditioner = amgcl::amg<Backend,
        amgcl::coarsening::smoothed_aggregation, amgcl::relaxation::spai0>;

    struct Implementation {
        virtual ~Implementation() {}
        virtual std::tuple<std::size_t, double> solve(
            const std::vector<double>& rhs, std::vector<double>& x) = 0;
        virtual const Backend::matrix& system_matrix() const = 0;
    };

    template <class IterativeSolver>
    struct ImplementationFor : Implementation {
        using AMG = amgcl::make_solver<Preconditioner, IterativeSolver>;
        AMG amg;

        static typename AMG::params parameters(double tolerance, int maxiter) {
            typename AMG::params p;
            p.solver.tol = tolerance;
            p.solver.maxiter = static_cast<std::size_t>(maxiter);
            p.solver.verbose = false;
            p.precond.coarse_enough = 64;
            p.precond.allow_rebuild = false;
            // Equal pre/post SPAI0 relaxation preserves symmetry for CG.
            p.precond.npre = p.precond.npost = 1;
            return p;
        }

        ImplementationFor(const Matrix& A, double tol, int maxiter)
            : amg(std::tie(A.n, A.ptr, A.col, A.val), parameters(tol, maxiter)) {}

        std::tuple<std::size_t, double> solve(
                const std::vector<double>& rhs, std::vector<double>& x) override {
            return amg(rhs, x);
        }

        const Backend::matrix& system_matrix() const override {
            return amg.system_matrix();
        }
    };

public:
    explicit Solver(Matrix A, const std::string& method = "cg",
                    double tol = 1e-8, int maxiter = 1000, int threads = 1)
        : n_(A.n), tolerance_(tol) {
        A.validate();
        if (method != "cg" && method != "bicgstab")
            throw std::invalid_argument("solver must be 'cg' or 'bicgstab'");
        if (!std::isfinite(tol) || tol <= 0)
            throw std::invalid_argument("tol must be finite and positive");
        if (maxiter < 1)
            throw std::invalid_argument("maxiter must be a positive integer");
        ThreadGuard guard(threads);
        try {
            if (method == "cg")
                implementation_.reset(new ImplementationFor<amgcl::solver::cg<Backend>>(
                    A, tol, maxiter));
            else
                implementation_.reset(new ImplementationFor<amgcl::solver::bicgstab<Backend>>(
                    A, tol, maxiter));
        } catch (const std::runtime_error& error) {
            throw std::runtime_error(std::string("AMGCL setup failed; check that the matrix ") +
                "is nonsingular and suitable for the unpivoted coarse solver. AMGCL: " +
                error.what());
        }
    }

    Solver(const Solver&) = delete;
    Solver& operator=(const Solver&) = delete;

    std::size_t size() const { return n_; }

    // A solver reuses mutable AMG workspaces: do not call it concurrently.
    SolveResult solve(const std::vector<double>& rhs, int threads = 1) {
        if (rhs.size() != n_)
            throw std::invalid_argument("Length of b must equal nrow(A)");
        double scale = 0;
        for (double value : rhs) {
            if (!std::isfinite(value))
                throw std::invalid_argument("b must contain only finite values");
            scale = std::max(scale, std::abs(value));
        }
        SolveResult result;
        result.x.assign(n_, 0.0);
        ThreadGuard guard(threads);
        result.threads = guard.threads();
        if (scale != 0) {
            // AMGCL treats sufficiently small RHS norms as zero. Scaling also
            // avoids overflow when computing the RHS inner product.
            std::vector<double> normalized(rhs.size());
            for (std::size_t j = 0; j < rhs.size(); ++j)
                normalized[j] = rhs[j] / scale;
            std::tie(result.iterations, result.error) =
                implementation_->solve(normalized, result.x);
            for (double& value : result.x) {
                value *= scale;
                if (!std::isfinite(value))
                    throw std::runtime_error("AMGCL returned a non-finite solution; check the matrix");
            }
            if (!std::isfinite(result.error))
                throw std::runtime_error("AMGCL returned a non-finite error estimate; check the matrix");
        }
        result.residual = residual(rhs, result.x);
        result.converged = result.residual <= tolerance_;
        return result;
    }

private:
    std::size_t n_;
    double tolerance_;
    std::unique_ptr<Implementation> implementation_;

    double residual(const std::vector<double>& rhs, const std::vector<double>& x) const {
        // AMG owns the original system matrix at the finest level, including
        // when the hierarchy consists only of its direct coarse solver.
        const Backend::matrix& A = implementation_->system_matrix();
        using Norms = std::pair<long double, long double>;
        const auto accumulate = [&](std::size_t begin, std::size_t end) {
            Norms norms(0, 0);
            for (std::size_t row = begin; row < end; ++row) {
                long double ax = 0;
                for (std::ptrdiff_t k = A.ptr[row]; k < A.ptr[row + 1]; ++k)
                    ax += static_cast<long double>(A.val[k]) * x[A.col[k]];
                norms.first = std::hypot(norms.first, static_cast<long double>(rhs[row]) - ax);
                norms.second = std::hypot(norms.second, static_cast<long double>(rhs[row]));
            }
            return norms;
        };
        Norms norms(0, 0);
        if (n_ < 32768) {
            norms = accumulate(0, n_);
        } else {
            // Fixed blocks give the same reduction order for any thread count.
            // Hypot within and between blocks avoids overflow/underflow from
            // squaring the residual or RHS, including on platforms where
            // long double has the same range as double. Workers use no R API.
            const std::size_t block_size = 4096;
            const std::ptrdiff_t blocks = static_cast<std::ptrdiff_t>(
                n_ / block_size + (n_ % block_size != 0));
            std::vector<Norms> partials(static_cast<std::size_t>(blocks));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(omp_get_max_threads() > 1)
#endif
            for (std::ptrdiff_t block = 0; block < blocks; ++block) {
                const std::size_t begin = static_cast<std::size_t>(block) * block_size;
                partials[block] = accumulate(begin, std::min(n_, begin + block_size));
            }
            for (const Norms& partial : partials) {
                norms.first = std::hypot(norms.first, partial.first);
                norms.second = std::hypot(norms.second, partial.second);
            }
        }
        const double value = static_cast<double>(norms.second == 0 ? norms.first : norms.first / norms.second);
        if (!std::isfinite(value))
            throw std::runtime_error("The computed residual is non-finite; check the matrix");
        return value;
    }
};

} // namespace rcppamgcl

#endif
