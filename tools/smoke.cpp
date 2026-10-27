// Standalone C++ integration check: no R headers or R runtime are required.
// Usage: smoke [on|off]. Built and run in both modes by GitHub Actions.
#include <RcppAMGCL.h>
#include <iostream>

namespace {
void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main(int argc, char **argv) {
    try {
#ifdef _OPENMP
        omp_set_dynamic(1);
        omp_set_num_threads(2);
#endif
        const auto before = rcppamgcl::thread_info();
        if (argc > 1) {
            const std::string mode(argv[1]);
            require(mode == "on" || mode == "off", "Expected mode must be on or off");
            require(before.openmp == (mode == "on"), "Unexpected OpenMP build mode");
        }
        int actual_threads = 1;
#ifdef _OPENMP
        {
            rcppamgcl::ThreadGuard guard(2);
            actual_threads = 0;
#pragma omp parallel reduction(+:actual_threads)
            { actual_threads += 1; }
            require(actual_threads == guard.threads(), "Unexpected OpenMP team size");
        }
#endif
        for (const std::string method : {"cg", "bicgstab"}) {
            rcppamgcl::Matrix A;
            A.n = 100;
            A.ptr.push_back(0);
            for (std::ptrdiff_t row = 0; row < 100; ++row) {
                if (row) { A.col.push_back(row - 1); A.val.push_back(-1); }
                A.col.push_back(row); A.val.push_back(3);
                if (row < 99) { A.col.push_back(row + 1); A.val.push_back(-1); }
                A.ptr.push_back(A.col.size());
            }
            rcppamgcl::Solver solver(A, method, 1e-10, 1000, 1);
            require(rcppamgcl::thread_info().current_threads == before.current_threads,
                    "Setup did not restore the previous OpenMP setting");
            for (double scale : {1.0, 1e-200, 1e200}) {
                for (int requested : {1, 2}) {
                    auto result = solver.solve(std::vector<double>(100, scale), requested);
                    require(result.converged && result.residual < 1e-10,
                            "Scaled RHS did not converge");
                    require(result.threads == std::min(requested, before.max_threads),
                            "Requested thread count was not capped correctly");
                    require(rcppamgcl::thread_info().current_threads == before.current_threads,
                            "Solve did not restore the previous OpenMP setting");
                }
            }
            auto zero = solver.solve(std::vector<double>(100, 0), 2);
            require(zero.converged && zero.residual == 0 && zero.iterations == 0,
                    "Zero RHS did not return the trivial solution");
            bool caught = false;
            try { solver.solve(std::vector<double>(1, 1), 2); }
            catch (const std::invalid_argument&) { caught = true; }
            require(caught, "Invalid RHS length was accepted");
        }

        // Nonsymmetric CSC detects accidental transposition in the adapter.
        const int csc_p[] = {0, 2, 5, 7};
        const int csc_i[] = {0, 1, 0, 1, 2, 1, 2};
        const double csc_x[] = {4, -2, -1, 4, -2, -1, 4};
        auto nonsym = rcppamgcl::Matrix::from_csc(3, csc_p, 4, csc_i, csc_x, 7);
        rcppamgcl::Solver nonsym_solver(nonsym, "bicgstab", 1e-12, 1000, 2);
        auto result = nonsym_solver.solve({2, 3, 8}, 2);
        require(result.converged && result.residual < 1e-12,
                "Nonsymmetric CSC solve did not converge");
        for (int j = 0; j < 3; ++j)
            require(std::abs(result.x[j] - (j + 1)) < 1e-12,
                    "Nonsymmetric CSC solve returned the wrong solution");

        // This fails inside AMG setup, after the temporary thread settings
        // have been installed, so it exercises exception unwinding.
        rcppamgcl::Matrix singular;
        singular.n = 2;
        singular.ptr = {0, 2, 4};
        singular.col = {0, 1, 0, 1};
        singular.val = {1, 1, 1, 1};
        bool caught = false;
        try { rcppamgcl::Solver bad(singular, "cg", 1e-8, 1000, 1); }
        catch (const std::runtime_error&) { caught = true; }
        require(caught, "Singular coarse factorization unexpectedly succeeded");
        require(rcppamgcl::thread_info().current_threads == before.current_threads,
                "Setup exception did not restore the previous OpenMP setting");
#ifdef _OPENMP
        require(omp_get_dynamic() == 1, "Dynamic OpenMP setting was not restored");
#endif
        std::cout << "C++ header checks passed; OpenMP=" << before.openmp
                  << "; actual team=" << actual_threads << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
