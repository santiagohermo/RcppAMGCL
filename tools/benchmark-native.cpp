// Diagnostic benchmark, deliberately outside CRAN tests and examples.
// Build with benchmark-native.ps1, which instruments a temporary header copy.
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <new>

namespace allocation_measurement {
// Counts requested live C++ allocation payloads, not resident memory, allocator
// bookkeeping, malloc allocations, OpenMP runtime memory, or setup peak memory.
std::atomic<std::size_t> live_bytes{0};
struct alignas(std::max_align_t) Header { std::size_t size; };

void *allocate(std::size_t size) {
    if (size > std::numeric_limits<std::size_t>::max() - sizeof(Header))
        throw std::bad_alloc();
    const std::size_t payload = size == 0 ? 1 : size;
    void *storage = std::malloc(sizeof(Header) + payload);
    if (!storage) throw std::bad_alloc();
    Header *header = static_cast<Header *>(storage);
    header->size = size;
    live_bytes.fetch_add(size, std::memory_order_relaxed);
    return header + 1;
}

void release(void *pointer) noexcept {
    if (!pointer) return;
    Header *header = static_cast<Header *>(pointer) - 1;
    live_bytes.fetch_sub(header->size, std::memory_order_relaxed);
    std::free(header);
}
}

void *operator new(std::size_t size) { return allocation_measurement::allocate(size); }
void *operator new[](std::size_t size) { return allocation_measurement::allocate(size); }
void operator delete(void *pointer) noexcept { allocation_measurement::release(pointer); }
void operator delete[](void *pointer) noexcept { allocation_measurement::release(pointer); }
void operator delete(void *pointer, std::size_t) noexcept { allocation_measurement::release(pointer); }
void operator delete[](void *pointer, std::size_t) noexcept { allocation_measurement::release(pointer); }

// Only the temporary, instrumented header refers to this benchmark-owned value.
namespace rcppamgcl_benchmark { double residual_seconds = -1; }
#include <RcppAMGCL.h>

namespace {
using Clock = std::chrono::steady_clock;

double elapsed(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

rcppamgcl::Matrix poisson3d(std::size_t side) {
    const std::size_t plane = side * side;
    const std::size_t n = plane * side;
    const std::size_t nnz = 7 * n - 6 * plane;
    rcppamgcl::Matrix A;
    A.n = n;
    A.ptr.resize(n + 1);
    A.col.resize(nnz);
    A.val.resize(nnz);
    std::size_t position = 0;
    auto append = [&](std::size_t column, double value) {
        A.col[position] = static_cast<std::ptrdiff_t>(column);
        A.val[position] = value;
        ++position;
    };
    for (std::size_t row = 0; row < n; ++row) {
        A.ptr[row] = static_cast<std::ptrdiff_t>(position);
        const std::size_t x = row % side;
        const std::size_t y = (row / side) % side;
        const std::size_t z = row / plane;
        if (z) append(row - plane, -1);
        if (y) append(row - side, -1);
        if (x) append(row - 1, -1);
        append(row, 6);
        if (x + 1 < side) append(row + 1, -1);
        if (y + 1 < side) append(row + side, -1);
        if (z + 1 < side) append(row + plane, -1);
    }
    A.ptr[n] = static_cast<std::ptrdiff_t>(position);
    if (position != nnz) throw std::runtime_error("Incorrect Poisson matrix nnz");
    return A;
}

// Independent verification directly applies the known 3D stencil; it does not
// access solver internals or retain a second sparse matrix.
long double stencil_value(const std::vector<double>& x, std::size_t row,
                          std::size_t side) {
    const std::size_t plane = side * side;
    long double value = 6.0L * x[row];
    if (row / plane) value -= x[row - plane];
    if ((row / side) % side) value -= x[row - side];
    if (row % side) value -= x[row - 1];
    if (row % side + 1 < side) value -= x[row + 1];
    if ((row / side) % side + 1 < side) value -= x[row + side];
    if (row / plane + 1 < side) value -= x[row + plane];
    return value;
}

double verify(const std::vector<double>& x, const std::vector<double>& rhs,
              std::size_t side) {
    long double squared_error = 0;
    long double squared_rhs = 0;
    for (std::size_t row = 0; row < x.size(); ++row) {
        const long double difference = static_cast<long double>(rhs[row]) -
            stencil_value(x, row, side);
        squared_error += difference * difference;
        squared_rhs += static_cast<long double>(rhs[row]) * rhs[row];
    }
    return static_cast<double>(std::sqrt(squared_error / squared_rhs));
}
}

int main(int argc, char **argv) {
    try {
        if (argc != 4) throw std::runtime_error("Usage: benchmark-native SIDE THREADS REPETITIONS");
        const int side_arg = std::atoi(argv[1]);
        const int threads = std::atoi(argv[2]);
        const int repetitions = std::atoi(argv[3]);
        if (side_arg < 4 || side_arg > 1000 || threads < 1 || repetitions < 1)
            throw std::runtime_error("Invalid benchmark arguments");
        const std::size_t side = static_cast<std::size_t>(side_arg);
        const std::size_t n = side * side * side;
        const std::size_t nnz = 7 * n - 6 * side * side;
        const auto info = rcppamgcl::thread_info();
        const std::size_t before = allocation_measurement::live_bytes.load();
        rcppamgcl::Matrix A = poisson3d(side);
        const std::size_t input_csr_bytes = allocation_measurement::live_bytes.load() - before;
        const auto setup_begin = Clock::now();
        // Stack allocation excludes the changed size of the outer Solver object.
        // Move the input so the snapshot measures persistent solver allocations.
        rcppamgcl::Solver solver(std::move(A), "cg", 1e-8, 1000, threads);
        const double setup_seconds = elapsed(setup_begin);
        const std::size_t solver_live_bytes = allocation_measurement::live_bytes.load() - before;

        std::vector<double> reference(n);
        for (std::size_t row = 0; row < n; ++row)
            reference[row] = std::sin(static_cast<double>(row + 1) / 17) +
                std::cos(static_cast<double>(row + 1) / 113);
        std::vector<double> rhs(n);
        for (std::size_t row = 0; row < n; ++row)
            rhs[row] = static_cast<double>(stencil_value(reference, row, side));
        {
            const auto warm = solver.solve(rhs, threads);
            if (!warm.converged) throw std::runtime_error("Warmup did not converge");
        }
        std::cout << "n,nnz,requested_threads,effective_threads,repetition,input_csr_bytes,"
                     "solver_live_cpp_bytes,setup_seconds,solve_seconds,residual_seconds,"
                     "residual_fraction,iterations,residual,verified_residual,openmp,max_threads\n";
        std::cout << std::setprecision(12);
        for (int repetition = 1; repetition <= repetitions; ++repetition) {
            const auto solve_begin = Clock::now();
            const auto result = solver.solve(rhs, threads);
            const double solve_seconds = elapsed(solve_begin);
            const double residual_seconds = rcppamgcl_benchmark::residual_seconds;
            const double independently_verified = verify(result.x, rhs, side);
            if (residual_seconds < 0)
                throw std::runtime_error("Temporary header timing instrumentation is missing");
            if (!result.converged || !std::isfinite(independently_verified) ||
                    independently_verified > 1e-8)
                throw std::runtime_error("Measured solve failed independent residual verification");
            std::cout << n << ',' << nnz << ',' << threads << ',' << result.threads << ','
                      << repetition << ',' << input_csr_bytes << ',' << solver_live_bytes << ','
                      << setup_seconds << ',' << solve_seconds << ',' << residual_seconds << ','
                      << residual_seconds / solve_seconds << ',' << result.iterations << ','
                      << result.residual << ',' << independently_verified << ','
                      << info.openmp << ',' << info.max_threads << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
