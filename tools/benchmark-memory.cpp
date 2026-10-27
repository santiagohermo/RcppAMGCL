// Private benchmark helper, excluded from the R package build with tools/.
// Compile a copy in a temporary directory with: R CMD SHLIB benchmark-memory.cpp
// After dyn.load(), .Call("RcppAMGCL_benchmark_peak_rss") returns process peak
// resident memory in BYTES, or NA_real_ if the platform/API is unavailable.
// .Call("RcppAMGCL_benchmark_time") returns monotonic high-resolution seconds;
// subtract two readings to measure elapsed time independently of wall-clock time.
// Call both exports on R's main thread. Peak memory is a process-lifetime
// high-water mark including R, libraries, inputs, and setup temporaries;
// it cannot be reset.

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
// Version 2 resolves GetProcessMemoryInfo to the Kernel32 export, avoiding
// an additional psapi import-library link requirement on modern Windows.
#ifndef PSAPI_VERSION
#define PSAPI_VERSION 2
#endif
#include <windows.h>
#include <psapi.h>
#elif defined(__linux__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

#include <chrono>

#include <R.h>
#include <Rinternals.h>

extern "C" SEXP RcppAMGCL_benchmark_peak_rss() {
    double bytes = NA_REAL;
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS counters = {};
    counters.cb = sizeof(counters);
    if (GetProcessMemoryInfo(GetCurrentProcess(), &counters, sizeof(counters)))
        bytes = static_cast<double>(counters.PeakWorkingSetSize);
#elif defined(__linux__) || defined(__APPLE__)
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0 && usage.ru_maxrss >= 0) {
        bytes = static_cast<double>(usage.ru_maxrss);
#ifdef __linux__
        bytes *= 1024.0;
#endif
    }
#endif
    return Rf_ScalarReal(bytes);
}

extern "C" SEXP RcppAMGCL_benchmark_time() {
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
    return Rf_ScalarReal(seconds);
}
