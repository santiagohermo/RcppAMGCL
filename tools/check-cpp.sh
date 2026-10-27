#!/usr/bin/env bash
# Build a true standalone executable: no R headers or R runtime linkage.
# The CI compilers are GCC/Rtools and Apple clang with Homebrew libomp.
set -euo pipefail
mode="${1:-${RCPPAMGCL_OPENMP:-off}}"
case "$mode" in on|off) ;; *) echo "Expected on or off" >&2; exit 1 ;; esac
read -r -a compiler <<< "$(R CMD config CXX)"
# R >= 4.6 no longer exposes CXX14 variables. Select the standard directly for
# this GCC/Clang-only standalone compatibility test, using R's current compiler.
standard=(-std=c++14)
flags=()
if [[ "$mode" == on ]]; then
    if [[ "$(uname -s)" == Darwin ]]; then
        omp_prefix="$(brew --prefix libomp)"
        flags=(-Xclang -fopenmp "-I$omp_prefix/include" "-L$omp_prefix/lib"
               "-Wl,-rpath,$omp_prefix/lib" -lomp)
    else
        flags=(-fopenmp)
    fi
else
    flags=(-DRCPPAMGCL_DISABLE_OPENMP)
fi
mkdir -p .local
"${compiler[@]}" "${standard[@]}" -O2 -Iinst/include tools/smoke.cpp \
    "${flags[@]}" -o .local/cpp-smoke
OMP_THREAD_LIMIT=2 .local/cpp-smoke "$mode"
OMP_THREAD_LIMIT=1 .local/cpp-smoke "$mode"
