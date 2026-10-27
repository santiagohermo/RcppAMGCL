# Native allocation and residual diagnostics

On the measured 1,000,000-row, 6,940,000-nonzero three-dimensional Poisson
problem, the solver retains **400,604,084 bytes of tracked C++ allocations**.
At eight threads, its median residual calculation takes **0.02527 seconds**
within a **0.24928-second solve**. These diagnostics complement the R-facing
setup, solve, and total-time benchmarks in the [benchmark overview](README.md).

## Measurement scope

Measurements use Windows, Rtools GCC 14.3.0, C++14, `-O2`, and OpenMP on a
machine reporting capacity for 16 threads. Each requested thread count runs
in a fresh process, constructs one reusable CG solver, performs one warmup
solve, and measures three more solves. Every solve takes 13 iterations and
achieves an independently verified relative residual of approximately 6.23e-9
against a tolerance of 1e-8. Verification applies the known Poisson stencil
outside the solve timer.

The allocation tracker counts requested live C++ `new`/`new[]` payloads after
setup, with the input CSR moved into a stack-allocated solver and before
allocating the RHS, reference, or solution. It excludes R memory, direct
`malloc` allocations, OpenMP runtime memory, allocator bookkeeping, and the
outer solver object's size. It measures persistent allocation payload,
**not peak resident memory or peak setup memory**.

The allocation audit gives 519,644,092 bytes when a separate input CSR copy
is retained and 400,604,084 bytes when the solver uses AMG's matrix. Their
**119,040,008-byte (113.525 MiB)** difference exactly matches the 64-bit CSR
arrays: `8 * (n + 1) + 8 * nnz + 8 * nnz`. All four thread counts agree.

## Residual and solve time

The helper inserts a timer around the actual residual call in a temporary
header copy; production headers contain no timing hooks. The solve timer
includes normalization, iteration, rescaling, and residual calculation.
The residual uses long-double accumulation and `hypot`, with fixed 4,096-row
blocks and a fixed final reduction order.

| Threads | Residual (s) | Solve only (s) | Residual share of solve |
| ---: | ---: | ---: | ---: |
| 1 | 0.14837 | 0.94374 | 15.72% |
| 2 | 0.06868 | 0.56379 | 12.27% |
| 4 | 0.03431 | 0.33758 | 10.19% |
| 8 | 0.02527 | 0.24928 | 10.14% |

Times are medians of three measurements; shares are medians of the paired
residual-time/solve-time ratios. Raw `total_seconds` values add the measured
setup time to each solve time. Each configuration has one setup observation
and three warmed solve observations, so these derived totals are **not three
fresh one-off measurements**. Setup is excluded from the table above.

The reference audit with a serial residual reports 0.13584 seconds for the
residual and 0.38035 seconds for the solve at eight threads. The corresponding
parallel figures are about 5.38 times faster for the residual and 34.5% lower
for the solve. These are instrumented measurements from one machine, including
allocation-tracker overhead; small timing differences should not be treated
as general performance guarantees. Builds and other benchmarks were kept out
of the measurement periods.

## Reproduction and audit data

Run from the repository root with an OpenMP-capable compiler:

```powershell
./tools/benchmark-native.ps1 -Revision 6cb9d30 -Label baseline
./tools/benchmark-native.ps1 -Revision 8b96cca -Label memory
./tools/benchmark-native.ps1 -Revision 771186f -Label parallel
```

The defaults are side length 100 and three repetitions. `-BuildOnly` and
`-RunOnly` separate compilation from measurement. Exact header and benchmark
source hashes accompany the raw data; snapshots were verified against the
named commits.

- [Retained CSR audit](native-baseline.csv), source `6cb9d30`, [metadata](native-baseline.metadata.txt).
- [Shared matrix audit](native-memory.csv), source `8b96cca`, [metadata](native-memory.metadata.txt).
- [Parallel residual measurements](native-parallel.csv), source `771186f`, [metadata](native-parallel.metadata.txt).

The helper records runtime OpenMP settings with `-RunOnly`. The first two
datasets recorded `OMP_THREAD_LIMIT` and `OMP_NUM_THREADS`, both unset during
measurement; their additional runtime settings and start timestamps were
not recorded.
