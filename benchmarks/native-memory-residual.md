# Native memory and residual measurements

These Windows measurements use Rtools GCC 14.3.0, C++14, `-O2`, and OpenMP.
They compare the original implementation (`6cb9d30`), removal of the retained
CSR copy (`8b96cca`), and the parallel residual implementation (`771186f`).
The exact original header SHA256 and benchmark source SHA256 accompany each
CSV file. The latter two were measured from working-tree snapshots, then
verified to match the committed header blobs. The machine reported capacity
for 16 OpenMP threads.

The workload is a 100 by 100 by 100 three-dimensional Poisson matrix:
1,000,000 rows and 6,940,000 nonzeros. Each requested thread count (1, 2, 4, 8)
runs in a fresh process, constructs one reusable CG solver, warms it with one
solve, and measures three more solves. Every solve took 13 iterations and
achieved an independently verified relative residual of approximately
6.23e-9 against a tolerance of 1e-8. Verification applies the known Poisson
stencil independently of the stored sparse matrix and is outside the timer.

## Persistent memory

The standalone benchmark overrides C++ `new`, `new[]`, and corresponding
deallocation operators to count requested live allocation payloads. It moves
the input CSR into a stack-allocated solver, then records the live allocation
increase after setup and before allocating the RHS, reference, or solution.
Using a stack-allocated outer solver excludes its changed object size from
this comparison. All four thread counts gave the same figures.

| Implementation | Live C++ solver allocation payload |
| --- | ---: |
| Original | 519,644,092 bytes |
| Reuse AMG's matrix | 400,604,084 bytes |
| Reuse matrix and parallel residual | 400,604,084 bytes |

The measured saving is **119,040,008 bytes (113.525 MiB)**, exactly the removed
64-bit CSR arrays: `8 * (n + 1) + 8 * nnz + 8 * nnz`. This is 22.9% of the
original tracked persistent payload. These numbers exclude R memory, direct
`malloc` allocations, OpenMP runtime memory, and allocator bookkeeping. They
measure persistent C++ allocations, not resident process memory or peak setup
memory; the input matrix and AMG's copy can both exist during setup.

## Residual cost

The helper instruments only a temporary header copy, placing a steady-clock
timer immediately around the actual residual calculation. Production headers
contain no timing hooks. The total solve timer includes normalization, solver
execution, rescaling, and the independent residual calculation.

The table reports medians of three measurements. Residual shares are medians
of the paired residual-time/solve-time ratios from individual repetitions.

| Threads | Original residual (s) | Original share | Matrix reuse residual (s) | Matrix reuse share | Parallel residual (s) | Parallel share |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.13686 | 14.90% | 0.14287 | 15.09% | 0.14837 | 15.72% |
| 2 | 0.13842 | 22.39% | 0.13742 | 22.56% | 0.06868 | 12.27% |
| 4 | 0.13480 | 29.87% | 0.13718 | 31.14% | 0.03431 | 10.19% |
| 8 | 0.13508 | 35.68% | 0.13584 | 36.04% | 0.02527 | 10.14% |

| Threads | Original total solve (s) | Matrix reuse total solve (s) | Parallel residual total solve (s) |
| ---: | ---: | ---: | ---: |
| 1 | 0.92023 | 0.94663 | 0.94374 |
| 2 | 0.66174 | 0.61085 | 0.56379 |
| 4 | 0.45018 | 0.44377 | 0.33758 |
| 8 | 0.38468 | 0.38035 | 0.24928 |

At eight threads, parallelizing the residual reduced its median time by about
5.38-fold and the median total native solve time by 34.5% relative to matrix
reuse alone. The parallel implementation retains long-double accumulation
and `hypot`, using fixed 4,096-row blocks and a fixed final reduction order.
One-thread differences are small relative to the observed timing variability.
These instrumented measurements include allocation-tracker overhead and were
collected sequentially on one machine; they are diagnostic results, not a
general performance guarantee. The existing R benchmarks measure the package
interface separately.

## Reproduction and evidence

Run from the repository root with an OpenMP-capable compiler:

```powershell
./tools/benchmark-native.ps1 -Revision 6cb9d30 -Label baseline
./tools/benchmark-native.ps1 -Revision 8b96cca -Label memory
./tools/benchmark-native.ps1 -Revision 771186f -Label parallel
```

The default side length is 100 and the default repetition count is three.
`-BuildOnly` and `-RunOnly` separate compilation from timed runs. No other
package builds or benchmarks ran concurrently with these measurements.

Raw measurements and source/environment metadata:

- [Original CSV](native-baseline.csv) and [metadata](native-baseline.metadata.txt).
- [Matrix reuse CSV](native-memory.csv) and [metadata](native-memory.metadata.txt).
- [Parallel residual CSV](native-parallel.csv) and [metadata](native-parallel.metadata.txt).

The helper now records run-time OpenMP environment settings even with
`-RunOnly`. The first two measurement stages predated that extra metadata;
their build-time `OMP_THREAD_LIMIT` and `OMP_NUM_THREADS` values were recorded
and remained unset during measurement. Their other run-time settings and
start timestamps were not recorded.
