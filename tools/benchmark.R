#!/usr/bin/env Rscript
# Outside tests/examples: benchmarks deliberately request up to eight threads.
# Usage: Rscript tools/benchmark.R [output.csv] [repetitions=3] [minimum_n=160000]
library(RcppAMGCL)
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args) >= 1L) args[[1L]] else "benchmarks/results-local.csv"
repetitions <- if (length(args) >= 2L) as.integer(args[[2L]]) else 3L
minimum_n <- if (length(args) >= 3L) as.integer(args[[3L]]) else 160000L
stopifnot(length(repetitions) == 1L, !is.na(repetitions), repetitions >= 1L,
          length(minimum_n) == 1L, !is.na(minimum_n), minimum_n >= 100000L)

poisson <- function(dimension) {
    side <- ceiling(minimum_n^(1 / dimension))
    n <- as.integer(side^dimension)
    nodes <- seq_len(n)
    source <- target <- vector("list", dimension)
    for (axis in seq_len(dimension)) {
        stride <- side^(axis - 1L)
        source[[axis]] <- nodes[((nodes - 1L) %/% stride) %% side < side - 1L]
        target[[axis]] <- source[[axis]] + stride
    }
    source <- unlist(source, use.names = FALSE)
    target <- unlist(target, use.names = FALSE)
    Matrix::sparseMatrix(i = c(nodes, source, target), j = c(nodes, target, source),
                         x = c(rep(2 * dimension, n), rep(-1, 2L * length(source))),
                         dims = c(n, n))
}

info <- amgcl_info()
print(info)
if (!info$openmp) message("Serial build: all requested thread counts will use one thread.")
rows <- list()
for (dimension in c(2L, 3L)) {
    A <- poisson(dimension)
    n <- nrow(A)
    reference <- sin(seq_len(n) / 17) + cos(seq_len(n) / 113)
    b <- as.numeric(A %*% reference)
    for (threads in c(1L, 2L, 4L, 8L)) {
        # Warm up runtime/allocator outside the measured repetitions.
        warm <- amgcl_setup(A, threads = threads)
        invisible(amgcl_solve(warm, b, threads = threads))
        rm(warm)
        for (repetition in seq_len(repetitions)) {
            invisible(gc())
            setup_seconds <- unname(system.time({
                solver <- amgcl_setup(A, solver = "cg", tol = 1e-8, threads = threads)
            })[["elapsed"]])
            solve_seconds <- unname(system.time({
                result <- amgcl_solve(solver, b, threads = threads)
            })[["elapsed"]])
            stopifnot(result$converged, is.finite(result$residual), result$residual < 1e-7)
            row <- data.frame(matrix = paste0("poisson", dimension, "d"), n = n,
                              nnz = length(A@x), solver = "cg",
                              requested_threads = threads, effective_threads = result$threads,
                              repetition = repetition, setup_seconds = setup_seconds,
                              solve_seconds = solve_seconds, iterations = result$iterations,
                              residual = result$residual, converged = result$converged,
                              openmp = info$openmp, max_threads = info$max_threads)
            rows[[length(rows) + 1L]] <- row
            print(row, row.names = FALSE)
            rm(solver)
        }
    }
    rm(A)
}
results <- do.call(rbind, rows)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(results, output, row.names = FALSE)
writeLines(c(capture.output(sessionInfo()), "", capture.output(Sys.info()), "",
             capture.output(info), "",
             paste("OMP_THREAD_LIMIT:", Sys.getenv("OMP_THREAD_LIMIT", unset = "unset")),
             paste("OMP_NUM_THREADS:", Sys.getenv("OMP_NUM_THREADS", unset = "unset"))),
           paste0(output, ".session.txt"))
message("Wrote measurements to ", output)
