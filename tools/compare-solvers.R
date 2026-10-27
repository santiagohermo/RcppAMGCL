#!/usr/bin/env Rscript
# Repository-only benchmark; deliberately permits more than CRAN's two threads.
# Rscript tools/compare-solvers.R output.csv [sizes=10000,40000,100000]
#     [repetitions=3] [right_hand_sides=10] [benchmark-memory DLL path]
# Compile tools/benchmark-memory.cpp once with R CMD SHLIB. The DLL may also
# be specified with RCPPAMGCL_BENCHMARK_MEMORY_DLL. No timing includes that build.
# Each configuration runs in a fresh process, with serial BLAS settings, one
# warmup setup/solve, and a 300-second timeout. Every setup uses a fresh dsCMatrix
# with an empty factor cache. Construction and external validation are untimed.
# AMGCL's internal residual calculation remains part of its public solve time.
# RSS is the process lifetime high-water mark, including inputs, warmup, R,
# factors, workspaces, results, and validation; it is not factor-only memory.

positive_integer <- function(value, name) {
    number <- suppressWarnings(as.numeric(value))
    if (length(number) != 1L || is.na(number) || !is.finite(number) ||
        number < 1 || number > .Machine$integer.max || number != floor(number)) {
        stop(name, " must be a positive integer", call. = FALSE)
    }
    as.integer(number)
}

dimensions <- function(minimum_n, dimension) {
    side <- as.integer(ceiling(minimum_n^(1 / dimension)))
    # Protect exact cubes against floating-point root rounding.
    while (side > 1L && (side - 1)^dimension >= minimum_n) side <- side - 1L
    while (side^dimension < minimum_n) side <- side + 1L
    n <- side^dimension
    if (n > .Machine$integer.max) stop("Requested matrix is too large")
    list(side = side, n = as.integer(n))
}

poisson <- function(minimum_n, dimension) {
    shape <- dimensions(minimum_n, dimension)
    nodes <- seq_len(shape$n)
    source <- target <- vector("list", dimension)
    for (axis in seq_len(dimension)) {
        stride <- shape$side^(axis - 1L)
        source[[axis]] <- nodes[((nodes - 1L) %/% stride) %% shape$side < shape$side - 1L]
        target[[axis]] <- source[[axis]] + stride
    }
    source <- unlist(source, use.names = FALSE)
    target <- unlist(target, use.names = FALSE)
    Matrix::sparseMatrix(i = c(nodes, source), j = c(nodes, target),
                         x = c(rep(2 * dimension, shape$n), rep(-1, length(source))),
                         dims = c(shape$n, shape$n), symmetric = TRUE)
}

fresh_matrix <- function(A) {
    # Construct a distinct S4 object and a fresh factor list for every setup.
    # Sharing immutable sparse coefficient slots avoids an unrelated data copy.
    result <- methods::new("dsCMatrix", Dim = A@Dim, Dimnames = A@Dimnames,
                           p = A@p, i = A@i, x = A@x, uplo = A@uplo,
                           factors = list())
    stopifnot(length(result@factors) == 0L)
    result
}

new_row <- function(config, repetition) {
    shape <- dimensions(config$minimum_n, config$dimension)
    data.frame(case_id = config$id, matrix = paste0("poisson", config$dimension, "d"),
               minimum_n = config$minimum_n, side = shape$side, n = shape$n,
               matrix_class = "dsCMatrix", nnz = NA_real_,
               method = config$method, requested_threads = config$threads,
               effective_threads = NA_integer_, repetition = repetition,
               rhs_count = config$rhs_count, completed_rhs = 0L,
               factor_class = NA_character_, setup_seconds = NA_real_,
               first_solve_seconds = NA_real_, total_seconds = NA_real_,
               solve_all_seconds = NA_real_, batch_total_seconds = NA_real_,
               max_residual = NA_real_, max_solution_error = NA_real_,
               max_reported_residual = NA_real_, max_iterations = NA_integer_,
               checkpoint_peak_rss_mib = NA_real_, peak_rss_mib = NA_real_,
               peak_rss_complete = FALSE, openmp = NA, max_threads = NA_integer_,
               thread_limit = NA_integer_, child_exit_status = NA_integer_,
               status = "not_run", message = "", stringsAsFactors = FALSE)
}

append_message <- function(old, new) {
    if (nzchar(old)) paste(old, new, sep = "; ") else new
}

validate <- function(A, rhs, solution, reference) {
    if (length(solution) != length(rhs) || any(!is.finite(solution))) {
        return(c(residual = Inf, solution_error = Inf))
    }
    rhs_scale <- max(abs(rhs))
    if (rhs_scale == 0) rhs_scale <- 1
    normalized_rhs <- rhs / rhs_scale
    difference <- as.double(A %*% (solution / rhs_scale)) - normalized_rhs
    denominator <- sqrt(sum(normalized_rhs^2))
    residual <- sqrt(sum(difference^2))
    if (denominator > 0) residual <- residual / denominator
    reference_scale <- max(abs(reference))
    if (reference_scale == 0) reference_scale <- 1
    normalized_reference <- reference / reference_scale
    solution_error <- sqrt(sum((solution / reference_scale - normalized_reference)^2))
    reference_norm <- sqrt(sum(normalized_reference^2))
    if (reference_norm > 0) solution_error <- solution_error / reference_norm
    c(residual = residual, solution_error = solution_error)
}

check_accuracy <- function(metrics) {
    all(is.finite(metrics)) && metrics[["residual"]] <= 1e-8 &&
        metrics[["solution_error"]] <= 1e-5
}

compiler_info <- function() {
    executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "R.exe" else "R")
    unlist(lapply(c("CXX", "CXXFLAGS"), function(variable) {
        result <- tryCatch(suppressWarnings(system2(executable,
            c("CMD", "config", variable), stdout = TRUE, stderr = TRUE)),
            error = function(error) conditionMessage(error))
        paste(variable, paste(result, collapse = " "), sep = " = ")
    }))
}

write_metadata <- function(path, config, info) {
    header <- system.file("include", "RcppAMGCL.h", package = "RcppAMGCL")
    helper_source <- file.path(dirname(config$script), "benchmark-memory.cpp")
    hashes <- tools::md5sum(c(script = config$script, wrapper_header = header,
                             helper_source = helper_source, memory_dll = config$memory_dll))
    source_commit <- tryCatch(suppressWarnings(system2("git",
        c("-C", shQuote(dirname(dirname(config$script))), "rev-parse", "HEAD"),
        stdout = TRUE, stderr = TRUE)), error = function(error) conditionMessage(error))
    session <- sessionInfo()
    blas <- c(session$BLAS, extSoftVersion()[intersect(c("BLAS", "LAPACK"), names(extSoftVersion()))])
    blas <- unique(blas[!is.na(blas) & nzchar(blas)])
    if (!length(blas)) {
        blas <- if (.Platform$OS.type == "windows") {
            paste("Windows R BLAS:", file.path(R.home("bin"), "Rblas.dll"),
                  "(sessionInfo/extSoftVersion do not report its vendor)")
        } else "BLAS vendor/path not reported by sessionInfo or extSoftVersion"
    }
    environment <- Sys.getenv(c("OMP_THREAD_LIMIT", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"))
    writeLines(c("RcppAMGCL versus Matrix::Cholesky", paste("Case:", config$id),
                 paste("UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
                 paste("Repository commit:", paste(source_commit, collapse = " ")),
                 "", "Protocol:",
                 "Fresh process per matrix/size/method/thread setting; timeout 300 seconds.",
                 "Warmup: one fresh setup and first RHS solve, including validation.",
                 "Every measured setup starts with a fresh dsCMatrix and empty factors slot.",
                 "Matrix::Cholesky uses perm=TRUE, LDL=FALSE, super=NA; BLAS env requests one thread.",
                 "AMGCL uses CG, tol=1e-8, explicit requested threads; internal residual is timed.",
                 "RHS generation, matrix construction, and external validation are outside timers.",
                 "Total = setup + first solve; batch total = setup + sum of all RHS solves.",
                 "External validation: relative residual <=1e-8; relative solution error <=1e-5.",
                 "peak_rss_mib is each configuration's final process high-water RSS, repeated on its rows.",
                 "RSS includes R, inputs, warmup, all repetitions, workspace, and external validation.",
                 "R's configured compiler is recorded below; it need not identify a prebuilt package's compiler.",
                 "", capture.output(session), "", capture.output(blas),
                 "", capture.output(Sys.info()[c("sysname", "release", "version", "machine")]),
                 "", capture.output(info), "", capture.output(environment),
                 "", compiler_info(), "", "Artifact MD5 hashes:", capture.output(hashes)), path)
}

run_child <- function(config_file) {
    config <- readRDS(config_file)
    .libPaths(config$libraries)
    rows <- do.call(rbind, lapply(seq_len(config$repetitions), function(i) new_row(config, i)))
    flush <- function() write.csv(rows, config$result_file, row.names = FALSE)
    flush()
    peak_rss <- function() NA_real_
    initialization_error <- tryCatch({
        if (!requireNamespace("RcppAMGCL", quietly = TRUE) ||
            !requireNamespace("Matrix", quietly = TRUE)) stop("RcppAMGCL and Matrix must be installed")
        dll <- dyn.load(config$memory_dll)
        peak_rss <- function() {
            value <- .Call("RcppAMGCL_benchmark_peak_rss", PACKAGE = dll[["name"]])
            if (length(value) != 1L || !is.finite(value) || value <= 0) {
                stop("The memory helper did not return a positive finite process peak RSS")
            }
            as.double(value) / 1024^2
        }
        now <- function() .Call("RcppAMGCL_benchmark_time", PACKAGE = dll[["name"]])
        stopifnot(is.finite(now()))
        invisible(peak_rss())
        info <- RcppAMGCL::amgcl_info()
        write_metadata(config$session_file, config, info)
        A <- poisson(config$minimum_n, config$dimension)
        stopifnot(inherits(A, "dsCMatrix"), length(A@factors) == 0L)
        indexes <- seq_len(nrow(A))
        references <- lapply(seq_len(config$rhs_count), function(j) {
            sin(indexes / (17 + j)) + cos(indexes / (113 + j))
        })
        rhs <- lapply(references, function(x) as.double(A %*% x))
        rows$nnz <- 2 * length(A@x) - nrow(A)
        rows$openmp <- info$openmp
        rows$max_threads <- info$max_threads
        rows$thread_limit <- info$thread_limit

        setup <- function(matrix) {
            if (config$method == "AMGCL_CG") {
                RcppAMGCL::amgcl_setup(matrix, solver = "cg", tol = 1e-8, threads = config$threads)
            } else {
                Matrix::Cholesky(matrix, perm = TRUE, LDL = FALSE, super = NA)
            }
        }
        solve_rhs <- function(factor, b) {
            if (config$method == "AMGCL_CG") {
                result <- RcppAMGCL::amgcl_solve(factor, b, threads = config$threads)
                list(x = as.double(result$x), converged = result$converged,
                     residual = result$residual, iterations = result$iterations,
                     threads = result$threads)
            } else {
                list(x = as.double(Matrix::solve(factor, b, system = "A")), converged = TRUE,
                     residual = NA_real_, iterations = NA_integer_, threads = 1L)
            }
        }

        # Keep warmup references out of subsequent measured repetitions.
        local({
            warm <- setup(fresh_matrix(A))
            result <- solve_rhs(warm, rhs[[1L]])
            metrics <- validate(A, rhs[[1L]], result$x, references[[1L]])
            if (!result$converged || !check_accuracy(metrics)) {
                stop(sprintf("Warmup validation failed: residual=%.17g; solution_error=%.17g; converged=%s",
                             metrics[["residual"]], metrics[["solution_error"]], result$converged))
            }
        })

        run_repetition <- function(repetition) {
            row <- rows[repetition, , drop = FALSE]
            solve_times <- rep(NA_real_, config$rhs_count)
            residuals <- errors <- reported <- rep(NA_real_, config$rhs_count)
            iterations <- rep(NA_integer_, config$rhs_count)
            row$status <- "ok"
            tryCatch({
                matrix <- fresh_matrix(A)
                started <- now()
                factor <- setup(matrix)
                row$setup_seconds <- now() - started
                row$factor_class <- paste(class(factor), collapse = "/")
                for (j in seq_len(config$rhs_count)) {
                    started <- now()
                    result <- solve_rhs(factor, rhs[[j]])
                    solve_times[j] <- now() - started
                    row$completed_rhs <- j
                    row$effective_threads <- result$threads
                    # This multiplication and accuracy check are outside the solve timer.
                    metrics <- validate(A, rhs[[j]], result$x, references[[j]])
                    residuals[j] <- metrics[["residual"]]
                    errors[j] <- metrics[["solution_error"]]
                    reported[j] <- result$residual
                    iterations[j] <- result$iterations
                    if (!result$converged || !check_accuracy(metrics)) {
                        row$status <- "validation_failed"
                        row$message <- append_message(row$message, sprintf(
                            "RHS %d: residual=%.17g; solution_error=%.17g; converged=%s",
                            j, residuals[j], errors[j], result$converged))
                    }
                }
            }, error = function(error) {
                row$status <<- "error"
                row$message <<- append_message(row$message, conditionMessage(error))
            })
            row$first_solve_seconds <- solve_times[1L]
            row$total_seconds <- row$setup_seconds + row$first_solve_seconds
            if (all(is.finite(solve_times))) row$solve_all_seconds <- sum(solve_times)
            row$batch_total_seconds <- row$setup_seconds + row$solve_all_seconds
            if (any(!is.na(residuals))) row$max_residual <- max(residuals, na.rm = TRUE)
            if (any(!is.na(errors))) row$max_solution_error <- max(errors, na.rm = TRUE)
            if (any(!is.na(reported))) row$max_reported_residual <- max(reported, na.rm = TRUE)
            if (any(!is.na(iterations))) row$max_iterations <- max(iterations, na.rm = TRUE)
            row$checkpoint_peak_rss_mib <- peak_rss()
            row
        }
        for (repetition in seq_len(config$repetitions)) {
            invisible(gc())
            rows[repetition, ] <- run_repetition(repetition)
            flush()
        }
        NULL
    }, error = identity)
    if (inherits(initialization_error, "error")) {
        pending <- rows$status == "not_run"
        rows$status[pending] <- "error"
        rows$message[pending] <- conditionMessage(initialization_error)
        # Also retain an error after completed samples, for example a failed
        # final memory checkpoint, rather than silently returning success.
        if (!any(pending)) {
            rows$status[nrow(rows)] <- "error"
            rows$message[nrow(rows)] <- append_message(rows$message[nrow(rows)],
                                                       conditionMessage(initialization_error))
        }
    }
    final_peak <- tryCatch(peak_rss(), error = function(error) NA_real_)
    rows$peak_rss_mib <- final_peak
    rows$peak_rss_complete <- is.finite(final_peak)
    rows$child_exit_status <- if (all(rows$status == "ok")) 0L else 1L
    flush()
    print(rows[, c("case_id", "repetition", "total_seconds", "batch_total_seconds",
                   "max_residual", "max_solution_error", "peak_rss_mib", "status", "message")],
          row.names = FALSE)
    if (any(rows$status != "ok")) quit(save = "no", status = 1L)
    invisible(NULL)
}

run_parent <- function(args, script) {
    if (length(args) < 1L || length(args) > 5L) {
        stop("Usage: Rscript tools/compare-solvers.R output.csv [sizes] [repetitions] [rhs_count] [memory_dll]",
             call. = FALSE)
    }
    output <- args[[1L]]
    sizes <- if (length(args) >= 2L) strsplit(args[[2L]], ",", fixed = TRUE)[[1L]] else
        c("10000", "40000", "100000")
    sizes <- unique(vapply(sizes, positive_integer, integer(1), name = "size"))
    repetitions <- if (length(args) >= 3L) positive_integer(args[[3L]], "repetitions") else 3L
    rhs_count <- if (length(args) >= 4L) positive_integer(args[[4L]], "rhs_count") else 10L
    memory_dll <- if (length(args) >= 5L) args[[5L]] else
        Sys.getenv("RCPPAMGCL_BENCHMARK_MEMORY_DLL",
                   unset = file.path(dirname(script), paste0("benchmark-memory", .Platform$dynlib.ext)))
    if (!file.exists(memory_dll)) {
        stop("Compile tools/benchmark-memory.cpp with R CMD SHLIB and provide its DLL path as argument five ",
             "or RCPPAMGCL_BENCHMARK_MEMORY_DLL", call. = FALSE)
    }
    memory_dll <- normalizePath(memory_dll, winslash = "/", mustWork = TRUE)
    inherited_limit <- Sys.getenv("OMP_THREAD_LIMIT", unset = "")
    limit <- if (nzchar(inherited_limit)) min(8L, positive_integer(inherited_limit, "OMP_THREAD_LIMIT")) else 8L
    Sys.setenv(OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1",
               VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1", OMP_THREAD_LIMIT = as.character(limit))
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    output <- normalizePath(output, winslash = "/", mustWork = FALSE)
    cases_dir <- Sys.getenv("RCPPAMGCL_BENCHMARK_CASE_DIR", unset =
        file.path(dirname(script), "..", ".local", "compare-solvers", paste0(basename(output), ".cases")))
    dir.create(cases_dir, recursive = TRUE, showWarnings = FALSE)
    cases_dir <- normalizePath(cases_dir, winslash = "/", mustWork = TRUE)
    executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    configs <- list()
    methods <- data.frame(method = c(rep("AMGCL_CG", 4L), "Matrix_Cholesky"),
                           threads = c(1L, 2L, 4L, 8L, 1L))
    for (dimension in c(2L, 3L)) for (size in sizes) {
        for (method_index in seq_len(nrow(methods))) {
            method <- methods$method[method_index]
            threads <- methods$threads[method_index]
            id <- sprintf("poisson%dd-n%d-%s-t%d", dimension, size, method, threads)
            configs[[length(configs) + 1L]] <- list(id = id, dimension = dimension, minimum_n = size,
                method = method, threads = threads, repetitions = repetitions, rhs_count = rhs_count,
                script = script, libraries = .libPaths(), memory_dll = memory_dll,
                result_file = file.path(cases_dir, paste0(id, ".csv")),
                session_file = file.path(cases_dir, paste0(id, ".session.txt")))
        }
    }
    accumulated <- list()
    for (i in seq_along(configs)) {
        config <- configs[[i]]
        config_file <- file.path(cases_dir, paste0(config$id, ".rds"))
        log_file <- file.path(cases_dir, paste0(config$id, ".log"))
        saveRDS(config, config_file)
        # Avoid mistaking an earlier run's output for a timed-out child result.
        if (file.exists(config$result_file)) unlink(config$result_file)
        message(sprintf("[%d/%d] %s", i, length(configs), config$id))
        child_error <- ""
        status <- tryCatch(suppressWarnings(system2(executable,
            c("--vanilla", shQuote(script), "--child", shQuote(config_file)),
            stdout = log_file, stderr = log_file, timeout = 300)),
            error = function(error) { child_error <<- conditionMessage(error); 1L })
        rows <- tryCatch(read.csv(config$result_file, stringsAsFactors = FALSE,
                                  colClasses = c(message = "character")), error = function(error) NULL)
        if (is.null(rows) || nrow(rows) != repetitions) {
            rows <- do.call(rbind, lapply(seq_len(repetitions), function(j) new_row(config, j)))
        }
        pending <- rows$status == "not_run"
        if (any(pending)) {
            rows$status[pending] <- if (identical(as.integer(status), 124L)) "timeout" else "child_failed"
            rows$message[pending] <- if (nzchar(child_error)) child_error else
                paste("Child exit status", status, "- see", basename(log_file))
        }
        if (status != 0L && all(rows$status == "ok")) {
            rows$status[nrow(rows)] <- "child_failed"
            rows$message[nrow(rows)] <- paste("Child exited with status", status, "after writing samples")
        }
        if (status == 124L) {
            rows$peak_rss_mib <- NA_real_
            rows$peak_rss_complete <- FALSE
        }
        rows$child_exit_status <- as.integer(status)
        accumulated[[i]] <- rows
        write.csv(do.call(rbind, accumulated), output, row.names = FALSE)
    }
    results <- do.call(rbind, accumulated)
    failed <- results$status != "ok" | results$child_exit_status != 0L
    report <- c(paste("Configurations:", length(configs)), paste("Samples:", nrow(results)),
                paste("Failed samples:", sum(failed)), paste("Child timeout: 300 seconds; OMP_THREAD_LIMIT:", limit),
                "Total and batch totals are within-repetition sums, not sums of independent medians.",
                paste("Per-case logs, session/BLAS/compiler metadata, and artifact hashes:", cases_dir),
                "BLAS environment requests one thread; effective vendor thread counts are not independently queried.",
                "Memory is whole-process high-water RSS including warmup and validation, not incremental factor memory.")
    if (any(failed)) report <- c(report, "", capture.output(print(results[failed,
        c("case_id", "repetition", "status", "message")], row.names = FALSE)))
    writeLines(report, paste0(output, ".report.txt"))
    sessions <- vapply(configs, function(config) config$session_file, character(1))
    existing_sessions <- sessions[file.exists(sessions)]
    if (length(existing_sessions)) {
        writeLines(c("Representative fresh child session; per-configuration settings and measurements are in the CSV.",
                     "All child sessions/logs are retained at:", cases_dir, "",
                     readLines(existing_sessions[[1L]], warn = FALSE)), paste0(output, ".session.txt"))
    }
    message("Wrote ", output, "; failed samples: ", sum(failed))
    if (any(failed)) quit(save = "no", status = 1L)
    invisible(results)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L && identical(args[[1L]], "--child")) {
    if (length(args) != 2L) stop("Internal child mode requires one config RDS path")
    run_child(args[[2L]])
} else {
    script_arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(script_arg) != 1L) stop("Run this script with Rscript")
    script <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
    run_parent(args, script)
}
