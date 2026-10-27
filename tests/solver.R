library(RcppAMGCL)

expect_error <- function(expr) {
    result <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(result, "error"))
    invisible(result)
}

# 1024 unknowns exceed the wrapper's coarse-grid threshold, exercising
# aggregation, interpolation, smoothing, and the coarse solve.
side <- 32L
n <- side * side
ii <- seq_len(n)
horizontal <- ii[ii %% side != 0L]
vertical <- ii[ii <= n - side]
A <- Matrix::sparseMatrix(
    i = c(ii, horizontal, horizontal + 1L, vertical, vertical + side),
    j = c(ii, horizontal + 1L, horizontal, vertical + side, vertical),
    x = c(rep(4, n), rep(-1, 2L * (length(horizontal) + length(vertical)))),
    dims = c(n, n)
)
stopifnot(inherits(A, "dgCMatrix"))
x_true <- sin(ii / 13) + cos(ii / 29)
b <- as.numeric(A %*% x_true)

info <- amgcl_info()
stopifnot(is.logical(info$openmp), length(info$openmp) == 1L,
          info$max_threads >= 1, info$thread_limit >= 1,
          info$current_threads >= 1, nzchar(info$amgcl_version))
if (!info$openmp) stopifnot(info$max_threads == 1L)

check_result <- function(result, matrix, rhs, reference, threads) {
    independently_computed <- sqrt(sum((as.numeric(matrix %*% result$x) - rhs)^2))
    if (any(rhs != 0)) independently_computed <- independently_computed / sqrt(sum(rhs^2))
    stopifnot(is.list(result), is.numeric(result$x), length(result$x) == n,
              result$converged, result$iterations >= 0L,
              is.finite(result$error), result$error >= 0,
              isTRUE(all.equal(result$residual, independently_computed,
                               tolerance = 1e-10)),
              result$residual < 1e-8,
              max(abs(result$x - reference)) < 1e-6,
              result$threads == min(threads, info$max_threads),
              identical(amgcl_info()$current_threads, info$current_threads))
}

for (solver in c("cg", "bicgstab")) {
    one <- amgcl_solve(A, b, solver = solver, tol = 1e-10, threads = 1L)
    check_result(one, A, b, x_true, 1L)
    reusable <- amgcl_setup(A, solver = solver, tol = 1e-10, threads = 2L)
    stopifnot(inherits(reusable, "amgcl_solver"),
              identical(amgcl_info()$current_threads, info$current_threads))
    check_result(amgcl_solve(reusable, b, threads = 2L), A, b, x_true, 2L)
    b2 <- as.numeric(A %*% (2 * x_true))
    check_result(amgcl_solve(reusable, b2), A, b2, 2 * x_true, 1L)
    zero <- amgcl_solve(reusable, numeric(n))
    stopifnot(zero$converged, zero$residual == 0, all(zero$x == 0))
    restored <- unserialize(serialize(reusable, NULL))
    expect_error(amgcl_solve(restored, b))
}

# BiCGStab is also exercised with a strictly diagonally dominant nonsymmetric A.
B <- A + Matrix::sparseMatrix(i = ii, j = ii, x = rep(0.5, n), dims = c(n, n)) +
    Matrix::sparseMatrix(i = horizontal, j = horizontal + 1L,
                         x = rep(0.1, length(horizontal)), dims = c(n, n))
bb <- as.numeric(B %*% x_true)
check_result(amgcl_solve(B, bb, solver = "bicgstab", tol = 1e-10, threads = 2L),
             B, bb, x_true, 2L)

# Setup takes an owning copy; later mutation or collection of the R matrix is safe.
copy <- A
owned <- amgcl_setup(copy, tol = 1e-10)
copy@x[] <- 7
rm(copy)
invisible(gc())
check_result(amgcl_solve(owned, b), A, b, x_true, 1L)

exhausted <- suppressWarnings(amgcl_solve(A, b, tol = 1e-14, maxiter = 1L))
stopifnot(!exhausted$converged, exhausted$iterations <= 1L,
          exhausted$residual > 1e-14)

expect_error(amgcl_setup(as.matrix(A)))
expect_error(amgcl_setup(A[, -1L]))
expect_error(amgcl_setup(Matrix::sparseMatrix(i = integer(), j = integer(),
                                            x = numeric(), dims = c(0L, 0L))))
expect_error(amgcl_setup(A, solver = "unknown"))
for (value in list(0, -1, NA_real_, Inf, c(1e-8, 1e-9))) {
    expect_error(amgcl_setup(A, tol = value))
}
for (value in list(0, -1, 1.5, NA_real_, Inf, c(1L, 2L))) {
    expect_error(amgcl_setup(A, threads = value))
    expect_error(amgcl_setup(A, maxiter = value))
    expect_error(amgcl_solve(owned, b, threads = value))
}
expect_error(amgcl_solve(owned, b[-1L]))
expect_error(amgcl_solve(owned, replace(b, 1L, NA_real_)))
expect_error(amgcl_solve(owned, replace(b, 1L, Inf)))
expect_error(amgcl_solve(owned, as.character(b)))
invalid <- A
invalid@x[1L] <- NaN
expect_error(amgcl_setup(invalid))
invalid <- A
invalid@i[1L] <- n
expect_error(amgcl_setup(invalid))
invalid <- A
invalid@p[length(invalid@p)] <- 0L
expect_error(amgcl_setup(invalid))

null_solver <- owned
null_solver$pointer <- new("externalptr")
expect_error(amgcl_solve(null_solver, b))
foreign_solver <- owned
foreign_solver$pointer <- getDLLRegisteredRoutines("RcppAMGCL")[[".Call"]][[1L]]$address
expect_error(amgcl_solve(foreign_solver, b))
stopifnot(identical(amgcl_info()$current_threads, info$current_threads))

# OpenMP reads its environment at process startup. A fresh child tests a hard
# one-thread limit and preservation of the calling task's two-thread setting.
local({
    old <- Sys.getenv(c("OMP_THREAD_LIMIT", "OMP_NUM_THREADS"), unset = NA_character_)
    on.exit({
        for (name in names(old)) {
            if (is.na(old[[name]])) Sys.unsetenv(name)
            else do.call(Sys.setenv, setNames(list(old[[name]]), name))
        }
    })
    Sys.setenv(OMP_THREAD_LIMIT = "1", OMP_NUM_THREADS = "2")
    child <- tempfile(fileext = ".R")
    on.exit(unlink(child), add = TRUE)
    libraries <- paste(capture.output(dput(.libPaths())), collapse = "\n")
    writeLines(c(paste0(".libPaths(", libraries, ")"),
                 "library(RcppAMGCL)",
                 "before <- amgcl_info()",
                 "A <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = c(2, 3, 4))",
                 "s <- amgcl_setup(A, threads = 2L)",
                 "z <- amgcl_solve(s, c(2, 3, 4), threads = 2L)",
                 "stopifnot(z$converged, z$threads == 1L, before$max_threads == 1L)",
                 "stopifnot(identical(before$current_threads, amgcl_info()$current_threads))",
                 "if (before$openmp) stopifnot(before$thread_limit == 1L, before$current_threads == 2L)"),
               child)
    output <- system2(file.path(R.home("bin"), "Rscript"),
                      c("--vanilla", shQuote(child)), stdout = TRUE, stderr = TRUE)
    status <- attr(output, "status")
    if (!is.null(status) && status != 0L) stop(paste(output, collapse = "\n"))
})
