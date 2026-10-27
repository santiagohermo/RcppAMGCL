library(RcppAMGCL)

expect_error <- function(expr, patterns = character()) {
    result <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(result, "error"))
    for (pattern in patterns) {
        stopifnot(grepl(pattern, conditionMessage(result), fixed = TRUE))
    }
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
    rhs <- as.numeric(rhs)
    independently_computed <- sqrt(sum((as.numeric(matrix %*% result$x) - rhs)^2))
    if (any(rhs != 0)) independently_computed <- independently_computed / sqrt(sum(rhs^2))
    stopifnot(is.list(result), is.numeric(result$x), length(result$x) == nrow(matrix),
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

# Symmetric compressed matrices store only one triangle. Expanding either
# triangle must preserve the matrix, including the natural Matrix-valued RHS.
for (triangle in c("U", "L")) {
    symmetric <- Matrix::forceSymmetric(A, uplo = triangle)
    stopifnot(inherits(symmetric, "dsCMatrix"), symmetric@uplo == triangle)
    symmetric_rhs <- symmetric %*% x_true
    check_result(amgcl_solve(symmetric, symmetric_rhs, tol = 1e-10, threads = 2L),
                 A, b, x_true, 2L)
    symmetric_solver <- amgcl_setup(symmetric, tol = 1e-10)
    check_result(amgcl_solve(symmetric_solver, symmetric_rhs), A, b, x_true, 1L)
}

# crossprod() naturally returns symmetric compressed storage without an
# explicit coercion by the caller.
design_n <- 128L
design <- Matrix::sparseMatrix(i = c(seq_len(design_n), seq_len(design_n - 1L)),
                                j = c(seq_len(design_n), 2:design_n),
                                x = c(rep(2, design_n), rep(-1, design_n - 1L)))
normal <- Matrix::crossprod(design)
normal_reference <- cos(seq_len(design_n) / 7)
normal_rhs <- normal %*% normal_reference
stopifnot(inherits(normal, "dsCMatrix"))
check_result(amgcl_solve(normal, normal_rhs, tol = 1e-10),
             normal, normal_rhs, normal_reference, 1L)
normal_solver <- amgcl_setup(normal, tol = 1e-10, threads = 2L)
check_result(amgcl_solve(normal_solver, normal_rhs, threads = 2L),
             normal, normal_rhs, normal_reference, 2L)

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

# Base and Matrix one-column RHS inputs work for both one-off and reusable
# solves. Integer numeric vectors are accepted without changing their values.
rhs_inputs <- list(matrix(b, ncol = 1L),
                   Matrix::Matrix(matrix(b, ncol = 1L), sparse = FALSE),
                   Matrix::Matrix(matrix(b, ncol = 1L), sparse = TRUE),
                   A %*% x_true)
for (rhs in rhs_inputs) {
    check_result(amgcl_solve(A, rhs, tol = 1e-10), A, b, x_true, 1L)
    check_result(amgcl_solve(owned, rhs, threads = 2L), A, b, x_true, 2L)
}
integer_reference <- as.integer(ii %% 7L)
integer_rhs <- as.integer(A %*% integer_reference)
check_result(amgcl_solve(A, integer_rhs, tol = 1e-10),
             A, integer_rhs, integer_reference, 1L)
check_result(amgcl_solve(owned, integer_rhs), A, integer_rhs, integer_reference, 1L)
integer_rhs_matrix <- matrix(integer_rhs, ncol = 1L)
check_result(amgcl_solve(A, integer_rhs_matrix, tol = 1e-10),
             A, integer_rhs, integer_reference, 1L)
check_result(amgcl_solve(owned, integer_rhs_matrix),
             A, integer_rhs, integer_reference, 1L)

printed <- capture.output(print_result <- withVisible(print(owned)))
stopifnot(length(printed) > 0L, any(nzchar(printed)),
          any(grepl(as.character(n), printed, fixed = TRUE)),
          any(grepl(owned$solver, printed, fixed = TRUE)),
          !print_result$visible, identical(print_result$value, owned))

exhausted <- amgcl_solve(A, b, tol = 1e-14, maxiter = 1L)
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
for (rhs in list(rep(TRUE, n), as.complex(b), matrix(TRUE, n, 1L),
                 matrix(as.complex(b), ncol = 1L), matrix(b, ncol = 2L),
                 matrix(b, nrow = 1L), array(b, dim = c(n, 1L, 1L)),
                 Matrix::Matrix(matrix(TRUE, n, 1L), sparse = FALSE),
                 Matrix::Matrix(matrix(TRUE, n, 1L), sparse = TRUE),
                 Matrix::Matrix(matrix(b, ncol = 2L), sparse = TRUE),
                 Matrix::Matrix(matrix(b, ncol = 2L), sparse = FALSE),
                 data.frame(b = b))) {
    expect_error(amgcl_solve(A, rhs))
    expect_error(amgcl_solve(owned, rhs))
}
for (sparse in c(FALSE, TRUE)) {
    invalid_rhs <- Matrix::Matrix(matrix(b, ncol = 1L), sparse = sparse)
    invalid_rhs@x[1L] <- Inf
    expect_error(amgcl_solve(owned, invalid_rhs))
}
expect_error(amgcl_solve(owned, b, solver = "cg"))
expect_error(amgcl_solve(owned, b, tol = 1e-10))
expect_error(amgcl_solve(owned, b, maxiter = 1000L))
invalid <- A
invalid@x[1L] <- NaN
expect_error(amgcl_setup(invalid))
invalid <- A
invalid@i[1L] <- n
expect_error(amgcl_setup(invalid))
invalid <- A
invalid@p[length(invalid@p)] <- 0L
expect_error(amgcl_setup(invalid))

# Validate symmetric storage before expansion; invalid row/triangle data must
# not be silently reinterpreted or repaired by a coercion.
for (triangle in c("U", "L")) {
    symmetric <- Matrix::forceSymmetric(A, uplo = triangle)
    invalid <- symmetric
    invalid@i[1L] <- n
    expect_error(amgcl_setup(invalid))
    invalid <- symmetric
    invalid@p[length(invalid@p)] <- 0L
    expect_error(amgcl_setup(invalid))
    invalid <- symmetric
    invalid@uplo <- "X"
    expect_error(amgcl_setup(invalid))
    invalid <- symmetric
    invalid@x[1L] <- NaN
    expect_error(amgcl_setup(invalid))
}
invalid <- Matrix::forceSymmetric(A, uplo = "U")
invalid@i[1L] <- 1L
expect_error(amgcl_setup(invalid))

missing_diagonal <- Matrix::sparseMatrix(i = c(1L, 2L, 2L), j = c(2L, 1L, 2L),
                                        x = c(1, 1, 2), dims = c(2L, 2L))
zero_diagonal <- Matrix::sparseMatrix(i = c(1L, 2L, 1L, 2L), j = c(1L, 1L, 2L, 2L),
                                     x = c(0, 1, 1, 2), dims = c(2L, 2L))
stopifnot(any(zero_diagonal@x == 0))
expect_error(amgcl_setup(missing_diagonal), "nonzero diagonal")
expect_error(amgcl_setup(zero_diagonal), "nonzero diagonal")

# This singular matrix has nonzero diagonals, so it reaches the coarse
# factorization and exercises native exception context and thread restoration.
singular <- Matrix::sparseMatrix(i = c(1L, 2L, 1L, 2L), j = c(1L, 1L, 2L, 2L),
                                x = rep(1, 4L), dims = c(2L, 2L))
expect_error(amgcl_setup(singular, threads = 2L),
             c("AMGCL setup failed", "Zero sum in skyline_lu factorization"))
stopifnot(identical(amgcl_info()$current_threads, info$current_threads))

null_solver <- owned
null_solver$pointer <- new("externalptr")
expect_error(amgcl_solve(null_solver, b))
foreign_solver <- owned
foreign_solver$pointer <- getDLLRegisteredRoutines("RcppAMGCL")[[".Call"]][[1L]]$address
expect_error(amgcl_solve(foreign_solver, b))
stopifnot(identical(amgcl_info()$current_threads, info$current_threads))

# OpenMP reads its environment at process startup. A fresh child tests a hard
# one-thread limit and preservation of the calling task's observed setting.
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
                 "if (before$openmp) stopifnot(before$thread_limit == 1L, before$current_threads >= 1L, before$current_threads <= 2L)"),
               child)
    output <- system2(file.path(R.home("bin"), "Rscript"),
                      c("--vanilla", shQuote(child)), stdout = TRUE, stderr = TRUE)
    status <- attr(output, "status")
    if (!is.null(status) && status != 0L) stop(paste(output, collapse = "\n"))
})
