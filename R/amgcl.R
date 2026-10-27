.amgcl_positive_integer <- function(x, name) {
    if (!is.numeric(x) || is.complex(x) || length(x) != 1L || is.na(x) ||
        !is.finite(x) || x < 1 || x > .Machine$integer.max || x != floor(x)) {
        stop(name, " must be a positive integer scalar.", call. = FALSE)
    }
    as.integer(x)
}

.amgcl_tolerance <- function(tol) {
    if (!is.numeric(tol) || is.complex(tol) || length(tol) != 1L || is.na(tol) ||
        !is.finite(tol) || tol <= 0) {
        stop("tol must be a finite positive numeric scalar.", call. = FALSE)
    }
    as.double(tol)
}

.amgcl_matrix <- function(A) {
    if (inherits(A, "dsCMatrix")) {
        methods::validObject(A)
        # Expand the stored triangle without constructing a dense matrix.
        A <- methods::as(A, "generalMatrix")
    }
    if (!inherits(A, "dgCMatrix")) {
        stop("A must be a Matrix::dgCMatrix or Matrix::dsCMatrix.", call. = FALSE)
    }
    A
}

.amgcl_rhs <- function(b, n) {
    dimensions <- dim(b)
    if (!is.null(dimensions)) {
        if (length(dimensions) != 2L || dimensions[2L] != 1L ||
            !(is.matrix(b) || inherits(b, "dMatrix"))) {
            stop("b must be a numeric vector or a single-column numeric matrix.",
                 call. = FALSE)
        }
        if (inherits(b, "dMatrix")) {
            methods::validObject(b)
            b <- as.double(b)
        } else {
            # Retain the storage type until the numeric/complex checks below.
            dim(b) <- NULL
        }
    }
    if (!is.numeric(b) || is.complex(b) || length(b) != n ||
        anyNA(b) || any(!is.finite(b))) {
        stop("b must contain finite numeric values with length equal to nrow(A).",
             call. = FALSE)
    }
    as.double(b)
}

.amgcl_validate_solver <- function(object) {
    required <- c("pointer", "n", "solver", "tol", "maxiter")
    valid <- is.list(object) && all(required %in% names(object)) &&
        typeof(object$pointer) == "externalptr" &&
        is.numeric(object$n) && !is.complex(object$n) && length(object$n) == 1L &&
        !is.na(object$n) && is.finite(object$n) && object$n >= 1 &&
        object$n <= .Machine$integer.max && object$n == floor(object$n) &&
        is.character(object$solver) && length(object$solver) == 1L &&
        !is.na(object$solver) && object$solver %in% c("cg", "bicgstab") &&
        is.numeric(object$tol) && !is.complex(object$tol) && length(object$tol) == 1L &&
        !is.na(object$tol) && is.finite(object$tol) && object$tol > 0 &&
        is.numeric(object$maxiter) && !is.complex(object$maxiter) &&
        length(object$maxiter) == 1L &&
        !is.na(object$maxiter) && is.finite(object$maxiter) &&
        object$maxiter >= 1 && object$maxiter <= .Machine$integer.max &&
        object$maxiter == floor(object$maxiter)
    if (!valid) {
        stop("Invalid amgcl_solver object; create one with amgcl_setup().",
             call. = FALSE)
    }
    invisible(object)
}

amgcl_setup <- function(A, solver = c("cg", "bicgstab"), tol = 1e-8,
                        maxiter = 1000L, threads = 1L) {
    A <- .amgcl_matrix(A)
    solver <- match.arg(solver)
    tol <- .amgcl_tolerance(tol)
    maxiter <- .amgcl_positive_integer(maxiter, "maxiter")
    threads <- .amgcl_positive_integer(threads, "threads")
    pointer <- amgcl_setup_cpp(A, solver, tol, maxiter, threads)
    structure(list(pointer = pointer, n = nrow(A), solver = solver,
                   tol = tol, maxiter = maxiter), class = "amgcl_solver")
}

amgcl_solve <- function(A, b, solver = c("cg", "bicgstab"), tol = 1e-8,
                        maxiter = 1000L, threads = 1L) {
    threads <- .amgcl_positive_integer(threads, "threads")
    reusable <- inherits(A, "amgcl_solver")
    if (reusable) {
        .amgcl_validate_solver(A)
        if (!missing(solver) || !missing(tol) || !missing(maxiter)) {
            stop("solver, tol and maxiter are fixed by amgcl_setup() for a reusable solver.",
                 call. = FALSE)
        }
        object <- A
        n <- object$n
    } else {
        A <- .amgcl_matrix(A)
        n <- nrow(A)
        if (length(n) != 1L || is.na(n) || n < 1L) {
            stop("A must be a nonempty square matrix.", call. = FALSE)
        }
    }
    b <- .amgcl_rhs(b, n)
    if (!reusable) {
        object <- amgcl_setup(A, solver = match.arg(solver), tol = tol,
                              maxiter = maxiter, threads = threads)
    }
    amgcl_solve_cpp(object$pointer, b, threads)
}

amgcl_info <- function() {
    amgcl_info_cpp()
}

print.amgcl_solver <- function(x, ...) {
    .amgcl_validate_solver(x)
    cat("AMGCL reusable solver\n",
        "  Dimension: ", x$n, " x ", x$n, "\n",
        "  Method: ", x$solver, "; tolerance: ", format(x$tol),
        "; maximum iterations: ", x$maxiter, "\n", sep = "")
    invisible(x)
}
