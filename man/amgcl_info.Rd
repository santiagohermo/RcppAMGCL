\name{amgcl_info}
\alias{amgcl_info}
\title{Report AMGCL and OpenMP Capabilities}
\usage{amgcl_info()}
\description{
Report the bundled AMGCL version and OpenMP capabilities of the installed
package. Querying this information does not change OpenMP settings.
}
\value{
A list with components:
\describe{
  \item{openmp}{Logical; whether this package was compiled with OpenMP.}
  \item{max_threads}{Available thread capacity, capped by processor count
    and the OpenMP thread limit, and at least one.}
  \item{thread_limit}{OpenMP thread limit, reflecting \env{OMP_THREAD_LIMIT}.
    One for a serial build.}
  \item{current_threads}{Current calling-task OpenMP thread-count setting.
    One for a serial build.}
  \item{amgcl_version}{Version of the bundled AMGCL library.}
}
}
\details{
Parallelism is opt-in: \code{amgcl_setup()} and \code{amgcl_solve()} default
to one thread even if OpenMP is enabled. Set \env{OMP_THREAD_LIMIT} before
starting R to impose a process-wide OpenMP thread limit. Runtime resource
constraints can result in fewer worker threads than requested.
}
\seealso{\code{\link{amgcl_setup}}, \code{\link{amgcl_solve}}}
\examples{amgcl_info()}
