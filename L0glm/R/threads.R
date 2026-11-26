#' Query or set BLAS and OpenMP thread limits
#'
#' This helper is for code that coordinates several parallel operations in the
#' same R process. For an individual fit, use the fitting function's
#' `ncores` and `nthreads` arguments.
#'
#' @param blas,openmp A positive integer to set the corresponding limit.
#'   `NULL` or `NA_integer_` leaves that limit as it is.
#' @return A named list with the previous limits `blas` and `openmp`, returned
#'   invisibly when setting a limit. `NA_integer_` means the corresponding
#'   runtime control is unavailable. Pass the list back with
#'   `do.call(l0glmThreads, previous)` to restore the limits.
#' @details Runtime BLAS control uses the thread API of the BLAS loaded by R.
#'   Supported APIs are OpenBLAS, oneMKL, FlexiBLAS and BLIS configured with
#'   automatic threading. Other BLAS builds use their existing configuration.
#'   OpenMP control applies to the runtime used by L0glm when compiled with
#'   OpenMP. Call this helper from the main thread of the R process.
#' @md
#' @export
l0glmThreads <- function(blas = NULL, openmp = NULL) {
  blas <- L0glm_thread_argument(blas, "blas")
  openmp <- L0glm_thread_argument(openmp, "openmp")
  previous <- l0glm_threads_cpp(blas, openmp)
  if (blas == 0L && openmp == 0L) previous else invisible(previous)
}

L0glm_thread_argument <- function(value, name) {
  if (is.null(value) || identical(value, NA_integer_)) return(0L)
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 1 || value > .Machine$integer.max || value != floor(value)) {
    stop(name, " must be a positive integer, NULL or NA_integer_", call. = FALSE)
  }
  as.integer(value)
}
