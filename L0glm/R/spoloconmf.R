#' Sparse Poisson log-concave nonnegative matrix factorization
#'
#' Recover compounds with nonnegative spectra and log-concave elution profiles
#' from a matrix of Poisson observations. A sparse deconvolution fit supplies
#' candidate profiles; spectral clustering, local fusion and Poisson updates
#' recover broader or asymmetric compounds while retaining every outcome channel.
#'
#' @param X Nonnegative design matrix whose columns are contiguous, log-concave
#'   candidate profiles, for example shifted copies of a point-spread function
#'   from \code{\link{convmatrix}}. Dense and sparse matrices are supported.
#'   Columns are normalized to unit sum internally.
#' @param y Nonnegative count matrix, with scans in rows and outcome channels
#'   in columns, or a count vector for a single channel.
#' @param background Optional nonnegative background basis with \code{nrow(y)}
#'   rows, for example a low-df B-spline basis. Its coefficients are nonnegative
#'   and unpenalized. \code{NULL} fits without a background.
#' @param lambda Information-criterion preset (\code{"gic"}, \code{"bic"},
#'   \code{"aic"}, \code{"ebic"} or \code{"mbic"}), or a finite nonnegative
#'   numeric compound penalty. See Details for its scale.
#' @param ncores Number of OpenMP threads. \code{NULL} uses the physical cores
#'   available on the machine.
#' @param control Named list of optional controls: \code{broad_repeats = 2L},
#'   \code{full_sweeps = 1L}, \code{final_polishes = 2L},
#'   \code{native_sweeps = 8L}, \code{block_sweeps = 12L}, and
#'   \code{psf_width = NULL}. The last is the typical PSF full width at half
#'   maximum in scans; by default it is inferred from the middle design column.
#' @param verbose Print the time, objective change and compound count at each
#'   stage.
#'
#' @details
#' The mean is \eqn{\mu = HS + B A}, where each column of \eqn{H} is a
#' unit-sum elution profile, rows of \eqn{S} are spectra and \eqn{BA} is the
#' smooth background. Profiles have contiguous support and satisfy
#' \eqn{\Delta^2 \log H_{\cdot k} \leq 0} on their positive support.
#'
#' The selection objective is \eqn{-2\ell(y;\mu) + \lambda_R K}, with
#' \eqn{K} the number of compounds. An IC preset sets
#' \eqn{\lambda_R = m\lambda_{IC}}, where \eqn{m} is the number of outcome
#' channels. A numeric \code{lambda} specifies \eqn{\lambda_R} directly;
#' the initial grouped deconvolution then uses \eqn{\lambda_R/m}.
#' This compound-count criterion assigns one count per compound. Spectral
#' entries and profile complexity have no additional penalty.
#'
#' Clustering proposes fusions from normalized spectra and overlapping
#' profiles without fitting a model. Proposals are subsequently checked with
#' local Poisson updates, holding the rest of the fitted signal in an offset.
#' Every completed stage is accepted only when the objective on all original
#' observations improves. The stage ledger records rejected proposals too.
#' Sparse profiles remain sparse in the native spectral updates, which run
#' across channels using OpenMP.
#'
#' @return An object of class \code{"spoloconmf"} containing sparse profiles
#'   \code{H}, spectra \code{S}, background coefficients \code{alpha}, fitted
#'   means \code{mu}, the background basis, the objective and stage ledger
#'   \code{history}, timing information \code{seconds}, and resolved controls.
#'   \code{fitted()} returns the fitted mean matrix.
#' @seealso \code{\link{L0glm_sparse_deco}}, \code{\link{L0glm_blind_deco}},
#'   \code{\link{L0glm.fit}}
#' @example examples/spoloconmf.R
#' @export
spoloconmf <- function(X, y, background = NULL, lambda = "gic", ncores = NULL,
                      control = list(), verbose = FALSE) {
  started <- proc.time()[["elapsed"]]
  y <- as.matrix(y)
  if (!is.numeric(y) || !length(y) || any(!is.finite(y)) || any(y < 0) ||
      any(abs(y - round(y)) > 1e-7))
    stop("y must contain finite nonnegative Poisson counts.", call. = FALSE)
  storage.mode(y) <- "double"
  if (!is.matrix(X) && !inherits(X, "Matrix")) X <- as.matrix(X)
  values <- if (inherits(X, "sparseMatrix")) X@x else X
  if (!is.numeric(values) || nrow(X) != nrow(y) || !ncol(X) ||
      any(!is.finite(values)) || any(values < 0))
    stop("X must be a finite nonnegative design with nrow(y) rows.", call. = FALSE)
  if (inherits(X, "sparseMatrix")) X <- methods::as(X, "dgCMatrix")
  mass <- as.numeric(Matrix::colSums(X))
  if (any(mass <= 0)) stop("Each column of X must have positive mass.", call. = FALSE)
  for (j in seq_len(ncol(X))) {
    if (inherits(X, "dgCMatrix")) {
      at <- seq.int(X@p[j] + 1L, X@p[j + 1L]); positive <- X@x[at] > 0
      rows <- X@i[at][positive]; h <- X@x[at][positive]
    } else {
      rows <- which(X[, j] > 0); h <- as.numeric(X[rows, j])
    }
    if (any(diff(rows) != 1L) || any(diff(log(h), differences = 2L) > 1e-6))
      stop("Columns of X must have contiguous log-concave positive support.", call. = FALSE)
  }
  if (max(abs(mass - 1)) > 1e-12) {
    X <- if (inherits(X, "sparseMatrix")) X %*% Matrix::Diagonal(x = 1 / mass) else
      sweep(as.matrix(X), 2L, mass, "/")
  }
  B <- if (is.null(background)) matrix(0, nrow(y), 0L) else as.matrix(background)
  if (!is.numeric(B) || nrow(B) != nrow(y) || any(!is.finite(B)) || any(B < 0))
    stop("background must be a finite nonnegative basis with nrow(y) rows.", call. = FALSE)
  storage.mode(B) <- "double"
  criterion <- NULL
  if (is.character(lambda)) {
    criterion <- L0glm_ic_name(lambda)
    scalar <- preset.lambda(criterion, p = ncol(X), n = nrow(y))
    lambda <- ncol(y) * scalar
  } else {
    if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0)
      stop("lambda must be an IC preset or a finite nonnegative compound penalty.", call. = FALSE)
    scalar <- lambda / ncol(y)
  }
  defaults <- list(broad_repeats = 2L, full_sweeps = 1L, final_polishes = 2L,
                   native_sweeps = 8L, block_sweeps = 12L, psf_width = NULL)
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      any(!nzchar(names(control))) || anyDuplicated(names(control)))))
    stop("control must be a named list.", call. = FALSE)
  if (length(setdiff(names(control), names(defaults))))
    stop("Unknown spoloconmf control.", call. = FALSE)
  control <- utils::modifyList(defaults, control, keep.null = TRUE)
  for (key in setdiff(names(defaults), "psf_width")) {
    value <- control[[key]]
    minimum <- if (key %in% c("native_sweeps", "block_sweeps")) 1L else 0L
    L0glm_deco_length(value, key, minimum)
  }
  if (is.null(control$psf_width)) {
    h <- as.numeric(X[, max(1L, floor(ncol(X) / 2L))])
    control$psf_width <- sum(h >= max(h) / 2)
  }
  L0glm_deco_length(control$psf_width, "psf_width", 1L)
  if (is.null(ncores)) {
    ncores <- L0glm_detect_physical_cores()
    if (is.na(ncores)) ncores <- 1L
  }
  L0glm_deco_length(ncores, "ncores", 1L)
  ncores <- as.integer(ncores)
  d <- list(X = X, Y = y, B = B, lambda = scalar, lambda_rank = lambda,
            psf_width = control$psf_width)
  current <- spolo_initial(d, ncores)
  initial <- current
  initial_seconds <- current$seconds
  history <- list(data.frame(stage = "initial", seconds = initial_seconds, delta = 0,
                 accepted = TRUE, components = ncol(current$H), clustering_seconds = 0))
  clustering <- 0
  step <- function(label, expression) {
    elapsed <- system.time(candidate <- force(expression))[["elapsed"]]
    delta <- spolo_delta_cpp(y, current$mu, candidate$mu) +
      lambda * (ncol(candidate$H) - ncol(current$H))
    accepted <- is.finite(delta) && delta <= 1e-6
    cluster <- if ("clustering" %in% names(candidate$seconds)) candidate$seconds[["clustering"]] else 0
    clustering <<- clustering + cluster
    history[[length(history) + 1L]] <<- data.frame(stage = label, seconds = elapsed,
      delta = delta, accepted = accepted, components = ncol(candidate$H), clustering_seconds = cluster)
    if (accepted) current <<- candidate
    if (isTRUE(verbose)) cat(label, ": ", round(elapsed, 2), " s; objective change ",
      signif(delta, 5), "; ", if (accepted) "accepted" else "rejected", "; ",
      ncol(current$H), " compounds\n", sep = "")
    invisible(NULL)
  }
  if (ncol(current$H)) {
    step("fusion_1", spolo_fast(current, d, ncores, control$native_sweeps, control$block_sweeps))
    for (pass in seq_len(control$broad_repeats)) {
      step(paste0("broad_", pass), spolo_broad_polish(current, d, ncores))
      step(paste0("fusion_", pass + 1L), spolo_fast(current, d, ncores,
        control$native_sweeps, control$block_sweeps, background = pass == 1L))
    }
    if (control$full_sweeps) step("full_sweep", spolo_sweep(current, d, control$full_sweeps, ncores))
    step("completion", spolo_complete(current, d, ncores, control$block_sweeps))
    for (pass in seq_len(control$final_polishes)) {
      step(paste0("polish_broad_", pass), spolo_broad_polish(current, d, ncores))
      step(paste0("polish_spectra_", pass), spolo_spectral_polish(current, d, ncores, control$native_sweeps))
    }
  }
  labels <- if (nrow(current$S)) paste0("compound", seq_len(nrow(current$S))) else character()
  dimnames(current$S) <- list(labels, colnames(y))
  colnames(current$H) <- rownames(current$S)
  dimnames(current$mu) <- dimnames(y)
  current$history <- do.call(rbind, history)
  current$background <- B
  current$seconds <- c(initial = initial_seconds, clustering = clustering,
                       total = proc.time()[["elapsed"]] - started)
  current$objective <- list(criterion = criterion, lambda_rank = lambda, compounds = ncol(current$H),
    deviance = sum(poisson()$dev.resids(as.numeric(y), as.numeric(current$mu), 1)),
    loglik = sum(stats::dpois(y, current$mu, log = TRUE)),
    spectral_entry_penalty = 0, profile_complexity_penalty = 0)
  current$objective$IC_compounds <- -2 * current$objective$loglik + lambda * ncol(current$H)
  current$objective$rank_score <- current$objective$deviance + lambda * ncol(current$H)
  current$diagnostics <- list(fitting_calls_in_clustering = 0L,
    complete_delta = spolo_delta_cpp(y, initial$mu, current$mu) +
      lambda * (ncol(current$H) - ncol(initial$H)))
  current$control <- c(list(ncores = ncores), control)
  class(current) <- "spoloconmf"
  current
}

#' @rdname spoloconmf
#' @param x,object A fitted \code{spoloconmf} object.
#' @param ... Additional arguments to S3 methods.
#' @method print spoloconmf
#' @export
print.spoloconmf <- function(x, ...) {
  cat("Sparse Poisson log-concave NMF\n", ncol(x$H), " compounds, ", ncol(x$S),
      " channels; ", format(round(x$seconds[["total"]], 2)), " s\n",
      "Compound-count IC: ", format(x$objective$IC_compounds, digits = 8), "\n", sep = "")
  invisible(x)
}

#' @rdname spoloconmf
#' @export
fitted.spoloconmf <- function(object, ...) object$mu
