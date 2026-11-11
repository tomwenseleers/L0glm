# Canonical observation matrices for row-wise cross-validation.
L0glm_tune_matrix <- function(value, n, m, default, name) {
  if (is.null(value)) return(matrix(default, n, m))
  if (length(value) == 1L) return(matrix(as.numeric(value), n, m))
  if (is.null(dim(value))) value <- matrix(value, ncol = 1L)
  if (nrow(value) != n || !ncol(value) %in% c(1L, m))
    stop(name, " must have nrow(X) rows and one or ncol(y) columns.", call. = FALSE)
  matrix(as.numeric(value), n, m)
}

L0glm_cv_folds <- function(weights, nfolds, foldid) {
  n <- nrow(weights)
  active <- rowSums(weights) > 0
  if (!is.null(nfolds) && (!is.numeric(nfolds) || length(nfolds) != 1L ||
      !is.finite(nfolds) || nfolds < 2 || nfolds != floor(nfolds)))
    stop("nfolds must be an integer between 2 and the number of positive-weight rows.",
         call. = FALSE)
  if (!is.null(foldid)) {
    if (!is.numeric(foldid) || length(foldid) != n || anyNA(foldid) ||
        any(!is.finite(foldid) | foldid < 1 | foldid != floor(foldid)))
      stop("foldid must contain one positive integer per observation row.", call. = FALSE)
    labels <- sort(unique(foldid))
    if (!identical(as.integer(labels), seq_along(labels)))
      stop("foldid labels must be consecutive integers starting at 1.", call. = FALSE)
    if (is.null(nfolds)) nfolds <- length(labels)
    if (nfolds != length(labels))
      stop("nfolds must match the number of folds in foldid.", call. = FALSE)
  }
  if (!is.numeric(nfolds) || length(nfolds) != 1L || !is.finite(nfolds) ||
      nfolds < 2 || nfolds != floor(nfolds) || nfolds > sum(active))
    stop("nfolds must be an integer between 2 and the number of positive-weight rows.",
         call. = FALSE)
  nfolds <- as.integer(nfolds)
  if (is.null(foldid)) {
    foldid <- integer(n)
    foldid[active] <- sample(rep(seq_len(nfolds), length.out = sum(active)))
    foldid[!active] <- rep(seq_len(nfolds), length.out = sum(!active))
  }
  if (any(vapply(seq_len(nfolds), function(k) !any(active & foldid == k), logical(1))))
    stop("Every fold must contain at least one positive-weight observation row.", call. = FALSE)
  as.integer(foldid)
}

# Use one common full-data lambda grid and train each model independently.
# Dividing the total likelihood by total prior weight shows that a training fit
# needs lambda * (training weight / full weight) on the fitter's sum-loss scale.
L0glm_cross_validate <- function(args, design, lambdas, foldid, family, mu_floor,
                                pilot_bounds = NULL, screen.path = "independent",
                                reference = NULL, maxit = NULL) {
  nfolds <- max(foldid)
  fold_weights <- vapply(seq_len(nfolds), function(k)
    sum(args$weights[foldid == k, , drop = FALSE]), numeric(1))
  if (ncol(args$weights)==1L) fold_weights <- fold_weights * ncol(args$y)
  fractions <- 1 - fold_weights / sum(fold_weights)
  fold_lambdas <- outer(lambdas, fractions)
  deviance <- matrix(NA_real_, length(lambdas), nfolds,
                     dimnames = list(NULL, paste0("Fold", seq_len(nfolds))))
  starts <- pilots <- preprocessing <- vector("list", nfolds)
  path_order <- L0glm_cv_path_order(lambdas, screen.path, reference, maxit)
  candidates <- vector("list", nfolds)
  interaction_terms <- if(!is.null(args$interactions) && args$interactions>0L)vector('list',nfolds)else NULL
  preclustering <- if (isTRUE(args$precluster)) vector("list", nfolds) else NULL
  large <- as.double(nrow(args$X))*ncol(args$X) > 1e8
  if (large) invisible(gc())
  for (k in seq_len(nfolds)) {
    validation <- foldid == k
    training <- !validation
    train_args <- args
    train_args$X <- args$X[training, , drop = FALSE]
    train_args$y <- args$y[training, , drop = FALSE]
    train_args$weights <- args$weights[training, , drop = FALSE]
    train_args$offset <- args$offset[training, , drop = FALSE]
    initialization <- if (is.null(pilot_bounds)) list(start = args$start, pilot = NULL) else
      L0glm_pilot_start(train_args, pilot_bounds * fractions[k])
    train_args$start <- initialization$start
    starts[k] <- list(initialization$start)
    pilots[k] <- list(initialization$pilot)
    w <- args$weights[validation, , drop = FALSE]
    yy <- args$y[validation, , drop = FALSE]
    off <- L0glm_observation_broadcast(args$offset[validation, , drop = FALSE])
    for (i in path_order) {
      if (!is.null(maxit)) train_args$maxit <- maxit[i]
      fit <- tryCatch(do.call(l0glmFit, c(train_args, list(lambda = fold_lambdas[i, k]))),
        error = function(e) stop(sprintf("CV fit failed in fold %d at lambda %.4g: %s",
                                         k, lambdas[i], conditionMessage(e)), call. = FALSE))
      if (i == path_order[1L]) {
        preprocessing[[k]] <- fit$preprocessing
        if (!is.null(preclustering)) preclustering[[k]] <- fit$preclustering
        train_args <- L0glm_precluster_reuse(train_args, fit)
        train_args <- L0glm_path_screen_reuse(train_args, fit, screen.path)
        candidates[k] <- list(train_args$control.solver$.screen.candidates)
      }
      active <- which(Matrix::rowSums(abs(fit$coefficients)) != 0)
      if(!is.null(interaction_terms)) {
        if(is.null(interaction_terms[[k]]))interaction_terms[[k]] <- vector('list',length(lambdas))
        interaction_terms[[k]][i] <- list(fit$interactions$terms)
      }
      eta <- if(!is.null(fit$interactions)||inherits(args$X,'L0glm_data'))predict(fit,args$X[validation,,drop=FALSE],
        type='link',offset=off) else as.matrix(design[validation, active, drop = FALSE] %*%
        fit$coefficients[active, , drop = FALSE]) + off
      mu <- family$linkinv(eta)
      dim(mu) <- dim(eta)
      if (family$family %in% c("poisson", "quasipoisson", "Gamma", "inverse.gaussian"))
        mu <- pmax(mu, mu_floor)
      # Family C kernels require equal-length inputs, including shared weights.
      loss <- vapply(seq_len(ncol(yy)),function(j) {
        ww <- L0glm_observation_column(w,j); keep <- ww>0
        sum(family$dev.resids(yy[keep,j],mu[keep,j],ww[keep]))
      },numeric(1))
      deviance[i, k] <- if (anyNA(loss)) Inf else sum(loss) / fold_weights[k]
    }
    rm(train_args, fit, initialization)
    if (large) invisible(gc())
  }
  proportions <- fold_weights / sum(fold_weights)
  mean_deviance <- as.numeric(deviance %*% proportions)
  se <- rep(NA_real_, length(lambdas))
  finite <- is.finite(mean_deviance)
  if (any(finite)) {
    centered <- sweep(deviance[finite, , drop = FALSE], 1L, mean_deviance[finite], "-")
    se[finite] <- sqrt(as.numeric(centered^2 %*% proportions) / (nfolds - 1L))
  }
  result <- list(lambda = lambdas, mean = mean_deviance, se = se, deviance = deviance,
       foldid = foldid, nfolds = nfolds, fold_weights = fold_weights,
       starts = starts, pilots = pilots, preprocessing = preprocessing,
       lambda_train = fold_lambdas)
  if (!is.null(preclustering)) result$preclustering <- preclustering
  if(!is.null(interaction_terms))result$interaction_terms <- interaction_terms
  if (screen.path == "preset") result$screen.candidates <- candidates
  result
}

L0glm_cv_path_order <- function(lambdas, screen.path, reference, maxit = NULL) {
  indices <- seq_along(lambdas)
  if (screen.path != "preset") return(indices)
  anchor <- which.min(abs(log(lambdas) - log(reference)))
  if (!is.null(maxit)) {
    tied <- which(lambdas == lambdas[anchor])
    anchor <- tied[which.max(maxit[tied])]
  }
  c(anchor, indices[-anchor])
}

# A scalar follows the ordinary path; a vector defines a CV candidate grid.
L0glm_tune_iterations <- function(maxit, do_cv) {
  if (!is.numeric(maxit) || !length(maxit) || anyNA(maxit) ||
      any(!is.finite(maxit) | maxit < 1 | maxit > .Machine$integer.max | maxit != floor(maxit)))
    stop("maxit must contain positive finite integers.", call. = FALSE)
  maxit <- sort(unique(as.integer(maxit)))
  if (length(maxit) > 1L && !do_cv)
    stop('A vector of maxit values requires tune = "cv".', call. = FALSE)
  maxit
}

# The one-SE rule keeps its largest-lambda meaning. Resolve ties at that
# lambda by the largest iteration cap; it is a deterministic convention.
L0glm_cv_choice <- function(cv, lambdas, maxit) {
  best <- which.min(cv$mean)
  eligible <- which(is.finite(cv$mean) & cv$mean <= cv$mean[best] + cv$se[best])
  list(best = best,
       one_se = eligible[order(-lambdas[eligible], -maxit[eligible])][1L])
}
