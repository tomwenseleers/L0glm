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
                                pilot_bounds = NULL) {
  nfolds <- max(foldid)
  fold_weights <- vapply(seq_len(nfolds), function(k)
    sum(args$weights[foldid == k, , drop = FALSE]), numeric(1))
  fractions <- 1 - fold_weights / sum(fold_weights)
  fold_lambdas <- outer(lambdas, fractions)
  deviance <- matrix(NA_real_, length(lambdas), nfolds,
                     dimnames = list(NULL, paste0("Fold", seq_len(nfolds))))
  starts <- pilots <- preprocessing <- vector("list", nfolds)
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
    keep <- w > 0
    yy <- args$y[validation, , drop = FALSE][keep]
    for (i in seq_along(lambdas)) {
      fit <- tryCatch(do.call(L0glm.fit, c(train_args, list(lambda = fold_lambdas[i, k]))),
        error = function(e) stop(sprintf("CV fit failed in fold %d at lambda %.4g: %s",
                                         k, lambdas[i], conditionMessage(e)), call. = FALSE))
      if (i == 1L) preprocessing[[k]] <- fit$preprocessing
      active <- which(Matrix::rowSums(abs(fit$coefficients)) != 0)
      eta <- as.matrix(design[validation, active, drop = FALSE] %*%
        fit$coefficients[active, , drop = FALSE]) +
        args$offset[validation, , drop = FALSE]
      mu <- family$linkinv(eta)
      dim(mu) <- dim(eta)
      if (family$family %in% c("poisson", "quasipoisson", "Gamma", "inverse.gaussian"))
        mu <- pmax(mu, mu_floor)
      loss <- family$dev.resids(yy, mu[keep], w[keep])
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
  list(lambda = lambdas, mean = mean_deviance, se = se, deviance = deviance,
       foldid = foldid, nfolds = nfolds, fold_weights = fold_weights,
       starts = starts, pilots = pilots, preprocessing = preprocessing,
       lambda_train = fold_lambdas)
}
