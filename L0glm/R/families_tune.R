# Extended-family tuning is isolated from the ordinary-GLM path. Each fit starts
# independently, including theta and preprocessing inside every training fold.
L0glm_ext_cv <- function(args, lambdas, foldid, kind, mu_floor, pilot_bounds) {
  response <- args$y
  mass <- switch(kind, nb = rowSums(args$weights),
    multinomial = args$weights[, 1L] * rowSums(response))
  K <- max(foldid)
  fold_weights <- vapply(seq_len(K), function(k) sum(mass[foldid == k]), numeric(1))
  if (any(fold_weights <= 0) || any(fold_weights >= sum(mass)))
    stop("Each CV fold and its complement need positive likelihood mass.", call. = FALSE)
  fractions <- 1 - fold_weights/sum(mass)
  lambda_train <- outer(lambdas, fractions)
  loss <- matrix(NA_real_, length(lambdas), K)
  starts <- pilots <- preprocessing <- theta <- vector("list", K)
  for (k in seq_len(K)) {
    train <- foldid != k
    a <- args
    for (name in c("X", "y", "weights", "offset")) a[[name]] <- args[[name]][train, , drop = FALSE]
    initialization <- if (is.null(pilot_bounds)) list(start = args$start, pilot = NULL) else
      L0glm_pilot_start(a, pilot_bounds * fractions[k])
    a$start <- initialization$start
    starts[k] <- list(a$start); pilots[k] <- list(initialization$pilot)
    theta[[k]] <- vector("list", length(lambdas))
    for (i in seq_along(lambdas)) {
      fit <- tryCatch(do.call(L0glm.fit, c(a, list(lambda = lambda_train[i, k]))),
        error = function(e) stop(sprintf("CV fold %d, lambda %.4g: %s", k, lambdas[i], conditionMessage(e)), call. = FALSE))
      if (i == 1L) preprocessing[[k]] <- fit$preprocessing
      theta[[k]][i] <- list(fit$theta)
      loss[i, k] <- L0glm.metrics(fit, args$X[!train, , drop = FALSE],
        response[!train, , drop = FALSE], args$weights[!train, , drop = FALSE],
        args$offset[!train, , drop = FALSE], mu_floor)$mean_deviance

    }
  }
  proportions <- fold_weights/sum(fold_weights)
  means <- as.numeric(loss %*% proportions)
  centered <- sweep(loss, 1L, means, "-")
  se <- sqrt(as.numeric(centered^2 %*% proportions)/(K - 1L))
  list(lambda = lambdas, mean = means, se = se, deviance = loss, foldid = foldid,
    nfolds = K, fold_weights = fold_weights, lambda_train = lambda_train,
    starts = starts, pilots = pilots, preprocessing = preprocessing, theta = theta,
    scoring = "held-out deviance per weighted response/trial")
}

L0glm_ext_tune <- function(X, y, family, weights, offset, dispersion, tune,
                            nlambda, prop_range, coefs.true, printbest, plot,
                            mu_floor, nfolds, foldid, lambda_min, start, dots) {
  if (length(nlambda) != 1L || !is.finite(nlambda) || nlambda < 2 || nlambda != floor(nlambda) ||
      length(prop_range) != 1L || !is.finite(prop_range) || prop_range <= 0 || prop_range > 1 ||
      length(lambda_min) != 1L || !is.finite(lambda_min) || lambda_min <= 0 ||
      length(mu_floor) != 1L || !is.finite(mu_floor) || mu_floor < 0)
    stop("Invalid nlambda, prop_range, lambda_min or mu_floor.", call. = FALSE)
  if ("IC" %in% names(dots)) stop("IC was renamed to tune.", call. = FALSE)
  if ("lambda" %in% names(dots)) stop("Use lambda_min to control the tuning grid.", call. = FALSE)
  pilot_start <- identical(start, "pilot")
  if (is.character(start) && !pilot_start)
    stop('start must be NULL, numeric, or "pilot".', call. = FALSE)
  kind <- L0glm_ext_kind(family)
  n <- nrow(X)
  response <- switch(kind, nb = L0glm_ext_nb_response(y, family, n),
    multinomial = L0glm_ext_multinomial_response(y, family, n))
  weights <- L0glm_ext_matrix(weights, n, if (kind == "nb") response$m else 1L, 1, "weights")
  offset <- L0glm_ext_matrix(offset, n, response$m, 0, "offset")
  if (any(!is.finite(weights) | weights < 0) || !any(weights > 0))
    stop("Weights must be finite, nonnegative and not all zero.", call. = FALSE)
  intercept <- if (is.null(dots$intercept)) kind == "multinomial" else isTRUE(dots$intercept)
  p <- L0glm_predictor_count(X, intercept, dots$always.include)
  n_ic <- sum(rowSums(weights) > 0)
  if (n_ic < 1) stop("No positive likelihood mass.", call. = FALSE)
  do_cv <- tune == "cv"
  if (do_cv) foldid <- L0glm_cv_folds(weights, nfolds, foldid) else
    if (!is.null(foldid)) stop("foldid is only used for CV.", call. = FALSE)
  lambda_preset <- preset.lambda(if (do_cv) "gic" else tune, max(1, p), max(2, n_ic))
  reference <- max(lambda_min, lambda_preset)
  args <- c(list(X = X, y = response$y, family = family, weights = weights,
    offset = offset, dispersion = dispersion, start = if (pilot_start) NULL else start), dots)
  cv_args <- args
  pilot_bounds <- if (pilot_start) c(lambda_min,
    max(lambda_min, preset.lambda("gic", max(1, p), max(2, n_ic)))) else NULL
  initialization <- if (pilot_start) L0glm_pilot_start(args, pilot_bounds) else list(start = start, pilot = NULL)
  args$start <- initialization$start
  fit_lambda <- function(lambda) do.call(L0glm.fit, c(args, list(lambda = lambda)))
  preset_fit <- fit_lambda(reference)
  is_zero <- function(fit) fit$metrics$selected == 0
  lambda_max <- 2 * reference
  endpoint <- fit_lambda(lambda_max)
  zero_prior <- is.null(dots$prior.mean) || all(dots$prior.mean == 0)
  if (zero_prior) for (attempt in seq_len(30L)) {
    if (is_zero(endpoint)) break
    lambda_max <- 2 * lambda_max
    endpoint <- fit_lambda(lambda_max)
  }
  lambda_nonzero <- NA_real_
  if (zero_prior && !is_zero(preset_fit) && is_zero(endpoint)) {
    lambda_nonzero <- reference
    for (attempt in seq_len(16L)) {
      if (lambda_max/lambda_nonzero <= 1.05) break
      midpoint <- sqrt(lambda_max * lambda_nonzero)
      candidate <- fit_lambda(midpoint)
      if (is_zero(candidate)) { lambda_max <- midpoint; endpoint <- candidate } else lambda_nonzero <- midpoint
    }
  }
  if (!is_zero(endpoint)) warning("lambda_max is not a verified zero-support endpoint.", call. = FALSE)
  r <- log10(reference)
  low <- r + prop_range * (log10(lambda_min) - r)
  high <- r + prop_range * (log10(lambda_max) - r)
  up <- max(2L, ceiling(nlambda/2)); down <- nlambda - up
  grid <- if (nlambda == 2L || low == r) 10^seq(low, high, length.out = nlambda) else
    10^c(seq(low, r, length.out = down + 1L)[seq_len(down)], seq(r, high, length.out = up))
  if (nlambda > 2L || low == r) grid[which.min(abs(log10(grid) - r))] <- reference
  if (prop_range == 1) grid[c(1L, length(grid))] <- c(lambda_min, lambda_max)
  if (anyDuplicated(grid)) stop("prop_range is too small for distinct lambda values.", call. = FALSE)
  cv <- if (do_cv) L0glm_ext_cv(cv_args, grid, foldid, kind, mu_floor, pilot_bounds) else NULL
  if (do_cv && !any(is.finite(cv$mean))) stop("No finite CV score.", call. = FALSE)
  cv_best <- if (do_cv) which.min(cv$mean) else NULL
  lambda_1se <- if (do_cv) max(grid[is.finite(cv$mean) & cv$mean <= cv$mean[cv_best] + cv$se[cv_best]]) else NULL
  rows <- coefficients <- theta <- vector("list", length(grid))
  fit_best <- fit_1se <- NULL
  best_score <- Inf; best <- NA_integer_
  # Keep coefficient paths, not one large design matrix for every lambda.
  for (i in seq_along(grid)) {
    lambda <- grid[i]
    fit <- if (lambda == reference) preset_fit else if (lambda == lambda_max) endpoint else fit_lambda(lambda)
    met <- L0glm.metrics(fit, mu_floor = mu_floor)
    met$support_size_estimated <- met$selected
    if (!is.null(coefs.true)) {
      b <- as.matrix(fit$coefficients); truth <- as.matrix(coefs.true)
      if (!identical(dim(b), dim(truth))) stop("coefs.true must have returned coefficient dimensions.", call. = FALSE)
      if (kind == "multinomial") truth <- truth - rowMeans(truth)
      keep <- setdiff(seq_len(nrow(b)), fit$always.include)
      a <- rowSums(b[keep, , drop = FALSE] != 0) > 0
      t <- rowSums(truth[keep, , drop = FALSE] != 0) > 0
      TP <- sum(a & t); FP <- sum(a & !t); TN <- sum(!a & !t); FN <- sum(!a & t)
      met$F1_score <- if (2*TP + FP + FN > 0) 2*TP/(2*TP + FP + FN) else NA_real_
      met$MCC <- support_mcc(TP, FP, TN, FN)
    }
    # A group-L0 criterion charges selected rows, not each class contrast.
    # Mandatory and estimated nuisance terms are a common additive constant.
    charge <- if (fit$grouped) met$selected else met$nonzero
    score <- if (do_cv) cv$mean[i] else -2 * met$loglik + lambda_preset * (charge + met$mandatory + met$nuisance)
    met$tuning_score <- score
    rows[[i]] <- cbind(lambda = lambda, met)
    coefficients[[i]] <- fit$coefficients; theta[i] <- list(fit$theta)
    if (is.finite(score) && score < best_score) { best <- i; best_score <- score; fit_best <- fit }
    if (do_cv && lambda == lambda_1se) fit_1se <- fit
  }
  metrics <- do.call(rbind, rows)
  if (do_cv) { metrics$cv_deviance <- cv$mean; metrics$cv_se <- cv$se }
  if (is.na(best)) stop("No finite tuning score.", call. = FALSE)
  result <- list(tune = tune, kind = kind, always.include = preset_fit$always.include, p = p,
    start = initialization$start, pilot = initialization$pilot,
    fit_args = args[setdiff(names(args), c("X", "y"))], preprocessing = fit_best$preprocessing,
    metrics = metrics, lambda_best = grid[best], lambda_1se = lambda_1se,
    lambda_preset = lambda_preset, lambda_reference = reference, lambdas = grid,
    lambda_min = lambda_min, lambda_max = lambda_max, lambda_nonzero = lambda_nonzero,
    endpoint_zero = is_zero(endpoint), mu_floor = mu_floor,
    coefficients = coefficients, theta = theta, fit_best = fit_best, fit_1se = fit_1se,
    cv = cv, criterion_support = if (preset_fit$grouped) "rows" else "coefficients")
  if (kind == "multinomial") result$coordinates <- "clr"
  if (printbest) message(sprintf("Best lambda by %s = %.4g", toupper(tune), grid[best]))
  if (plot) L0glm_plot_tuning(result, coefs.true)
  result
}
