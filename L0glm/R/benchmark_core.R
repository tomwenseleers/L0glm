L0glm_benchmark_methods <- function(method) {
  choices <- c("L0glm", "abess", "L0Learn", "MCP", "SCAD", "lasso",
    "relaxed_lasso", "adaptive_lasso", "adaptive_relaxed_lasso")
  if (!is.character(method) || !length(method) || anyNA(method))
    stop("method must contain method names or 'all'.", call. = FALSE)
  key <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  if (any(key(method) == "all")) return(choices)
  aliases <- c(glmnet = "lasso", l0learn = "L0Learn", alasso = "adaptive_lasso")
  for (i in seq_along(method)) if (key(method[i]) %in% names(aliases)) method[i] <- aliases[[key(method[i])]]
  index <- match(key(method), key(choices))
  if (anyNA(index)) stop("Unknown method: ", paste(method[is.na(index)], collapse = ", "), call. = FALSE)
  unique(choices[index])
}

L0glm_benchmark <- function(sim, family, method, tune, nlambda, lambda, support.size,
  nfolds, foldid, weights, offset, dispersion, intercept, nonnegative, options,
  approximation, adaptive.gamma, adaptive.epsilon, seed, mu_floor, plot, grouped) {
  if (!is.list(sim) || is.null(sim$X) || is.null(sim$y))
    stop("sim must contain X and y.", call. = FALSE)
  X <- sim$X; y <- sim$y
  n <- nrow(X); p <- ncol(X)
  entries <- if (inherits(X, "sparseMatrix")) {
    if ("x" %in% methods::slotNames(X)) X@x else numeric()
  } else X
  if (is.null(n) || is.null(p) || n < 3 || p < 1 || any(!is.finite(entries)))
    stop("X must be a finite numeric design matrix with at least three rows.", call. = FALSE)
  if (length(L0glm_intercepts(X)))
    stop("Supply X without an intercept column and use intercept=TRUE.", call. = FALSE)
  if (is.null(family)) family <- sim$family
  if (is.null(family)) family <- gaussian()
  if (is.character(family)) {
    name <- family
    if (length(name) == 1 && tolower(name) %in% c("cox", "cox.ph", "multinomial")) {
      L0glm_benchmark_require("mgcv")
      family <- if (tolower(name) == "multinomial") mgcv::multinom(K=
        (if (is.matrix(y)) ncol(y) else nlevels(factor(y))) - 1L) else mgcv::cox.ph()
    } else if (length(name) != 1 || !name %in% c("gaussian", "poisson", "binomial", "Gamma", "inverse.gaussian", "quasipoisson", "quasibinomial"))
      stop("Supply a standard GLM family object.", call. = FALSE)
    else family <- get(name, asNamespace("stats"))(link = if (name %in% c("binomial", "quasibinomial")) "logit" else "identity")
  }
  if (is.function(family)) family <- family()
  if (!inherits(family, "family")) stop("Supply a GLM family object.", call.=FALSE)
  kind <- L0glm_ext_kind(family)
  if (identical(kind, "nb")) stop("Estimated-theta NB comparisons are not yet supported.", call.=FALSE)
  if (is.null(kind) && inherits(y, "Surv")) stop("Use family='cox' for Surv responses.", call.=FALSE)
  response <- if (identical(kind, "multinomial")) L0glm_ext_multinomial_response(y, family, n) else
    if (identical(kind, "cox")) L0glm_ext_cox_response(y, n) else NULL
  if (is.null(response)) {
    y <- as.matrix(y); m <- ncol(y)
    if (nrow(y) != n || !is.numeric(y) || any(!is.finite(y)))
      stop("X and y must be finite, aligned numeric observation matrices.", call.=FALSE)
  } else {y <- response$y; m <- response$m}
  if (is.null(kind)) kind <- "glm"
  if (!is.logical(nonnegative) || length(nonnegative) != 1 || is.na(nonnegative))
    stop("nonnegative must be TRUE or FALSE.", call. = FALSE)
  if (kind == "multinomial" && nonnegative) {
    warning("nonnegative is ignored for multinomial CLR models.", call.=FALSE)
    nonnegative <- FALSE
  }
  if (is.null(intercept)) intercept <- kind != "cox" &&
    (!nonnegative || family$family %in% c("binomial", "quasibinomial"))
  if (!is.logical(intercept) || length(intercept) != 1 || is.na(intercept))
    stop("intercept must be TRUE or FALSE.", call. = FALSE)
  if (kind == "cox" && intercept) stop("A Cox intercept is not identifiable; use intercept=FALSE.", call.=FALSE)
  if (!(identical(grouped, "auto") || identical(grouped, TRUE) || identical(grouped, FALSE)))
    stop("grouped must be 'auto', TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(nlambda) || length(nlambda) != 1 || !is.finite(nlambda) || nlambda < 2 || nlambda != floor(nlambda))
    stop("nlambda must be an integer of at least 2.", call. = FALSE)
  tune <- L0glm_ic_name(tune, allow_cv = TRUE)
  approximation <- match.arg(approximation, c("none", "wls", "wls_grouped"))
  if (any(!is.finite(c(adaptive.gamma, adaptive.epsilon, mu_floor))) ||
      length(adaptive.gamma) != 1 || adaptive.gamma <= 0 || length(adaptive.epsilon) != 1 ||
      adaptive.epsilon <= 0 || length(mu_floor) != 1 || mu_floor < 0)
    stop("Invalid adaptive.gamma, adaptive.epsilon or mu_floor.", call. = FALSE)
  if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed) || seed < 0 || seed != floor(seed))
    stop("seed must be a nonnegative integer.", call. = FALSE)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  rng <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (had_seed) assign(".Random.seed", rng, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed)
  methods <- L0glm_benchmark_methods(method)
  all_methods <- any(gsub("[^a-z0-9]", "", tolower(method)) == "all")
  if (is.null(weights)) weights <- sim$weights
  if (is.null(offset)) offset <- sim$offset
  weights <- L0glm_tune_matrix(weights, n, if (kind == "glm") m else 1L, 1, "weights")
  offset <- L0glm_tune_matrix(offset, n, m, 0, "offset")
  if (any(!is.finite(weights) | weights < 0) || any(colSums(weights) <= 0) || any(!is.finite(offset)))
    stop("Scoring weights must be nonnegative with positive channel totals; offsets must be finite.", call. = FALSE)
  if (kind == "multinomial") offset <- offset - rowMeans(offset)
  dispersion <- if (kind == "glm") L0glm_dispersion(dispersion, family, m) else NULL
  mass <- if (kind == "cox") weights[, 1L] * response$status else
    if (kind == "multinomial") weights[, 1L] * rowSums(y) else rowSums(weights)
  if (sum(mass) <= 0) stop("No positive likelihood mass (observed events for Cox).", call.=FALSE)
  if (tune == "cv") foldid <- L0glm_cv_folds(matrix(mass, n, 1L), nfolds, foldid) else
    if (!is.null(foldid)) stop("foldid is only used with tune='cv'.", call. = FALSE)
  if (is.null(support.size)) support.size <- unique(as.integer(round(seq(0,
    min(p, n - 1L - intercept), length.out = nlambda))))
  if (!is.numeric(support.size) || !length(support.size) || anyNA(support.size) ||
      any(!is.finite(support.size) | support.size < 0 | support.size > p | support.size != floor(support.size)))
    stop("support.size must contain integers between zero and ncol(X).", call. = FALSE)
  support.size <- sort(unique(support.size))
  if (!is.list(options)) stop("options must be a named list.", call. = FALSE)
  if (length(options) && !all(vapply(options, is.list, logical(1)))) {
    if (length(methods) != 1L) stop("Use a named backend argument list for each method.", call. = FALSE)
    options <- stats::setNames(list(options), methods)
  }
  if (length(options)) {
    if (is.null(names(options)) || any(!nzchar(names(options)))) stop("Name each options list by method.", call. = FALSE)
    names(options) <- vapply(names(options), function(z) L0glm_benchmark_methods(z)[1L], character(1))
    if (anyDuplicated(names(options))) stop("Duplicated method options.", call. = FALSE)
  }
  validate_grid <- function(z) {
    if (!is.numeric(z) || !length(z) || any(!is.finite(z) | z < 0))
      stop("lambda must contain finite nonnegative numbers.", call. = FALSE)
    sort(unique(z), decreasing = TRUE)
  }
  if (is.list(lambda)) {
    if (is.null(names(lambda)) || any(!nzchar(names(lambda)))) stop("Name method-specific lambda sequences.", call. = FALSE)
    names(lambda) <- vapply(names(lambda), function(z) L0glm_benchmark_methods(z)[1L], character(1))
    if (anyDuplicated(names(lambda))) stop("Duplicated method lambda sequences.", call. = FALSE)
    lambda <- lapply(lambda, validate_grid)
  } else if (!is.null(lambda)) lambda <- validate_grid(lambda)
  truth <- sim$beta_true
  if (!is.null(truth)) {
    if (length(truth) == p * m && intercept) truth <- rbind(0, matrix(truth, p, m))
    if (length(truth) != (p + intercept) * m) stop("beta_true has incompatible dimensions.", call. = FALSE)
    truth <- matrix(as.numeric(truth), p + intercept, m)
    if (kind == "multinomial") truth <- truth - rowMeans(truth)
  }
  data <- list(X=X, design=if (intercept) cbind(Intercept=1, X) else X,
    y=y, m=m, kind=kind, response=response, family=family, weights=weights, offset=offset,
    dispersion=dispersion, intercept=intercept, nonnegative=nonnegative,
    grouped=grouped, nlambda=as.integer(nlambda), support.size=support.size,
    approximation=approximation, adaptive.gamma=adaptive.gamma,
    adaptive.epsilon=adaptive.epsilon, options=options, seed=seed, mu_floor=mu_floor,
    tune=tune)
  if (kind == "cox") data$cox.cache <- L0glm_ext_cox_cache(response, weights[, 1L])
  if (all_methods) {
    supported <- vapply(methods, function(name) tryCatch({
      L0glm_benchmark_setup(name, data)
      TRUE
    }, L0glm_benchmark_unavailable=function(e) FALSE), logical(1))
    methods <- methods[supported]
    if (!length(methods)) stop("No installed method supports the supplied arguments.", call.=FALSE)
  }
  fits <- tuning <- vector("list", length(methods)); names(fits) <- names(tuning) <- methods
  timings <- vector("list", length(methods))
  for (i in seq_along(methods)) {
    name <- methods[i]; warnings <- character(); start <- proc.time()[[3L]]
    path_seconds <- cv_seconds <- 0
    result <- tryCatch(withCallingHandlers({
      grid <- if (name == "abess") support.size else if (is.list(lambda)) lambda[[name]] else lambda
      path_start <- proc.time()[[3L]]
      path <- L0glm_benchmark_path(name, data, grid)
      full_scores <- lapply(path$coefficients, function(B) if (is.null(B)) NULL else
        L0glm.metrics(L0glm_benchmark_fit(B, data), mu_floor=mu_floor))
      path_seconds <- proc.time()[[3L]] - path_start
      cv <- NULL
      if (tune == "cv") {
        cv_start <- proc.time()[[3L]]
        cv <- L0glm_benchmark_cv(name, data, path, foldid)
        cv_seconds <- proc.time()[[3L]] - cv_start
        values <- cv$mean
      } else values <- vapply(full_scores, function(z) if (is.null(z)) Inf else z[[tune]], numeric(1))
      available <- !is.na(values) & values < Inf & !vapply(path$coefficients, is.null, logical(1))
      if (!any(available)) stop("No candidate has an available tuning score; inspect family, grid and predictions.")
      best <- which.min(replace(values, !available, Inf))
      fit <- L0glm_benchmark_fit(path$coefficients[[best]], data)
      candidates <- data.frame(parameter=path$grid, score=values,
        support_size=vapply(full_scores, function(z) if (is.null(z)) NA_real_ else
          if (is.null(z$support_size_estimated)) z$selected else z$support_size_estimated, numeric(1)))
      if (!is.null(cv)) candidates$se <- cv$se
      list(fit=fit, tuning=list(parameter=path$parameter, best=path$grid[best],
        index=best, candidates=candidates, cv=cv, adaptive=path$adaptive,
        scope=path$scope, approximation=path$approximation, backend=path$backend))
    }, warning=function(w) {warnings <<- unique(c(warnings, conditionMessage(w))); invokeRestart("muffleWarning")},
    message=function(w) {warnings <<- unique(c(warnings, conditionMessage(w))); invokeRestart("muffleMessage")}),
    error=function(e) e)
    elapsed <- proc.time()[[3L]] - start
    failed <- inherits(result, "condition")
    if (!failed) {fits[i] <- list(result$fit); tuning[i] <- list(result$tuning)}
    timings[[i]] <- data.frame(method=name, status=if (failed)
      if (inherits(result, "L0glm_benchmark_unavailable")) result$status else "error" else "ok",
      scope=if (failed) NA_character_ else result$tuning$scope,
      approximation=if (failed) NA_character_ else result$tuning$approximation,
      tune=tune, parameter=if (failed) NA_character_ else result$tuning$parameter,
      selected_parameter=if (failed) NA_real_ else result$tuning$best,
      tuning_score=if (failed) NA_real_ else result$tuning$candidates$score[result$tuning$index],
      cv_mean=if (failed || tune != "cv") NA_real_ else result$tuning$cv$mean[result$tuning$index],
      cv_se=if (failed || tune != "cv") NA_real_ else result$tuning$cv$se[result$tuning$index],
      seconds=elapsed, path_seconds=path_seconds, cv_seconds=cv_seconds,
      note=paste(unique(c(if (failed) conditionMessage(result), warnings)), collapse="; "))
  }
  answer <- structure(list(fits=fits, tuning=tuning, timings=do.call(rbind, timings),
    foldid=foldid, beta_true=truth, mu_floor=mu_floor, data=data, call=match.call()),
    class="sparseglm")
  if (any(answer$timings$status != "ok"))
    warning("Some benchmark methods were unavailable or failed; see L0glm.metrics(result)$status and $note.", call.=FALSE)
  if (plot) graphics::plot(answer)
  answer
}

# A common fitted-model wrapper scores every backend on the original data.
L0glm_benchmark_fit <- function(beta, data) {
  X <- data$design
  m <- data$m
  if (length(beta) != ncol(X) * m || any(!is.finite(beta)))
    stop("Backend returned incompatible or nonfinite coefficients.", call.=FALSE)
  beta <- matrix(as.numeric(beta), ncol(X), m,
    dimnames=list(colnames(X), if (data$kind == "multinomial") data$response$levels else
      if (data$kind == "cox") "linear_predictor" else colnames(data$y)))
  if (data$kind == "multinomial") beta <- beta - rowMeans(beta)
  beta <- Matrix::Matrix(beta, sparse=TRUE)
  beta <- methods::as(methods::as(beta, "generalMatrix"), "CsparseMatrix")
  if (data$kind != "glm") {
    result <- structure(list(coefficients=beta, X=X, y=data$y, family=data$family,
      kind=data$kind, weights=data$weights, offset=data$offset, response=data$response,
      theta.estimated=logical(), p=ncol(data$X), always.include=if (data$intercept) 1L else integer(),
      added_intercept=data$intercept, grouped=!identical(data$grouped, FALSE), init.lambda=NA_real_,
      n_ic=if (data$kind == "cox") sum(data$weights[, 1L] * data$response$status) else
        sum(data$weights[, 1L] > 0)), class=c("L0glm_extended_fit", "L0glm_fit"))
    if (data$kind == "cox") result$cox.cache <- data$cox.cache
    return(result)
  }
  scale <- if (is.null(data$dispersion)) rep(1, m) else abs(data$dispersion)
  structure(list(coefficients=beta, X=X, y=data$y, family=data$family,
    prior.weights=data$weights, offset=data$offset, dispersion=scale,
    dispersion.estimated=if (is.null(data$dispersion)) rep(TRUE, m) else data$dispersion < 0,
    p=ncol(data$X), always.include=if (data$intercept) 1L else integer(),
    added_intercept=data$intercept), class="L0glm_fit")
}

L0glm_benchmark_cv <- function(method, data, path, foldid) {
  K <- max(foldid); L <- length(path$grid)
  observation_mass <- if (data$kind == "cox") data$weights[, 1L] * data$response$status else
    if (data$kind == "multinomial") data$weights[, 1L] * rowSums(data$y) else rowSums(data$weights)
  mass <- vapply(seq_len(K), function(k) sum(observation_mass[foldid == k]), numeric(1))
  loss <- matrix(Inf, L, K); adaptive <- vector("list", K)
  grids <- matrix(NA_real_, L, K)
  for (k in seq_len(K)) {
    train <- foldid != k
    d <- data
    d$X <- data$X[train, , drop=FALSE]
    d$design <- data$design[train, , drop=FALSE]
    for (field in c("y", "weights", "offset")) d[[field]] <- data[[field]][train, , drop=FALSE]
    if (data$kind == "cox") {
      d$response <- L0glm_ext_cox_response(d$y, sum(train))
      d$cox.cache <- L0glm_ext_cox_cache(d$response, d$weights[, 1L])
    }
    d$grouped <- identical(path$scope, "grouped")
    grid <- path$grid
    if (method == "L0glm") grid <- grid * (1 - mass[k]/sum(mass))
    grids[, k] <- grid
    fold_path <- L0glm_benchmark_path(method, d, grid)
    adaptive[k] <- list(fold_path$adaptive)
    for (j in seq_len(L)) {
      B <- fold_path$coefficients[[j]]
      if (is.null(B)) next
      fit <- L0glm_benchmark_fit(B, d)
      if (data$kind == "cox") {
        eta <- as.numeric(data$X %*% B) + data$offset[, 1L]
        full <- L0glm_ext_cox_eta(eta, data$response, data$weights[, 1L], data$cox.cache)$loglik
        loss[j, k] <- -2 * (full - as.numeric(stats::logLik(fit)))/mass[k]
      } else loss[j, k] <- L0glm.metrics(fit, newdata=data$X[!train, , drop=FALSE],
          y=data$y[!train, , drop=FALSE], weights=data$weights[!train, , drop=FALSE],
          offset=data$offset[!train, , drop=FALSE], mu_floor=data$mu_floor)$deviance/mass[k]
    }
  }
  proportions <- mass/sum(mass)
  means <- as.numeric(loss %*% proportions)
  se <- rep(NA_real_, L); finite <- is.finite(means)
  centered <- sweep(loss[finite, , drop=FALSE], 1, means[finite], "-")
  se[finite] <- sqrt(as.numeric(centered^2 %*% proportions)/(K - 1))
  list(mean=means, se=se, deviance=loss, foldid=foldid, nfolds=K,
    fold_weights=mass, parameter_train=grids, adaptive=adaptive)
}

L0glm_benchmark_unavailable <- function(message, status="unsupported") {
  stop(structure(list(message=message, call=NULL, status=status),
    class=c("L0glm_benchmark_unavailable", "error", "condition")))
}

L0glm_benchmark_require <- function(package) {
  if (!requireNamespace(package, quietly=TRUE))
    L0glm_benchmark_unavailable(paste0("Install optional package '", package, "' for this method."), "missing_package")
}

L0glm_benchmark_align <- function(path, requested, tolerance=1e-8) {
  if (is.null(requested)) return(path)
  index <- vapply(requested, function(z) {
    distance <- abs(path$grid - z)/pmax(abs(path$grid), abs(z), 1e-300)
    if (!length(distance) || min(distance) > tolerance) NA_integer_ else which.min(distance)
  }, integer(1))
  path$coefficients <- lapply(index, function(i) if (is.na(i)) NULL else path$coefficients[[i]])
  path$grid <- requested
  path
}
