# Check capabilities before fitting, also used to expand method="all".
L0glm_benchmark_setup <- function(method, data) {
  kind <- data$kind
  if (kind != "glm" && method == "L0Learn")
    L0glm_benchmark_unavailable("L0Learn has no native Cox or multinomial loss.")
  if (kind == "multinomial" && method %in% c("MCP", "SCAD"))
    L0glm_benchmark_unavailable("The MCP/SCAD backends have no native joint multinomial loss.")
  gaussian <- data$family$family == "gaussian" && data$family$link == "identity"
  binomial <- data$family$family == "binomial" && data$family$link == "logit"
  poisson <- data$family$family == "poisson" && data$family$link == "log"
  ordinis_native <- gaussian || binomial || poisson
  native <- switch(method, abess=kind != "glm" || gaussian || binomial || poisson,
    L0Learn=gaussian || binomial, MCP=kind == "cox" || ordinis_native,
    SCAD=kind == "cox" || ordinis_native, TRUE)
  approximation <- "none"
  if (method != "L0glm" && data$approximation %in% c("wls", "wls_grouped") &&
      data$family$family == "poisson" && data$family$link == "identity") {
    if (data$approximation == "wls_grouped" && data$m > 1L && !identical(data$grouped, FALSE) &&
        method %in% c("abess", "lasso", "relaxed_lasso", "adaptive_lasso", "adaptive_relaxed_lasso")) {
      row_w <- rowMeans(data$weights)/(rowMeans(data$y) + .1)
      data$weights <- matrix(row_w, nrow(data$X), data$m)
      approximation <- "identity_poisson_row_wls"
    } else {
      data$weights <- data$weights/(data$y + .1)
      approximation <- "identity_poisson_wls"
    }
    data$family <- gaussian()
    gaussian <- TRUE
  } else if (!native) L0glm_benchmark_unavailable(paste(method, "does not natively support this family/link;",
    "identity Poisson can opt into approximation='wls'."))
  m <- data$m
  common <- m == 1L || all(data$weights == data$weights[, 1L])
  if (method == "abess" && approximation == "identity_poisson_wls" && m > 1L &&
      !identical(data$grouped, FALSE) && !common)
    L0glm_benchmark_unavailable(paste("ABESS cannot use channel-specific weights for grouped Poisson WLS;",
      "use approximation='wls_grouped' for shared row weights, or grouped=FALSE for independent fits."))
  can_group <- method == "L0glm" || (kind == "multinomial") || (gaussian && common &&
    method %in% c("abess", "lasso", "relaxed_lasso", "adaptive_lasso", "adaptive_relaxed_lasso"))
  grouped <- m > 1L && !identical(data$grouped, FALSE) && can_group
  if (m > 1L && identical(data$grouped, TRUE) && !can_group)
    L0glm_benchmark_unavailable(paste(method, "has no native grouped adapter for these channel weights and family."))
  scope <- if (m == 1L) "single" else if (grouped) "grouped" else
    if (kind == "multinomial") "ungrouped" else "independent"
  if (kind == "multinomial" && method == "abess" && !grouped)
    L0glm_benchmark_unavailable("Multinomial ABESS selects whole feature rows; it has no CLR-entry selection mode.")
  if (method != "L0glm" && method != "L0Learn" && data$nonnegative && data$intercept &&
      !data$family$family %in% c("binomial", "quasibinomial"))
    L0glm_benchmark_unavailable(paste(method, "cannot constrain its native intercept; use intercept=FALSE or nonnegative=FALSE."))
  if (method == "abess" && !gaussian && any(data$offset != 0))
    L0glm_benchmark_unavailable("ABESS has no offset argument for this non-Gaussian model.")
  if (method == "abess" && !any(data$support.size <= min(ncol(data$X), nrow(data$X) - 1L - data$intercept)))
    L0glm_benchmark_unavailable("No ABESS support.size candidate is feasible for these data.")
  if (kind == "cox" && method %in% c("abess", "MCP", "SCAD")) {
    if (anyDuplicated(data$y[, 1L]))
      L0glm_benchmark_unavailable(paste(method, "Cox tie handling differs from Breslow; this adapter requires untied times."))
    if (any(data$weights != 1) || (method != "abess" &&
        (data$nonnegative || any(data$offset != 0))))
      L0glm_benchmark_unavailable(paste(method, "Cox requires unit weights; MCP/SCAD also require signed coefficients and zero offsets."))
  }
  if (method == "L0Learn" && !gaussian &&
      (!binomial || any(data$weights != 1) || any(data$offset != 0) || any(!data$y %in% c(0, 1))))
    L0glm_benchmark_unavailable("L0Learn logistic paths require binary outcomes, unit weights and zero offsets.")
  backend <- if (method == "L0glm") "L0glm" else
    if (kind == "cox" && method %in% c("MCP", "SCAD")) "ncvreg" else
      switch(method, abess="abess", L0Learn="L0Learn", MCP="ordinis", SCAD="ordinis", "glmnet")
  if (backend == "glmnet" && kind == "glm" && !gaussian && !binomial && !poisson && !data$intercept) {
    mu <- data$family$linkinv(data$offset)
    if (any(!is.finite(mu)) || !data$family$validmu(mu) || any(data$family$variance(mu) <= 0))
      L0glm_benchmark_unavailable("glmnet requires a valid null mean; identity Poisson with a zero offset needs approximation='wls'.")
  }
  if (backend != "L0glm") L0glm_benchmark_require(backend)
  list(data=data, grouped=grouped, joint=grouped || kind == "multinomial",
    scope=scope, approximation=approximation, backend=backend)
}

# Shared path adapters: no backend information criteria are used for selection.
L0glm_benchmark_path <- function(method, data, grid = NULL) {
  original <- data
  setup <- L0glm_benchmark_setup(method, data)
  data <- setup$data
  kind <- data$kind; m <- data$m
  grouped <- setup$grouped; joint <- setup$joint
  options <- data$options[[method]]
  if (is.null(options)) options <- list()
  locked <- c("X", "x", "y", "weights", "weight", "offset", "family", "dispersion",
    "lambda", "lambdaGrid", "support.size", "intercept", "fit.intercept", "nonnegative",
    "lows", "highs", "lower.limits", "upper.limits", "penalty", "relax", "penalty.factor",
    "loss", "excludeFirstK", "maxSuppSize", "nLambda", "nlambda", "tune", "tune.type",
    "tune.path", "foldid", "nfolds", "nFolds", "grouped", "standardize.response", "alpha",
    "type.multinomial", "cox.ties")
  if (length(options) && (is.null(names(options)) || any(!nzchar(names(options))) || anyDuplicated(names(options))))
    stop("Backend options must have unique argument names.", call.=FALSE)
  if (any(names(options) %in% locked))
    stop("Use shared benchmark arguments instead of overriding: ", paste(intersect(names(options), locked), collapse=", "), call.=FALSE)
  adaptive <- NULL
  if (grepl("adaptive", method)) {
    pilot_options <- data$options$L0glm
    if (is.null(pilot_options)) pilot_options <- list()
    fixed <- list(X=original$X, y=original$y, weights=original$weights, offset=original$offset,
      family=original$family, dispersion=original$dispersion, intercept=original$intercept,
      nonnegative=original$nonnegative, start=NULL, lambda="gic", maxit=1,
      grouped=FALSE, relaxed=FALSE, thresh=0, parallelize=FALSE, screen="none")
    pilot <- do.call(L0glm.fit, utils::modifyList(pilot_options, fixed, keep.null=TRUE))
    B <- as.matrix(pilot$coefficients)
    if (data$intercept) B <- B[-1L, , drop=FALSE]
    magnitude <- if (joint) sqrt(rowMeans(B^2)) else abs(B)
    factor <- pmax(magnitude, data$adaptive.epsilon)^(-data$adaptive.gamma)
    if (any(!is.finite(factor))) stop("Adaptive weights overflowed; increase adaptive.epsilon.", call.=FALSE)
    adaptive <- list(coefficients=pilot$coefficients, penalty.factor=factor,
      gamma=data$adaptive.gamma, epsilon=data$adaptive.epsilon)
  }
  if (method == "L0glm") {
    args <- utils::modifyList(list(X=data$X, y=data$y, family=data$family,
      weights=data$weights, offset=data$offset, dispersion=data$dispersion,
      intercept=data$intercept, nonnegative=data$nonnegative, grouped=grouped,
      parallelize=FALSE), options, keep.null=TRUE)
    if (is.null(grid) && !grepl("^quasi", data$family$family)) {
      fit <- do.call(L0glm.tune, c(args, list(tune=if (data$tune == "cv") "gic" else data$tune, nlambda=data$nlambda,
        plot=FALSE, printbest=FALSE, mu_floor=data$mu_floor)))
      path <- list(grid=fit$lambdas, coefficients=fit$coefficients)
    } else {
      if (is.null(grid)) {
        preset <- preset.lambda("gic", max(1, ncol(data$X)), nrow(data$X))
        grid <- exp(seq(log(max(2, preset)*100), log(2), length.out=data$nlambda))
      }
      path <- list(grid=grid, coefficients=lapply(grid, function(value)
        do.call(L0glm.fit, c(args, list(lambda=value)))$coefficients))
    }
    path$parameter <- "lambda"
  } else {
    run <- function(d, g, channel=NULL) L0glm_benchmark_backend(method, d, g, options,
      if (is.null(adaptive)) NULL else if (joint) adaptive$penalty.factor else
        as.matrix(adaptive$penalty.factor)[, channel], grouped)
    if (joint) path <- run(data, grid) else {
      channel_data <- lapply(seq_len(m), function(k) {
        d <- data
        if (kind == "glm") for (field in c("y", "weights", "offset")) d[[field]] <- data[[field]][, k, drop=FALSE]
        d$m <- 1L
        d
      })
      paths <- lapply(seq_len(m), function(k) run(channel_data[[k]], grid, k))
      if (is.null(grid) && m > 1L && method != "abess") {
        positive <- unlist(lapply(paths, `[[`, "grid")); positive <- positive[positive > 0]
        if (!length(positive)) stop("Backend returned no positive path values.")
        grid <- exp(seq(log(max(positive)), log(min(positive)), length.out=data$nlambda))
        grid <- unique(grid)
        paths <- lapply(seq_len(m), function(k) run(channel_data[[k]], grid, k))
      }
      path <- paths[[1L]]
      path$coefficients <- lapply(seq_along(path$grid), function(j) {
        columns <- lapply(paths, function(z) z$coefficients[[j]])
        if (any(vapply(columns, is.null, logical(1)))) NULL else do.call(cbind, columns)
      })
    }
    path$parameter <- if (method == "abess") "support.size" else "lambda"
  }
  path$backend <- setup$backend
  path$scope <- setup$scope; path$approximation <- setup$approximation; path$adaptive <- adaptive
  L0glm_benchmark_align(path, grid)
}

L0glm_benchmark_backend <- function(method, data, grid, options, factor, grouped) {
  X <- data$X; y <- data$y; n <- nrow(X); p <- ncol(X); m <- data$m
  kind <- data$kind
  w <- data$weights[, 1L]; off <- data$offset
  family <- data$family
  gaussian <- family$family == "gaussian" && family$link == "identity"
  binomial <- family$family == "binomial" && family$link == "logit"
  add_intercept <- function(B, a) if (data$intercept) rbind(a, B) else B
  if (method == "abess") {
    sizes <- if (is.null(grid)) data$support.size else grid
    sizes <- sizes[sizes <= min(p, n - 1L - data$intercept)]
    if (!length(sizes)) return(list(grid=numeric(), coefficients=list()))
    yy <- if (gaussian) y - off else y
    weighted_intercept <- gaussian && any(w != 1) && data$intercept
    if (gaussian && any(w != 1)) {
      # Row scaling gives ABESS the exact WLS loss even in versions whose
      # Gaussian selection loss does not consistently apply native weights.
      X <- X * sqrt(w); yy <- yy * sqrt(w); w <- rep(1, n)
      if (weighted_intercept) {
        if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(p))
        X <- cbind(Intercept=sqrt(data$weights[, 1L]), X)
        colnames(X) <- make.unique(colnames(X))
        sizes <- sizes + 1L
        options$always.include <- unique(c(1L, options$always.include + 1L))
        if (!is.null(options$init.active.set)) options$init.active.set <- unique(c(1L, options$init.active.set + 1L))
        if (!is.null(options$group.index)) options$group.index <- c(1L, options$group.index + 1L)
      }
    }
    # ABESS 0.4.11's fallback for unnamed mgaussian responses calls colnames
    # with an invalid do.NULL argument. Supply names without altering the data.
    if (gaussian && m > 1L && is.null(colnames(yy)))
      colnames(yy) <- paste0("channel", seq_len(m))
    if (kind == "multinomial") {
      # Weighted class rows represent multinomial counts exactly without
      # expanding individual trials. Their counts remain together in CV.
      index <- which(y > 0, arr.ind=TRUE)
      X <- X[index[, 1L], , drop=FALSE]
      w <- w[index[, 1L]] * y[index]
      yy <- factor(data$response$levels[index[, 2L]], levels=data$response$levels)
    }
    if (kind == "cox") {
      yy <- as.matrix(y)
    }
    fit <- do.call(abess::abess, utils::modifyList(list(x=X,
      y=if (kind != "glm") yy else if (m == 1L) as.numeric(yy) else yy,
      family=if (kind != "glm") kind else if (m > 1L) "mgaussian" else if (gaussian) "gaussian" else if (binomial) "binomial" else "poisson",
      weight=w, fit.intercept=data$intercept && !weighted_intercept,
      beta.low=if (data$nonnegative) 0 else -.Machine$double.xmax,
      support.size=sizes, tune.path="sequence", tune.type="bic", lambda=0,
      num.threads=1L, seed=data$seed), options))
    coefficients <- lapply(seq_along(fit$support.size), function(j) {
      B <- if (m == 1L) fit$beta[, j, drop=FALSE] else fit$beta[[j]]
      a <- if (m == 1L) fit$intercept[j] else fit$intercept[[j]]
      if (weighted_intercept) {a <- as.numeric(B[1L, ]); B <- B[-1L, , drop=FALSE]}
      if (kind == "multinomial") {B <- cbind(B, 0); a <- c(a, 0)}
      add_intercept(B, a)
    })
    return(L0glm_benchmark_align(list(grid=fit$support.size - as.integer(weighted_intercept),
      coefficients=coefficients), grid))
  }
  if (kind == "cox" && method %in% c("MCP", "SCAD")) {
    args <- list(X=as.matrix(X), y=y, penalty=method, nlambda=data$nlambda)
    if (!is.null(grid)) args$lambda <- grid
    fit <- do.call(ncvreg::ncvsurv, utils::modifyList(args, options))
    return(L0glm_benchmark_align(list(grid=fit$lambda,
      coefficients=lapply(seq_len(ncol(fit$beta)), function(j) fit$beta[, j, drop=FALSE])), grid))
  }
  if (method == "L0Learn") {
    # Explicit mandatory intercept permits exact Gaussian row weighting and
    # nonnegative intercept constraints, which its native intercept cannot impose.
    A <- if (data$intercept) cbind(Intercept=1, X) else X
    yy <- as.numeric(y)
    if (gaussian) {A <- A * sqrt(w); yy <- (yy - as.numeric(off)) * sqrt(w)}
    low <- rep(if (data$nonnegative) 0 else -Inf, ncol(A))
    if (data$intercept && binomial) low[1L] <- -Inf
    fit <- do.call(L0Learn::L0Learn.fit, utils::modifyList(list(x=A, y=yy,
      loss=if (gaussian) "SquaredError" else "Logistic", penalty="L0",
      nLambda=data$nlambda, maxSuppSize=min(ncol(A), max(data$support.size) + data$intercept),
      nGamma=1L, intercept=FALSE, excludeFirstK=as.integer(data$intercept), lows=low,
      lambdaGrid=if (is.null(grid)) list() else list(grid)), options))
    B <- fit$beta[[1L]]
    return(L0glm_benchmark_align(list(grid=fit$lambda[[1L]],
      coefficients=lapply(seq_len(ncol(B)), function(j) B[, j, drop=FALSE])), grid,
      tolerance=5e-6)) # L0Learn rounds returned lambdas to six significant digits.
  }
  if (method %in% c("MCP", "SCAD")) {
    # ordinis currently requires dense input at its documented solver boundary.
    # The released ordinis wrapper cannot accept a family list on recent R;
    # strings select its native Gaussian/logit/Poisson-log implementations.
    fam <- if (gaussian) "gaussian" else if (binomial) "binomial" else "poisson"
    fit <- do.call(getExportedValue("ordinis", "ordinis"), utils::modifyList(list(x=as.matrix(X),
      y=as.numeric(y), family=fam, weights=w, offset=as.numeric(off),
      penalty=tolower(method), lambda=if (is.null(grid)) numeric() else grid,
      nlambda=data$nlambda, intercept=data$intercept,
      lower.limits=rep(if (data$nonnegative) 0 else -Inf, p)), options))
    B <- fit$beta
    coefficients <- lapply(seq_len(ncol(B)), function(j)
      if (data$intercept) B[, j, drop=FALSE] else B[-1L, j, drop=FALSE])
    return(L0glm_benchmark_align(list(grid=fit$lambda, coefficients=coefficients), grid))
  }
  relaxed <- grepl("relaxed", method)
  if (is.null(factor)) factor <- rep(1, p)
  if (p == 1L) {X <- cbind(X, 0); factor <- c(factor, 1)}
  fam <- if (kind != "glm") kind else if (grouped) "mgaussian" else if (gaussian) "gaussian" else if (binomial) "binomial" else
    if (family$family == "poisson" && family$link == "log") "poisson" else family
  yy <- if (gaussian) y - off else y
  if (kind == "cox") {
    L0glm_benchmark_require("survival")
    time <- y[, 1L]
    # glmnet requires positive times. Cox likelihood depends only on ordering
    # and ties, so integer ranks preserve it exactly when zero times occur.
    if (any(time <= 0)) time <- match(time, sort(unique(time)))
    yy <- survival::Surv(time, y[, 2L])
  }
  args <- list(x=X,
    y=if (kind != "glm" || grouped) yy else as.numeric(yy), family=fam, weights=w,
    offset=if (gaussian || !any(off != 0)) NULL else if (kind == "multinomial") off else as.numeric(off), alpha=1,
    lambda=grid, nlambda=data$nlambda, intercept=data$intercept,
    lower.limits=if (data$nonnegative) 0 else -Inf, penalty.factor=factor,
    relax=relaxed, standardize.response=FALSE,
    type.multinomial=if (grouped) "grouped" else "ungrouped", cox.ties="breslow")
  if (kind == "cox") args$intercept <- NULL
  args <- utils::modifyList(args, options, keep.null=TRUE)
  active_refits <- relaxed && gaussian && p > n
  if (active_refits) args$relax <- FALSE
  if (is.null(grid)) {
    ratio <- args$lambda.min.ratio
    if (is.null(ratio)) ratio <- if (n < ncol(X)) .01 else 1e-4
    if (length(ratio) != 1L || !is.finite(ratio) || ratio <= 0 || ratio >= 1)
      stop("lambda.min.ratio must lie strictly between zero and one.", call.=FALSE)
    # Obtain lambda_max in glmnet's own units, including standardization and
    # adaptive penalty factors. The short near-null pilot avoids fitting a
    # second full path. An explicit grid prevents automatic path truncation.
    probe <- args
    probe$relax <- FALSE; probe$nlambda <- 3L; probe$lambda.min.ratio <- .99
    lambda_max <- do.call(glmnet::glmnet, probe)$lambda[1L]
    if (length(lambda_max) != 1L || !is.finite(lambda_max) || lambda_max <= 0)
      stop("glmnet could not determine a finite positive lambda_max.", call.=FALSE)
    lambda_max <- lambda_max * (1 + 1e-6) # keep the null endpoint above roundoff
    args$lambda <- exp(seq(log(lambda_max), log(lambda_max) + log(ratio), length.out=data$nlambda))
  }
  fit <- do.call(glmnet::glmnet, args)
  if (active_refits)
    return(L0glm_benchmark_align(L0glm_benchmark_relax_gaussian(fit, args, p), grid))
  if (relaxed) fit <- fit$relaxed
  if (is.null(fit)) stop("Backend did not return a relaxed path.", call.=FALSE)
  coefficients <- lapply(seq_along(fit$lambda), function(j) {
    multivariate <- is.list(fit$beta)
    B <- if (multivariate) do.call(cbind, lapply(fit$beta, function(z) z[seq_len(p), j, drop=FALSE])) else
      fit$beta[seq_len(p), j, drop=FALSE]
    a <- if (multivariate) fit$a0[, j] else fit$a0[j]
    add_intercept(B, a)
  })
  L0glm_benchmark_align(list(grid=fit$lambda, coefficients=coefficients), grid)
}

# Refit selected columns so wide Gaussian relaxation never copies the entire
# design for each support. The likelihood, constraints and glmnet solver agree
# with its native relaxed fit; the default support limit is n-3 there too.
L0glm_benchmark_relax_gaussian <- function(fit, args, p) {
  maxp <- if (is.null(args$maxp)) nrow(args$x)-3L else args$maxp
  keep <- which(fit$df <= maxp)
  grouped <- is.list(fit$beta)
  m <- if (grouped) length(fit$beta) else 1L
  cache <- new.env(parent=emptyenv())
  coefficients <- lapply(keep, function(j) {
    B <- if (grouped) do.call(cbind, lapply(fit$beta, function(z) z[, j, drop=FALSE])) else
      fit$beta[, j, drop=FALSE]
    S <- which(Matrix::rowSums(abs(B)) > 0)
    key <- paste0("support:", paste(S, collapse=","))
    if (exists(key, cache, inherits=FALSE)) return(get(key, cache, inherits=FALSE))
    if (length(S)) {
      refit <- args
      refit$x <- args$x[, S, drop=FALSE]
      padded <- length(S)==1L
      if (padded) refit$x <- cbind(refit$x, 0)
      for (field in c("penalty.factor", "lower.limits", "upper.limits")) {
        z <- args[[field]]
        if (length(z) == p) z <- z[S]
        if (padded && field=="penalty.factor") z <- c(z, 1)
        refit[[field]] <- z
      }
      refit$lambda <- 0; refit$relax <- FALSE; refit$exclude <- NULL
      z <- do.call(glmnet::glmnet, refit)
      B[,] <- 0
      B[S, ] <- if (grouped) do.call(cbind, lapply(z$beta, function(b) b[seq_along(S), 1L])) else
        z$beta[seq_along(S), 1L]
      a <- if (grouped) z$a0[, 1L] else z$a0[1L]
    } else a <- if (grouped) fit$a0[, j] else fit$a0[j]
    answer <- if (args$intercept) rbind(a, B) else B
    assign(key, answer, cache)
    answer
  })
  list(grid=fit$lambda[keep], coefficients=coefficients)
}
