expect_model_criteria <- function(fit) {
  scores <- L0glm.metrics(fit, mu_floor = 0)
  likelihood <- logLik(fit)
  expect_s3_class(fit, "L0glm_fit")
  expect_s3_class(likelihood, "logLik")
  expect_identical(coef(fit), fit$coefficients)
  expect_equal(as.numeric(likelihood), scores$loglik)
  expect_equal(attr(likelihood, "df"), scores$df)
  expect_equal(attr(likelihood, "nobs"), nobs(fit))
  expected <- L0glm:::get.IC(scores$loglik, nobs(fit), max(1, fit$p), scores$df)
  for (ic in c("aic", "bic", "gic", "ebic", "mbic")) {
    value <- switch(ic, aic = AIC(fit), bic = BIC(fit), gic = GIC(fit),
                    ebic = EBIC(fit), mbic = MBIC(fit))
    expect_equal(value, scores[[ic]], info = ic)
    expect_equal(value, expected[[ic]], info = ic)
  }
  expect_false("hbic" %in% names(scores))
}

test_that("Gaussian likelihood methods use prior weights, offsets and profiled scale", {
  set.seed(1901)
  X <- matrix(rnorm(360), 120, 3)
  w <- runif(120, .5, 2); off <- runif(120, -.2, .2)
  y <- 1 + X %*% c(2, -1, .5) + off + rnorm(120, sd = .7/sqrt(w))
  ref <- glm(as.numeric(y) ~ X + offset(off), weights = w, family = gaussian())
  for (design in list(X, Matrix::Matrix(X, sparse = TRUE))) {
    fit <- L0glm.fit(design, y, weights = w, offset = off, intercept = TRUE,
      lambda = 0, thresh = 0, relaxed = FALSE, parallelize = FALSE, maxit = 100)
    expect_model_criteria(fit)
    expect_equal(as.numeric(coef(fit)), unname(coef(ref)), tolerance = 1e-8)
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = 1e-9)
    expect_equal(AIC(fit), AIC(ref), tolerance = 1e-9)
    expect_equal(BIC(fit), BIC(ref), tolerance = 1e-9)
    expect_equal(nobs(fit), 120)
    expect_true(fit$dispersion.estimated)
    expect_equal(as.numeric(fit$prior.weights), w)
    expect_equal(as.numeric(fit$offset), off)
  }
  expect_equal(rownames(AIC(fit, ref)), c("fit", "ref"))
  expect_equal(BIC(fit, ref)$BIC, c(BIC(fit), BIC(ref)))
  expect_equal(AIC(fit, k = 3), -2 * as.numeric(logLik(fit)) + 3 * attr(logLik(fit), "df"))
})

test_that("known and estimated multichannel scales retain their distinct likelihood dimensions", {
  set.seed(1902)
  X <- matrix(rnorm(240), 80, 3)
  w <- matrix(runif(160, .5, 2), 80, 2)
  off <- matrix(runif(160, -.1, .1), 80, 2)
  B <- matrix(c(1, -1, 0, 2, 0, 1), 3, 2)
  y <- X %*% B + off + matrix(rnorm(160), 80, 2)
  fit <- L0glm.fit(X, y, weights = w, offset = off, dispersion = c(2, -1),
    lambda = 0, thresh = 0, grouped = FALSE, relaxed = FALSE, parallelize = FALSE)
  expect_model_criteria(fit)
  mu <- as.matrix(X %*% coef(fit)) + off
  phi <- c(2, sum(w[, 2] * (y[, 2] - mu[, 2])^2)/80)
  expected <- sum(dnorm(y, mu, sqrt(rep(phi, each = 80)/w), log = TRUE))
  expect_equal(as.numeric(logLik(fit)), expected, tolerance = 1e-12)
  expect_equal(attr(logLik(fit), "df"), sum(coef(fit) != 0) + 1)
  expect_identical(fit$dispersion.estimated, c(FALSE, TRUE))
  expect_equal(nobs(fit), 80)
  expect_equal(L0glm.metrics(fit, beta_true = B)$within_row_TP,
    sum((coef(fit) != 0) & (B != 0)))
})

test_that("R and native Poisson methods retain channel-specific prior weights", {
  set.seed(1903)
  X <- matrix(runif(180), 60, 3)
  w <- matrix(runif(120, .5, 2), 60, 2); w[1, ] <- 0
  off <- matrix(runif(120, .5, 1), 60, 2)
  B <- matrix(c(3, 0, 1, 0, 2, 1), 3, 2)
  y <- matrix(rpois(120, X %*% B + off), 60, 2)
  for (design in list(X, Matrix::Matrix(X, sparse = TRUE))) {
    for (native in c(FALSE, TRUE)) {
      fit <- L0glm.fit(design, y, weights = w, offset = off,
        family = poisson(identity), nonnegative = TRUE, openmp = native,
        lambda = 2, grouped = FALSE, relaxed = FALSE, parallelize = FALSE)
      expect_model_criteria(fit)
      mu <- as.matrix(design %*% coef(fit)) + off
      expect_equal(as.numeric(logLik(fit)), sum(w * dpois(y, mu, log = TRUE)), tolerance = 1e-12)
      expect_equal(nobs(fit), 60) # Preserve ordinary-GLM criterion convention.
      expect_equal(as.matrix(fit$prior.weights), w)
      expect_equal(as.matrix(fit$offset), off)
      expect_false(any(fit$dispersion.estimated))
      expect_equal(L0glm.metrics(fit, beta_true = B)$MCC,
        L0glm:::get.metrics(coef(fit), design, y, weights = w, offset = off,
          family = poisson(identity), beta_true = B)$MCC)
    }
  }
  single <- L0nnpois(X, y[, 1], weights = w[, 1], offset = off[, 1],
    intercept = TRUE, parallelize = FALSE)
  expect_model_criteria(single)
  expect_true(is.na(L0glm.metrics(single)$within_row_F1_score))
  expect_equal(L0glm.metrics(single, newdata = X[1:10, ], y = y[1:10, 1],
    weights = w[1:10, 1], offset = off[1:10, 1])$loglik,
    sum(w[1:10, 1] * dpois(y[1:10, 1], as.numeric(cbind(1, X[1:10, ]) %*%
      coef(single)) + off[1:10, 1], log = TRUE)))
  expect_error(L0glm.metrics(single, newdata = X), "Supply y")
  expect_error(L0glm.metrics(single, newdata = X[, 1:2], y = y[, 1]), "dimensions")
  expect_error(L0glm.metrics(single, weights = rep(-1, 60)), "nonnegative")
  expect_error(L0glm.metrics(single, mu_floor = -1), "mu_floor")
  expect_error(L0glm.metrics(list()), "model returned")
})

test_that("binomial and quasi likelihoods follow the family conventions", {
  set.seed(1904)
  X <- matrix(rnorm(300), 100, 3)
  off <- runif(100, -.2, .2); trials <- sample(2:6, 100, TRUE)
  y <- rbinom(100, trials, plogis(.2 + X %*% c(.5, -.4, 0) + off))/trials
  fit <- L0glm.fit(X, y, family = binomial(), weights = trials, offset = off,
    lambda = 0, thresh = 0, parallelize = FALSE, relaxed = FALSE, maxit = 100)
  ref <- glm(y ~ X + offset(off), family = binomial(), weights = trials)
  expect_model_criteria(fit)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = 1e-8)
  expect_equal(AIC(fit), AIC(ref), tolerance = 1e-8)
  quasi <- L0glm.fit(X, y, family = quasibinomial(), weights = trials, offset = off,
    intercept = TRUE, lambda = 0, parallelize = FALSE, relaxed = FALSE, maxit = 100)
  expect_true(is.na(as.numeric(logLik(quasi))))
  for (criterion in list(AIC, BIC, GIC, EBIC, MBIC)) expect_true(is.na(criterion(quasi)))
})

test_that("zero mandatory coefficients and excluded starts keep the correct model dimensions", {
  set.seed(1905)
  X <- matrix(rnorm(120), 40, 3)
  y <- X[, 1] + rnorm(40)
  fit <- L0glm.fit(X, y, dispersion = 1, intercept = TRUE,
    always.include = 2, start = c(0, 1, 0, 0), lambda = 2,
    relaxed = FALSE, parallelize = FALSE)
  expect_model_criteria(fit)
  expect_equal(dim(coef(fit)), c(4L, 1L))
  expect_equal(as.numeric(coef(fit)[4, ]), 0)
  expect_equal(fit$p, 2)
  fit$coefficients[fit$always.include, ] <- 0
  expect_equal(attr(logLik(fit), "df"), sum(coef(fit) != 0) + length(fit$always.include))
})

test_that("continuous-family scoring retains fixed and estimated dispersion conventions", {
  set.seed(1907)
  X <- matrix(runif(180), 60, 3)
  w <- runif(60, .5, 2)
  phi <- .4
  for (family in list(Gamma(log), inverse.gaussian(log))) {
    mu <- exp(.2 + X[, 1])
    y <- if (family$family == "Gamma") rgamma(60, shape = 1/phi, scale = mu*phi) else
      statmod::rinvgauss(60, mean = mu, dispersion = phi)
    for (dispersion in list(phi, NULL)) {
      fit <- L0glm.fit(X, y, family = family, weights = w, dispersion = dispersion,
        intercept = TRUE, lambda = 2, relaxed = FALSE, parallelize = FALSE, maxit = 100)
      expect_model_criteria(fit)
      fitted_mean <- family$linkinv(as.numeric(cbind(1, X) %*% coef(fit)))
      if (!is.null(dispersion)) {
        expected <- if (family$family == "Gamma")
          sum(w * dgamma(y, shape = 1/phi, scale = fitted_mean*phi, log = TRUE)) else
          sum(w * statmod::dinvgauss(y, mean = fitted_mean, dispersion = phi, log = TRUE))
      } else {
        dev <- sum(family$dev.resids(y, fitted_mean, w))
        expected <- -.5 * (family$aic(y, w, fitted_mean, w, dev) - 2)
      }
      expect_equal(as.numeric(logLik(fit)), expected, tolerance = 1e-10)
      expect_equal(attr(logLik(fit), "df"), sum(coef(fit) != 0) + is.null(dispersion))
    }
  }
})

test_that("exact likelihood and optional scoring floors are distinguished", {
  fit <- L0nnpois(matrix(seq(.2, 1, length.out = 10), 10, 1), rep(1, 10),
    lambda = 2, parallelize = FALSE)
  fit$coefficients[,] <- 0
  expect_equal(as.numeric(logLik(fit)), -Inf)
  expect_equal(AIC(fit), Inf)
  scores <- L0glm.metrics(fit, mu_floor = 1e-4)
  likelihood <- logLik(fit, mu_floor = 1e-4)
  expect_true(is.finite(as.numeric(likelihood)))
  expect_equal(as.numeric(likelihood), scores$loglik)
  expect_equal(AIC(likelihood), scores$aic)
  expect_equal(GIC(fit, mu_floor = 1e-4), scores$gic)
  expect_equal(EBIC(fit, mu_floor = 1e-4), scores$ebic)
  expect_equal(MBIC(fit, mu_floor = 1e-4), scores$mbic)
})

test_that("special-family likelihoods and criteria inherit the common model interface", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("survival")
  set.seed(1906)
  X <- matrix(rnorm(240), 80, 3)
  y <- rnbinom(80, mu = exp(.5 + X[, 1]), size = 2)
  for (theta in c(2, -1)) {
    fit <- L0glm.fit(X, y, family = mgcv::nb(theta = theta), intercept = TRUE,
      lambda = 2, relaxed = FALSE, maxit = 100)
    expect_model_criteria(fit)
    expect_equal(attr(logLik(fit), "df"), sum(coef(fit) != 0) + (theta < 0))
  }
  counts <- t(vapply(seq_len(80), function(i) rmultinom(1, 4,
    c(.4, exp(.2 * X[i, 1]), .8)), numeric(3)))
  truth <- rbind(c(.2, -.1, -.1), c(1, -1, 0), c(0, 0, 0), c(0, .3, -.3))
  for (grouped in c(FALSE, TRUE)) {
    fit <- L0glm.fit(X, counts, family = mgcv::multinom(K = 2), intercept = TRUE,
      lambda = 0, grouped = grouped, relaxed = FALSE, maxit = 100)
    expect_model_criteria(fit)
    # A CLR row can have two nonzero entries; its sum constraint leaves one df.
    fit$coefficients[2, ] <- c(1, -1, 0)
    expect_model_criteria(fit)
    selectable <- setdiff(seq_len(nrow(coef(fit))), fit$always.include)
    active <- sum(Matrix::rowSums(abs(coef(fit)[selectable, , drop = FALSE])) > 0)
    nonzero <- sum(coef(fit)[selectable, , drop = FALSE] != 0)
    expected_df <- if (grouped) 2 * active + 2 else nonzero - active + 2
    expect_equal(attr(logLik(fit), "df"), expected_df)
    recovery <- L0glm.metrics(fit, beta_true = truth)
    expect_true(all(c("within_row_F1_score", "within_row_MCC") %in% names(recovery)))
    expect_equal(recovery$loglik, L0glm.metrics(fit)$loglik)
    expect_error(L0glm.metrics(fit, newdata = X), "Supply y")
  }
  event <- rexp(80, exp(X[, 1])); censor <- rexp(80, .4)
  response <- survival::Surv(pmin(event, censor), event <= censor)
  fit <- L0glm.fit(X, response, family = mgcv::cox.ph(), lambda = 2,
    relaxed = FALSE, maxit = 100)
  expect_model_criteria(fit)
  expect_equal(nobs(fit), sum(response[, 2]))
  expect_true(is.na(L0glm.metrics(fit, beta_true = c(1, 0, 0))$within_row_MCC))
})
