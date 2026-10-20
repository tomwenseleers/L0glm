# Sparse special families: GIC preset fits, GIC tuning and three-fold CV.
if (requireNamespace("mgcv", quietly=TRUE)) {
  set.seed(204)
  n <- 300L
  X_special <- matrix(rnorm(n*20), n, 20)
  folds_special <- sample(rep(1:3, length.out=n))
  family_par <- par(mgp=c(2, .7, 0))

  # Negative binomial: two nonzero slopes among 20 predictors.
  mu_nb <- exp(.4 + .8*X_special[, 1] - .6*X_special[, 2])
  y_nb <- rnbinom(n, mu=mu_nb, size=2)
  # nb() estimates theta; nb(theta=2) would keep it fixed at two.
  nb_fit <- L0glm.fit(X_special, y_nb, family=mgcv::nb(),
    intercept=TRUE, lambda="gic", relaxed=TRUE, maxit=100, tol=1e-8)
  nb_fit$theta
  coef(nb_fit)
  nb_gic <- L0glm.tune(X_special, y_nb, family=mgcv::nb(),
    intercept=TRUE, tune="gic", nlambda=5,
    relaxed=TRUE, maxit=100, tol=1e-8)
  nb_cv <- L0glm.tune(X_special, y_nb, family=mgcv::nb(),
    intercept=TRUE, tune="cv", nlambda=5, nfolds=3, foldid=folds_special,
    relaxed=TRUE, maxit=100, tol=1e-8)
  c(GIC=nb_gic$lambda_best, CV=nb_cv$lambda_best, CV_1se=nb_cv$lambda_1se)
  predict(nb_cv$fit_1se, X_special[1:5, , drop=FALSE], type="response")

  # Three-category multinomial: only two predictor rows affect probabilities.
  eta <- cbind(0, .9*X_special[, 1] - .7*X_special[, 2],
    -.8*X_special[, 1] + .6*X_special[, 2])
  probabilities <- exp(eta - apply(eta, 1, max))
  probabilities <- probabilities/rowSums(probabilities)
  y_multinomial <- factor(vapply(seq_len(n), function(i)
    sample.int(3, 1, prob=probabilities[i, ]), integer(1)), levels=1:3)
  multinomial_fit <- L0glm.fit(X_special, y_multinomial,
    family=mgcv::multinom(K=2), lambda="gic", grouped=TRUE,
    relaxed=TRUE, maxit=200, tol=1e-8)
  coef(multinomial_fit)               # original-scale, zero-sum CLR rows
  multinomial_gic <- L0glm.tune(X_special, y_multinomial,
    family=mgcv::multinom(K=2), tune="gic", nlambda=5, grouped=TRUE,
    relaxed=TRUE, maxit=200, tol=1e-8)
  multinomial_cv <- L0glm.tune(X_special, y_multinomial,
    family=mgcv::multinom(K=2), tune="cv", nlambda=5,
    nfolds=3, foldid=folds_special, grouped=TRUE,
    relaxed=TRUE, maxit=200, tol=1e-8)
  c(GIC=multinomial_gic$lambda_best, CV=multinomial_cv$lambda_best,
    CV_1se=multinomial_cv$lambda_1se)
  predict(multinomial_cv$fit_1se, X_special[1:5, , drop=FALSE], type="response")

  # Right-censored Cox model: two nonzero log-hazard effects, no intercept.
  if (requireNamespace("survival", quietly=TRUE)) {
    event_time <- rexp(n, rate=exp(.7*X_special[, 1] - .5*X_special[, 2]))
    censor_time <- rexp(n, rate=.4)
    y_cox <- survival::Surv(pmin(event_time, censor_time),
      event_time <= censor_time)
    cox_fit <- L0glm.fit(X_special, y_cox, family=mgcv::cox.ph(),
      lambda="gic", relaxed=TRUE, maxit=100, tol=1e-8)
    coef(cox_fit)
    cox_gic <- L0glm.tune(X_special, y_cox, family=mgcv::cox.ph(),
      tune="gic", nlambda=5, relaxed=TRUE, maxit=100, tol=1e-8)
    cox_cv <- L0glm.tune(X_special, y_cox, family=mgcv::cox.ph(),
      tune="cv", nlambda=5, nfolds=3, foldid=folds_special,
      relaxed=TRUE, maxit=100, tol=1e-8)
    c(GIC=cox_gic$lambda_best, CV=cox_cv$lambda_best, CV_1se=cox_cv$lambda_1se)
    predict(cox_cv$fit_1se, X_special[1:5, , drop=FALSE], type="risk")
  }
  par(family_par)
}
