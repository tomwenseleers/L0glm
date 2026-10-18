test_that("Cox tied weighted offset fits and survival predictions match coxph", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("survival")
  set.seed(1805)
  X <- matrix(rnorm(300 * 3), 300, 3); off <- runif(300, -.3, .3)
  event <- rexp(300, .1 * exp(X %*% c(.4, -.3, .2) + off)); censor <- rexp(300, .05)
  time <- round(pmin(event, censor), 1); status <- as.numeric(event <= censor)
  y <- survival::Surv(time, status); w <- sample(1:3, 300, TRUE)
  fit <- L0glm.fit(X, y, family = mgcv::cox.ph(), weights = w, offset = off,
    lambda = 0, maxit = 100, relaxed = FALSE, tol = 1e-9)
  ref <- survival::coxph(y ~ X + offset(off), weights = w, ties = "breslow",
    control = survival::coxph.control(eps = 1e-11))
  expect_equal(as.numeric(fit$coefficients), unname(coef(ref)), tolerance = 1e-7)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = 1e-7)
  base <- survival::basehaz(ref, centered = FALSE)
  # basehaz(centered=FALSE) retains the weighted mean of an explicit offset.
  expect_equal(tail(fit$baseline$cumhaz, 1), tail(base$hazard, 1) * exp(-weighted.mean(off, w)),
               tolerance = 1e-6)
  survival <- predict(fit, X[1:5, , drop = FALSE], offset = off[1:5], type = "survival",
                      times = c(0, 1, 10, 1000))
  expect_true(all(survival >= 0 & survival <= 1))
  expect_true(all(apply(survival, 1, function(z) all(diff(z) <= 0))))
  index <- sample.int(300)
  other <- L0glm.fit(Matrix::Matrix(X[index, ], sparse = TRUE), y[index, ],
    family = mgcv::cox.ph(), weights = w[index], offset = off[index],
    lambda = 0, maxit = 100, relaxed = FALSE, tol = 1e-9)
  expect_equal(as.matrix(other$coefficients), as.matrix(fit$coefficients), tolerance = 1e-7)
  expect_error(L0glm.fit(X, survival::Surv(rep(0, 300), time + 1, status),
                         family = mgcv::cox.ph()), "right censoring")
})

test_that("wide block IRLS and Cox kernel solves match coefficient-space quadratics", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("survival")
  set.seed(1807)
  n <- 18; p <- 30
  X <- matrix(rnorm(n * p), n, p)
  y <- factor(rep(1:3, length.out = n))
  args <- list(solver = "eigen", method = 1L, control.solver = control.solver.gen("eigen", 1L))
  data <- L0glm_ext_prepare(X, y, mgcv::multinom(K = 2), NULL, NULL,
    FALSE, NULL, FALSE, "none", 1, NULL, list())
  B <- data$B * .1; D <- matrix(3, p, 3)
  fit <- L0glm_ext_multinomial_working(data, B, D, list(), args)
  Q <- L0glm_ext_contrasts(3)
  eta <- X %*% B; P <- L0glm_ext_softmax(eta)$p; M <- P %*% Q
  H <- matrix(0,p*2,p*2)
  for(j in 1:2) for(k in 1:2) H[seq_len(p)+(j-1)*p,seq_len(p)+(k-1)*p] <-
    crossprod(X, X * (as.numeric(P %*% (Q[,j]*Q[,k])) - M[,j]*M[,k]))
  score <- as.numeric(crossprod(X, (data$y - P) %*% Q))
  ref <- matrix(solve(H + diag(3, p*2), H %*% as.numeric(B %*% Q) + score), p, 2) %*% t(Q)
  expect_equal(as.numeric(fit), as.numeric(ref), tolerance = 1e-8)
  data$Z <- Matrix::Matrix(X, sparse = TRUE)
  sparse <- L0glm_ext_multinomial_working(data, B, D, list(), args)
  expect_equal(sparse, fit, tolerance = 1e-8)
  y <- survival::Surv(rexp(n), rep(c(1, 1, 0), length.out = n))
  data <- L0glm_ext_prepare(X, y, mgcv::cox.ph(), NULL, NULL, FALSE,
                           NULL, FALSE, "none", 1, NULL, list())
  B <- matrix(rnorm(p, sd = .1), p, 1); D <- matrix(2, p, 1)
  proposal <- L0glm_ext_cox_working(data, B, D, list(), args)
  s <- L0glm_ext_cox_stats(data, B); A <- X[s$order, ]
  C <- rev(cumsum(rev(s$d/s$S0^2)))
  W <- diag(s$q) - tcrossprod(s$e) * matrix(C[outer(s$group, s$group, pmax)], n, n)
  H <- crossprod(A, W %*% A); score <- crossprod(A, s$score)
  ref <- solve(H + diag(2, p), H %*% B + score)
  expect_equal(as.numeric(proposal), as.numeric(ref), tolerance = 1e-8)
})
