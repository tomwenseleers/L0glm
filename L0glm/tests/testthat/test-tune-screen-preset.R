test_that("a fixed ISIS pool reproduces the preset and stays within each training fold", {
  set.seed(39)
  n <- 60L
  X <- matrix(rnorm(n * 800L), n)
  y <- X[, 5] * 3 - X[, 23] + rnorm(n)
  foldid <- rep(1:3, each=n/3)
  args <- list(X=X, y=y, family=gaussian(), nlambda=5L, prop_range=.2,
    screen.nthreads=1L, plot=FALSE, printbest=FALSE)
  independent <- do.call(l0glmTune, args)
  shared <- do.call(l0glmTune, c(args, list(screen.path="preset")))
  preset <- function(path) path$coefficients[[match(path$lambda_preset, path$lambdas)]]
  expect_identical(preset(shared), preset(independent))
  expect_true(length(shared$screen.candidates) > 0L)
  expect_true(all(vapply(shared$coefficients, function(B)
    all(which(as.numeric(B) != 0) %in% shared$screen.candidates), logical(1))))

  cv <- do.call(l0glmTune, c(args, list(tune="cv", foldid=foldid, screen.path="preset")))
  expect_true(all(is.finite(cv$cv$mean)))
  anchor <- which.min(abs(log(cv$lambdas) - log(cv$lambda_reference)))
  for (i in unique(foldid)) {
    keep <- foldid != i
    training <- l0glmFit(X[keep, , drop=FALSE], y[keep], family=gaussian(),
      lambda=cv$cv$lambda_train[anchor, i], screen.nthreads=1L)
    expect_identical(cv$cv$screen.candidates[[i]], training$screening$candidates)
  }
})
