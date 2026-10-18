test_that("uneven parallel chunks preserve channel inputs through selection and refitting", {
  ncores <- parallel::detectCores()
  skip_if(is.na(ncores) || ncores < 2L)
  set.seed(1823)
  n <- 60L; p <- 8L; m <- 3L
  X <- matrix(rnorm(n * p), n, p)
  truth <- matrix(0, p + 1L, m)
  truth[1, ] <- c(1, -1, 0.5)
  truth[2:4, ] <- rbind(c(2, -2, 1), c(0.7, 0, -0.7), c(0, 1, 0))
  weights <- matrix(runif(n * m, 0.4, 2), n, m)
  offset <- matrix(rnorm(n * m, sd = 0.2), n, m)
  prior <- matrix(rnorm((p + 1L) * m, sd = 0.01), p + 1L, m)
  start <- matrix(1, p + 1L, m)
  start[1, ] <- 0; start[2:3, 2] <- 0; start[7:9, ] <- 0
  y <- cbind(1, X) %*% truth + offset + matrix(rnorm(n * m, sd = 0.3), n, m)
  for (grouped in c(FALSE, TRUE)) {
    args <- list(X = X, y = y, family = gaussian(), weights = weights,
      offset = offset, prior.mean = prior, start = start, intercept = TRUE,
      always.include = 5L, dispersion = c(0.3, -2, -0.5), lambda = "gic",
      grouped = grouped, relaxed = TRUE, group.maxit = 3L, nthreads = 1L,
      openmp = FALSE, parallelize = FALSE)
    serial <- do.call(L0glm.fit, args)
    parallel <- do.call(L0glm.fit, modifyList(args, list(parallelize = TRUE, ncores = 2L)))
    for (name in c("coefficients", "lambdas", "weights", "z", "dispersion"))
      expect_identical(parallel[[name]], serial[[name]])
    expect_equal(dim(coef(serial)), c(p + 1L, m))
    mask <- start == 0; mask[c(1L, 6L), ] <- FALSE
    expect_true(all(as.matrix(coef(serial))[mask] == 0))
    expect_identical(serial$dispersion[1], 0.3)
    expect_true(any(serial$convergence$stage == "relaxed"))
    a <- serial$convergence; b <- parallel$convergence
    rownames(a) <- rownames(b) <- NULL
    expect_identical(a, b)
  }
})

test_that("empty multichannel supports retain their original sparse dimensions", {
  ncores <- parallel::detectCores()
  skip_if(is.na(ncores) || ncores < 2L)
  X <- matrix(seq_len(36), 12, 3)
  y <- matrix(2, 12, 3)
  args <- list(X = X, y = y, family = gaussian(), start = 0,
    grouped = TRUE, relaxed = TRUE, parallelize = FALSE, nthreads = 1L,
    openmp = FALSE)
  serial <- do.call(L0glm.fit, args)
  parallel <- do.call(L0glm.fit, modifyList(args, list(parallelize = TRUE, ncores = 2L)))
  expect_s4_class(coef(serial), "sparseMatrix")
  expect_equal(dim(coef(serial)), c(3L, 3L))
  expect_true(all(coef(serial) == 0))
  for (name in c("coefficients", "lambdas", "weights", "z", "dispersion"))
    expect_identical(parallel[[name]], serial[[name]])
})
