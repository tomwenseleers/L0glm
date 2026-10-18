test_that("unsupported Cox response structures are rejected explicitly", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("survival")
  X <- matrix(seq_len(30), 10, 3)
  y <- survival::Surv(seq_len(10), rep(1, 10))
  attr(y, "strata") <- rep(1:2, 5)
  expect_error(L0glm.fit(X, y, family = mgcv::cox.ph()), "strata")
  delayed <- survival::Surv(rep(0, 10), seq_len(10), rep(1, 10))
  expect_error(L0glm.fit(X, delayed, family = mgcv::cox.ph()), "right censoring")
})
