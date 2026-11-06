test_that("removed solver choices are rejected consistently", {
  expect_error(controlSolverGen("nnls"), "not available")
  expect_error(lmRidge(diag(3), 1:3, solver = "nnls"), "Invalid solver choice")
  for (fun in list(l0glmFit, l0nnPois))
    expect_error(fun(diag(3), 1:3, solver = "nnls", start = 0,
                     parallelize = FALSE), "Invalid solver choice")
})

