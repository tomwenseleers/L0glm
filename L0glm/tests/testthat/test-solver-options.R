test_that("removed solver choices are rejected consistently", {
  expect_error(control.solver.gen("nnls"), "not available")
  expect_error(lmridge(diag(3), 1:3, solver = "nnls"), "Invalid solver choice")
  for (fun in list(L0glm.fit, L0nnpois))
    expect_error(fun(diag(3), 1:3, solver = "nnls", start = 0,
                     parallelize = FALSE), "Invalid solver choice")
})

