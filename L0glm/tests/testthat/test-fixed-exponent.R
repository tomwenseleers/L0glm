test_that("fitting interfaces use a fixed adaptive-ridge exponent", {
  for (fun in list(L0glm.fit, L0nnpois, L0glm.fit.internal)) {
    expect_false("gamma" %in% names(formals(fun)))
    expect_error(fun(diag(2), c(2, 3), gamma = 2), "unused argument")
  }
})

test_that("special-family penalties use inverse squared magnitudes and row norms", {
  B <- matrix(c(2, -3, 0, -4, 5, 0), 3, 2)
  allowed <- matrix(TRUE, 3, 2)
  allowed[2, 2] <- FALSE
  for (grouped in c(FALSE, TRUE)) {
    denominator <- if (grouped) matrix(rowSums(B^2), 3, 2) else B^2
    expected <- 2 / denominator
    expected[!allowed] <- Inf
    expected[1, ] <- 0
    expect_equal(L0glm_ext_penalty(B, allowed, 1L, 8, grouped), expected)
    unpenalized <- L0glm_ext_penalty(B, allowed, 1L, 0, grouped)
    expect_equal(unpenalized[allowed], rep(0, sum(allowed)))
    expect_identical(unpenalized[!allowed], Inf)
  }
})
