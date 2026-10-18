test_that("within-row scores pool only correctly recovered nonmandatory rows", {
  truth <- rbind(c(1,0,1), c(0,1,0), c(1,0,0), c(0,0,0), c(1,1,1))
  estimate <- rbind(c(1,1,0), c(0,2,0), c(0,0,0), c(8,0,0), c(5,0,0))
  for (sparse in c(FALSE, TRUE)) {
    score <- function(B) get.metrics(Matrix::Matrix(B, sparse = sparse), diag(5), truth,
      beta_true = truth, always.include = 5L, dispersion = 1)
    result <- score(estimate)
    expect_equal(unlist(result[c("TP","FP","TN","FN")]), c(TP=2,FP=1,TN=0,FN=1))
    expect_equal(unlist(result[paste0("within_row_",c("TP","FP","TN","FN"))]),
      c(within_row_TP=2,within_row_FP=1,within_row_TN=2,within_row_FN=1))
    expect_equal(result$within_row_precision, 2/3)
    expect_equal(result$within_row_recall, 2/3)
    expect_equal(result$within_row_F1_score, 2/3)
    expect_equal(result$within_row_MCC, 1/3)
    altered <- estimate; altered[4:5, ] <- 100
    again <- score(altered)
    expect_equal(again[grep("^within_row_",names(again))],
      result[grep("^within_row_",names(result))])
    expect_equal(result$MAD_coefs, mean(abs(estimate-truth)))
    expect_equal(result$RMSE_coefs, sqrt(mean((estimate-truth)^2)))
    expect_equal(result$R2_coefs,
      summary(lm(as.vector(truth)~0+as.vector(estimate)))$r.squared)
  }
})

test_that("within-row scores handle missing truth, single channels and empty intersections", {
  for (truth in list(NULL, matrix(0,3,2))) {
    result <- get.metrics(matrix(1,3,2),diag(3),matrix(1,3,2),beta_true=truth,dispersion=1)
    expect_true(is.na(result$within_row_F1_score))
    expect_true(is.na(result$within_row_MCC))
    if (!is.null(truth)) expect_equal(result$within_row_TP, 0)
  }
  single <- get.metrics(c(1,0,1),diag(3),matrix(1,3,1),beta_true=c(1,0,0),dispersion=1)
  expect_true(all(is.na(unlist(single[grep("^within_row_",names(single))]))))
  zeros <- get.metrics(matrix(0,3,2),diag(3),matrix(0,3,2),
    beta_true=matrix(0,3,2),dispersion=1)
  expect_equal(zeros$RMSE_coefs,0)
  expect_true(is.na(zeros$R2_coefs))
  expect_false(any(c("MD_selcoefs","MAD_selcoefs","R2_selcoefs","MD_coefs",
    "sMAPE_coefs","bias_coefs","accuracy","specificity",
    "estimated_over_true_support_size") %in% names(zeros)))
  expect_error(get.metrics(c(1,2),diag(3),1:3),"dimensions")
  expect_error(get.metrics(matrix(0,3,2),diag(3),matrix(0,3,2),
    beta_true=matrix(1,2,3)),"dimensions")
})
