# Pool magnitudes before assigning a sign, so opposite channels cannot cancel.
L0glm_group_start <- function(coefs, nonnegative = FALSE, pooling = "mean_abs") {
  active <- as.numeric(Matrix::rowSums(coefs != 0)) > 0
  average <- as.numeric(Matrix::rowMeans(coefs))
  if (nonnegative) {
    seed <- average
    cancelled <- active & seed == 0
    seed[cancelled] <- as.numeric(Matrix::rowMeans(abs(coefs)))[cancelled]
  } else {
    magnitude <- if (pooling == "rms") sqrt(as.numeric(Matrix::rowMeans(coefs^2))) else
      as.numeric(Matrix::rowMeans(abs(coefs)))
    direction <- sign(average)
    cancelled <- which(active & direction == 0)
    for (j in cancelled) {
      row <- as.numeric(coefs[j, ])
      direction[j] <- sign(row[which.max(abs(row))])
    }
    seed <- direction * magnitude
  }
  list(vars = which(active),
       start = matrix(seed, nrow(coefs), ncol(coefs)))
}
