# A single common seed, built from independent pilot fits. Pool magnitudes
# before choosing a sign, so opposing fits cannot erase an eligible feature.
L0glm_pilot_start <- function(args, bounds) {
  lambdas <- unique(exp(seq(log(bounds[1L]), log(bounds[2L]), length.out = 10L)))
  lambdas[c(1L, length(lambdas))] <- bounds
  args$start <- NULL
  args$relaxed <- FALSE
  total <- absolute <- largest <- direction <- NULL
  timing <- numeric(length(lambdas))
  support <- integer(length(lambdas))
  for (i in seq_along(lambdas)) {
    timing[i] <- system.time(fit <- do.call(L0glm.fit,
      c(args, list(lambda = lambdas[i]))))[["elapsed"]]
    beta <- as.matrix(fit$coefficients)
    if (is.null(total)) {
      total <- absolute <- largest <- direction <- matrix(0, nrow(beta), ncol(beta),
                                                         dimnames = dimnames(beta))
    }
    total <- total + beta
    absolute <- absolute + abs(beta)
    larger <- abs(beta) > largest
    direction[larger] <- sign(beta[larger])
    largest[larger] <- abs(beta[larger])
    selectable <- setdiff(seq_len(nrow(beta)), fit$always.include)
    support[i] <- sum(rowSums(beta[selectable, , drop = FALSE] != 0) > 0)
  }
  direction[total != 0] <- sign(total[total != 0])
  seed <- direction * absolute / length(lambdas)
  list(start = seed, pilot = data.frame(lambda = lambdas, seconds = timing,
                                        support_size = support))
}
