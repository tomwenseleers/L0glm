#TODO: add argument to L0glm2 block = FALSE


block_lmridge <- function(X, y, lambdas,
                          coefs = rep(0, ncol(X)),
                          block.ids = NULL,
                          weights = NULL,
                          pen_weights = rep(1, ncol(X)),
                          lower = rep(-.Machine$double.xmax,ncol(X)),
                          upper = rep(.Machine$double.xmax,ncol(X)),
                          nonnegative = FALSE,
                          solver = "osqp",
                          method = NULL,
                          start = rep(0, ncol(X)),
                          control.solver = NULL,
                          control.fit = list(block.size = NULL, maxit = 10, tol = 1E-7),
                          verbose = FALSE){
  n <- length(y)
  p <- ncol(X)
  if(is.null(control.fit$block.size) || control.fit$block.size > p) control.fit$block.size <- p
  if(is.null(block.ids)) block.ids <- ceiling((1:p)/control.fit$block.size)
  if(length(lambda) == 1) lambda <- rep(lambda, p)

  Xw <- Xweights(X, weights)
  yw <- Yweights(y, weights)

  interc.ind <- which(grepl(colnames(Xw), pattern = "intercept", ignore.case = TRUE))
  if(max(block.ids) == 1){
    control.fit$maxit <- 1
  } else if(length(interc.ind) != 0){ # Rearrange blocks to account for intercepts
    if(nonnegative){ # fit intercepts at the end
      block.ids[interc.ind] <- max(block.ids) + 1
    } else { # start with fitting intercept
      block.ids <- block.ids + 1
      block.ids[interc.ind] <- 1
    }
  }
  blocks <- sort(unique(block.ids))

  iter <- 1
  while(iter <= control.fit$maxit){
    coefs0 <- coefs
    for(i in blocks){
      # Compute residuals
      if(length(blocks) > 1){
        A <- Xw[, block.ids == i,drop=F]
        is.zero <- Matrix::rowSums(A) == 0
        A <- A[!is.zero,,drop=F]
        b <- yw[!is.zero]
        b <- b - Xw[!is.zero, block.ids != i, drop = F] %*% coefs[block.ids != i]
      } else {
        A <- Xw
        b <- yw
      }
      # Compute coefficients
      coefs[block.ids == i] <- lmridge(X = A, y = as.vector(b), lambda = lambda[block.ids == i],
                                       weights = weights, pen_weights = pen_weights[block.ids == i],
                                      lower = lower[block.ids == i], upper = upper[block.ids == i], nonnegative = nonnegative,
                                      block = TRUE, solver = solver, method = method, start = start[block.ids == i],
                                      control.solver = control.solver, verbose = verbose)
      coefs[is.na(coefs) | is.infinite(coefs)] <- 0
    }
    conv <- sqrt(sum((coefs0 - coefs)^2))/sqrt(sum((coefs0)^2)) < control.fit$tol || all(coefs == 0)  # ||beta_{k-1} - beta_k||_2 / ||beta_{k-1}||_2
    if(conv) break
    iter <- iter + 1
  }

  if(control.fit$maxit > 1 && !conv) warning("Algorithm did not converge, maybe increase control.fit$maxit or control.fit$block.size")
    return(coefs)
}
