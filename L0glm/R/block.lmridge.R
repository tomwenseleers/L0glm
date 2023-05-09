#TODO: add argument to L0glm2 block = FALSE


block.lmridge <- function(X, y, lambdas,
                          coefs = rep(0, ncol(X)),
                          block.ids = NULL,
                          weights = NULL,
                          lower = rep(-.Machine$double.xmax,ncol(X)),
                          upper = rep(.Machine$double.xmax,ncol(X)),
                          nonnegative = FALSE,
                          solver = "osqp",
                          control.block = list(block.size = NULL, maxit = 10, tol = 1E-7),
                          method = NULL,
                          start = rep(0, ncol(X)),
                          control.solver = NULL,
                          verbose = FALSE){
  n <- nrow(X)
  p <- ncol(X)
  
  if(is.null(control.block$block.size) || control.block$block.size > p) control.block$block.size <- p
  if(is.null(block.ids)) block.ids <- ceiling((1:p)/control.block$block.size)
  if(length(lambda) == 1) lambda <- rep(lambda, p)
  if(is.null(control.solver)) control.solver <- do.call("control.solver.gen", list(solver = solver, method = method))
  
  

  Xw <- Xweights(X, weights)
  yw <- Yweights(y, weights)

  interc.ind <- which(grepl(colnames(Xw), pattern = "intercept", ignore.case = TRUE))
  if(max(block.ids) == 1){
    control.block$maxit <- 1
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
  while(iter <= control.block$maxit){
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
                                       weights = weights, lower = lower[block.ids == i], upper = upper[block.ids == i], nonnegative = nonnegative,
                                      block = TRUE, solver = solver, method = method, start = start[block.ids == i],
                                      control.solver = control.solver, verbose = verbose)
      coefs[is.na(coefs) | is.infinite(coefs)] <- 0
    }
    conv <- sqrt(sum((coefs0 - coefs)^2))/sqrt(sum((coefs0)^2)) < control.block$tol || all(coefs == 0)  # ||beta_{k-1} - beta_k||_2 / ||beta_{k-1}||_2
    if(conv) break
    iter <- iter + 1
  }

  if(control.block$maxit > 1 && !conv) warning("Algorithm did not converge, try increasing control.block$block.size or control.block$maxit")
    return(coefs)
}
