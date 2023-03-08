#' Model solvers
#'
#' @import nnls
#' @import Matrix
#' @import glmnet
#' @importFrom osqp osqpSettings osqp solve_osqp
#' @importFrom stats lm.fit
#' @importFrom utils flush.console
#'
#' @description





#' @describeIn L0glm
#'
#' Implementation of different solvers
#' It is called by the upper code \code{L0glm.iwls}.
#'
#' @export
lmridge <- function(X, y, lambda,
                    weights = NULL,
                    lower = rep(-.Machine$double.xmax,ncol(X)),
                    upper = rep(.Machine$double.xmax,ncol(X)),
                    nonnegative = FALSE,
                    solver = "osqp",
                    start = rep(0, ncol(X)),
                    #max iter and tol?
                    verbose){

  if (!is.null(weights)) {
    X <- X * sqrt(weights)
    y <- y * sqrt(weights)
  }
  issparse <- class(X) %in% c("dgCMatrix", "dsCMatrix")

  p <- ncol(X)

  # If nonnegative is true and no lower bound is provided, or the lower bound is larger than 0, it is clipped to 0
  if (nonnegative && all(lower < 0)) lower <- rep(0, p)

  # Row augmentation of the covariate matrix if lambda!=0 to add ridge penalty
  # Ref :
  # https://stats.stackexchange.com/a/164546/27340
  # https://stats.stackexchange.com/a/203699/27340
  if (sum(lambda) != 0) {
    if (issparse) {
    Xaug <- rbind(X, .sparseDiagonal(p, sqrt(lambda)) )
    yaug <- c(y, rep(0, p)) }  else {
      Xaug <- rbind(X, diag(sqrt(lambda),p) )
      yaug <- c(y, rep(0, p))}
  } else {
    Xaug <- X
    yaug <- y
  }

  n <- nrow(Xaug)

  if (solver == "osqp") {
    if (p<=n) {
      XtX = crossprod(Xaug, Xaug)
      Xty = crossprod(Xaug, yaug)
      #Obtain settings  and create model
      settings <- osqpSettings(verbose = FALSE, eps_abs = 1e-5, eps_rel = 1e-5,
                               max_iter = 4000, linsys_solver = 0L)
      if (sum(start)!=0) settings <- c(settings, warm_start = TRUE) else settings <- c(settings, warm_start = FALSE)

      model <- osqp(XtX, -Xty, .sparseDiagonal(p), l = lower, u = upper, pars = settings)

      if (!is.null(start)) model$WarmStart(x=start)

      #Solve model
      coefs <- model$Solve()$x
      coefs <- pmax(lower, pmin(coefs, upper)) # fitted coefficients sometimes go very slightly outside constraint zone due to numerical inaccuracies in solver - this is fixed here via clipping

    } else { # if p > n we use Woodbury identity
      # calculate (X'WX)^-1
      XtWX_inv <- solve(crossprod(X , weights))

      # calculate lambda*I + X * (X'WX)^-1 * X'
      lambdaI_plus_X_XtWX_inv_Xt <- lambda*diag(p) + crossprod(X, XtWX_inv)

      # calculate (lambda*I + X * (X'WX)^-1 * X')^-1
      lambdaI_plus_X_XtWX_inv_Xt_inv <- solve(lambdaI_plus_X_XtWX_inv_Xt)

      # calculate X' * (lambda*I + X * (X'WX)^-1 * X')^-1
      Xt_lambdaI_plus_X_XtWX_inv_Xt_inv <- tcrossprod(X, lambdaI_plus_X_XtWX_inv_Xt_inv)

      # calculate (X'WX + lambda*I)^-1 = (X'WX)^-1 - (X'WX)^-1 * X' * (lambda*I + X * (X'WX)^-1 * X')^-1 * X * (X'WX)^-1
      inv_XtWX_lambdaI <- XtWX_inv - XtWX_inv %*% Xt_lambdaI_plus_X_XtWX_inv_Xt_inv %*% X %*% XtWX_inv

      # calculate crossprod(X,W*y)
      XtWy <- crossprod(X,W*y)

      # set up quadratic programming problem
      pff <- .sparseDiagonal(ncol(X))
      model <- osqp(inv_XtWX_lambdaI, -XtWy, pff, l=lower, u=upper, pars=settings)

      # solve
      coefs = model$Solve()$x
      coefs = pmax(lower, pmin(coefs, upper) ) # fitted coefficients sometimes go very slightly outside constraint zone due to numerical inaccuracies in solver - this is fixed here via clipping
    }
  } else if (solver == "LSCG") {
    if (issparse) { coefs <- solve_sparse_LSCG(A = Xaug, b = yaug, x = start)$coefficients } else
                  { coefs <- solve_dense_LSCG(A = Xaug, b = yaug, x = start)$coefficients }
    coefs <- pmax(lower, pmin(coefs, upper)) # clip coefficients within allowed range
  } else if (solver == "QR") {
    if (issparse) { coefs <- solve_sparse_QR(A = Xaug, b = yaug)$coefficients } else
    { coefs <- solve_dense_QR(A = Xaug, b = yaug)$coefficients }
    coefs <- pmax(lower, pmin(coefs, upper)) # clip coefficients within allowed range
    } else if (solver == "glmnet") {
    model <- bigGlm(x = Xaug, y = yaug,
                    alpha = 0, lambda = 0, intercept = FALSE, standardize = FALSE, standardize.response = FALSE,
                    lower.limits = c(-Inf, lower), upper.limits = c(Inf, upper),
                    family = "gaussian", path = FALSE)
    coefs <- coef(model)[-1]
  } else if (solver == "nnls") {
    coefs <- nnls(A=Xaug, b=yaug)$x
    } else {
    stop("Invalid solver choice") }

  ###Old solvers###
  # Weighted LS update, with optional nonnegativity constraints
  # if(is.null(control.fit$block.size)) control.fit$block.size <- p
  # block.ids <- ceiling((1:p)/control.fit$block.size)
  # beta[!converged_set] <- block.fit(X = Xaug[good,!converged_set,drop=F]*sqrt(W[good]),
  #                                   y = zaug[good]*sqrt(W[good]),
  #                                   block.id = block.ids[!converged_set],
  #                                   nonnegative = nonnegative,
  #                                   control = control.fit)


  # TODO clean this
  # if(F){
  #   if(nonnegative){
  #     beta[!converged_set] <- nnls(A = Xaug[good,!converged_set,drop=F]*sqrt(W[good]),
  #                                  b = zaug[good]*sqrt(W[good]))$x
  #   } else {
  #     beta[!converged_set] <- lm.fit(x = as.matrix(Xaug[good,!converged_set,drop=F]*sqrt(W[good])),
  #                                    y = zaug[good]*sqrt(W[good]),
  #                                    tol = 1E-07)$coefficients
  #     beta[is.na(beta) | is.infinite(beta)] <- 0
  #   }
  # }
  return(coefs)
}
