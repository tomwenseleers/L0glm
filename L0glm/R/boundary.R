# A one-pass conditional boundary improvement after relaxed nonnegative fitting.
# Signed coefficients (including binomial intercepts) are never candidates.
L0glm_boundary <- function(X, y, beta, eta, weights, family, lower, offset=NULL) {
  original <- list(beta=as.numeric(beta), eta=as.numeric(eta), removed=0L)
  if (!any(lower == 0 & beta > 0)) return(original)
  if (issparse(X)) X <- as(X,"dgCMatrix")
  identity <- family$link == "identity" && family$family %in%
    c("gaussian", "poisson", "quasipoisson")
  if (identity) {
    code <- as.integer(family$family != "gaussian")
    positive <- code == 1L && !is.null(offset) && all(offset>=0) && all(beta>=0) &&
      if (issparse(X)) all(X@x>=0) else all(X>=0)
    answer <- if (issparse(X)) boundary_sparse_cpp(as(X,"dgCMatrix"),
      as.numeric(y),as.numeric(weights),as.numeric(beta),as.numeric(eta),lower,code,positive) else
      boundary_dense_cpp(as.matrix(X),as.numeric(y),as.numeric(weights),
        as.numeric(beta),as.numeric(eta),lower,code,positive)
  } else {
    answer <- original
    for (j in which(lower == 0 & beta > 0)) {
      if (issparse(X)) {
        at <- if (X@p[j] == X@p[j+1L]) integer() else seq.int(X@p[j]+1L,X@p[j+1L])
        rows <- X@i[at]+1L; x <- X@x[at]
      } else {
        x <- as.numeric(X[,j]); rows <- which(x != 0); x <- x[rows]
      }
      if (!length(rows)) {
        answer$beta[j] <- 0; answer$removed <- answer$removed+1L; next
      }
      next_eta <- answer$eta[rows]-x*answer$beta[j]
      change <- tryCatch(suppressWarnings({
        mu <- family$linkinv(next_eta)
        if (any(!is.finite(mu)) || (is.function(family$valideta) && !all(family$valideta(next_eta))) ||
            (is.function(family$validmu) && !all(family$validmu(mu)))) NULL else {
          good <- weights[rows]>0
          derivative <- family$mu.eta(next_eta[good])
          yy <- y[rows][good]; mm <- mu[good]
          score <- (mm-yy)*derivative/family$variance(mm)
          if (family$family %in% c("binomial","quasibinomial")) {
            score[mm==0 & yy==0] <- derivative[mm==0 & yy==0]
            score[mm==1 & yy==1] <- -derivative[mm==1 & yy==1]
          }
          slope <- sum(weights[rows][good]*x[good]*score)
          old_mu <- family$linkinv(answer$eta[rows][good])
          old <- sum(family$dev.resids(yy,old_mu,weights[rows][good]))
          new <- sum(family$dev.resids(yy,mm,weights[rows][good]))
          if (is.finite(slope) && slope>=0 && is.finite(new) &&
              is.finite(old) && new<=old) next_eta else NULL
        }
      }),error=function(e) NULL)
      if (!is.null(change)) {
        answer$beta[j] <- 0;answer$eta[rows] <- change
        answer$removed <- answer$removed+1L
      }
    }
  }
  # Keep the family's own numerical conventions, including custom links.
  if (answer$removed) {
    if (!is.null(offset)) answer$eta <- as.numeric(X%*%answer$beta+offset)
    value <- L0glm_irls_deviance(y,answer$eta,weights,family)
    if (!is.finite(value) || value > L0glm_irls_deviance(y,original$eta,weights,family))
      return(original)
  }
  answer
}

L0glm_ext_boundary <- function(data, result, validation=NULL) {
  if (data$kind == "multinomial" || !any(data$lower == 0)) return(result)
  before <- result$B
  if (data$kind == "nb") {
    eta <- as.matrix(data$Z%*%before+data$offset)
    for (j in seq_len(data$m)) {
      family <- unserialize(serialize(data$family,NULL))
      family$putTheta(log(result$state$theta[j]))
      polished <- L0glm_boundary(data$Z,data$y[,j],before[,j],eta[,j],
        data$weights[,j],family,data$lower,data$offset[,j])
      result$B[,j] <- polished$beta
    }
  } else if (data$kind == "cox") {
    stats <- L0glm_ext_cox_stats(data,before)
    eta <- stats$eta
    for (j in which(data$lower == 0 & before[,1L]>0)) {
      x <- as.numeric(data$Z[,j])
      trial <- tryCatch(L0glm_ext_cox_eta(eta-x*result$B[j,1L],
        data$response,data$weights[,1L],data$cox.cache,sorted=TRUE),error=function(e) NULL)
      if (!is.null(trial) && is.finite(trial$loglik) &&
          -sum(x*trial$score)>=0 && trial$loglik>=stats$loglik) {
        result$B[j,1L] <- 0;eta <- trial$eta;stats <- trial
      }
    }
  }
  if (!is.null(validation) && !identical(before,result$B)) {
    old <- L0glm_validation_loss(validation,L0glm_ext_original_coefficients(data,before),state=result$state)
    new <- L0glm_validation_loss(validation,L0glm_ext_original_coefficients(data,result$B),state=result$state)
    if (!is.finite(new) || new>old) result$B <- before
  }
  result
}
