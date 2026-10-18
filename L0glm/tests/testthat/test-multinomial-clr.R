test_that("CLR elementwise penalties and masks solve the exact equality quadratic", {
  skip_if_not_installed("mgcv")
  set.seed(1830)
  for (p in c(5L,30L)) {
    n <- 20L; C <- 4L
    X <- matrix(rnorm(n*p),n,p)
    y <- factor(rep(1:C,length.out=n))
    data <- L0glm_ext_prepare(X,y,mgcv::multinom(K=C-1L),NULL,NULL,
      FALSE,NULL,FALSE,"none",NULL,NULL,list())
    B <- matrix(0,p,C); D <- matrix(runif(p*C,.5,3),p,C)
    D[1,] <- 0; D[2,1] <- Inf; D[3,1:3] <- Inf
    data$prior <- matrix(rnorm(p*C,sd=.1),p,C); data$prior <- data$prior-rowMeans(data$prior)
    args <- list(solver="eigen",method=1L,control.solver=control.solver.gen("eigen",1L))
    candidate <- L0glm_ext_multinomial_working(data,B,D,list(),args)
    H <- kronecker(diag(1/C,C)-matrix(1/C^2,C,C),crossprod(X))
    score <- as.numeric(crossprod(X,data$y-1/C))
    finiteD <- replace(D,!is.finite(D),0)
    g <- matrix(H%*%as.numeric(candidate)-score,p,C)+finiteD*(candidate-data$prior)
    for(j in which(rowSums(is.finite(D))>1)) {
      keep <- is.finite(D[j,])
      expect_lt(max(abs(g[j,keep]-mean(g[j,keep]))),1e-8)
    }
    expect_lt(max(abs(rowSums(candidate))),1e-12)
    expect_true(all(candidate[!is.finite(D)]==0))
    expect_true(all(candidate[3,]==0))
    data$Z <- Matrix::Matrix(X*(abs(X)>.4),sparse=TRUE)
    sparse <- L0glm_ext_multinomial_working(data,B,D,list(),args)
    data$Z <- as.matrix(data$Z)
    dense <- L0glm_ext_multinomial_working(data,B,D,list(),args)
    expect_equal(sparse,dense,tolerance=1e-8)
  }
})

test_that("grouped and CLR elementwise L0 are category permutation invariant", {
  skip_if_not_installed("mgcv")
  set.seed(1831)
  n <- 350L; p <- 8L; C <- 3L
  X <- matrix(rnorm(n*p),n,p)
  B <- matrix(0,p+1,C); B[1:4,] <- rbind(c(.3,-.1,-.2),c(.9,-.9,0),c(0,.7,-.7),c(.6,.4,-1))
  P <- L0glm_ext_softmax(cbind(1,X)%*%B)$p
  Y <- t(vapply(seq_len(n),function(i)as.numeric(rmultinom(1,2,P[i,])),numeric(C)))
  colnames(Y) <- c("a","b","c")
  for(grouped in c(FALSE,TRUE)) {
    fit <- L0glm.fit(X,Y,family=mgcv::multinom(K=2),lambda=4,grouped=grouped,
      relaxed=TRUE,normalize="rms",maxit=200L,tol=1e-8)
    for(perm in list(c(3,1,2),c(2,3,1),c(3,2,1),c(2,1,3),c(1,3,2))) {
      other <- L0glm.fit(X,Y[,perm],family=mgcv::multinom(K=2),lambda=4,grouped=grouped,
        relaxed=TRUE,normalize="rms",maxit=200L,tol=1e-8)
      expect_equal(as.matrix(other$coefficients)[,order(perm)],as.matrix(fit$coefficients),tolerance=1e-6)
      expect_equal(other$objective,fit$objective,tolerance=1e-6)
    }
    contrasts <- L0glm_ext_contrasts
    rotation <- matrix(c(cos(.7),sin(.7),-sin(.7),cos(.7)),2,2)
    rotated <- testthat::with_mocked_bindings(
      L0glm.fit(X,Y,family=mgcv::multinom(K=2),lambda=4,grouped=grouped,
        relaxed=TRUE,normalize="rms",maxit=200L,tol=1e-8),
      L0glm_ext_contrasts=function(C) {
        Q <- contrasts(C)
        if(C==3L) Q%*%rotation else Q
      },.package="L0glm")
    expect_equal(rotated$coefficients,fit$coefficients,tolerance=1e-6)
    expect_equal(rotated$objective,fit$objective,tolerance=1e-6)
    expect_lt(max(abs(rowSums(fit$coefficients))),1e-12)
    expect_lt(max(abs(rowSums(predict(fit,type="link")))),1e-12)
    met <- fit$metrics
    expect_equal(met$df,if(grouped)2*(met$selected+1) else met$nonzero-met$selected+2)
    expect_equal(met$objective,-2*met$loglik+4*if(grouped)met$selected else met$nonzero)
    expect_true(all(fit$convergence$converged))
    expect_identical(fit$coordinates,"clr")
  }
})

test_that("binomial multinomial limit, CLR offsets and start values are consistent", {
  skip_if_not_installed("mgcv")
  set.seed(1832)
  X <- matrix(rnorm(250*4),250,4)
  eta <- .2+.7*X[,1]-.5*X[,2]
  y <- rbinom(250,1,plogis(eta)); off <- runif(250,-.2,.2)
  Y <- cbind(1-y,y)
  offset <- cbind(-off/2,off/2)
  fit <- L0glm.fit(X,Y,family=mgcv::multinom(K=1),offset=offset,
    lambda=0,relaxed=FALSE,maxit=100,tol=1e-9)
  ref <- glm(y~X+offset(off),family=binomial(),control=glm.control(epsilon=1e-12))
  expect_equal(as.numeric(fit$coefficients[,2]-fit$coefficients[,1]),unname(coef(ref)),tolerance=1e-7)
  expect_equal(fit$metrics$df,5)
  common <- matrix(runif(250),250,2)
  same <- L0glm.fit(X,Y,family=mgcv::multinom(K=1),offset=offset+common,
    lambda=0,relaxed=FALSE,maxit=100,tol=1e-9,start=matrix(1,5,2))
  expect_equal(fit$coefficients,same$coefficients,tolerance=1e-9)
  expect_equal(rowSums(same$offset),rep(0,250),tolerance=1e-13)
  expect_error(L0glm.fit(X,Y,family=mgcv::multinom(K=1),start=matrix(1,5,1)),"CLR column")
})
