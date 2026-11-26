# Shared data generation and one-fit worker for the solver-size benchmark.
solverBenchmarkProblem <- function(case, config) {
  if (case$Experiment == "deconvolution") {
    sim <- L0glm::simul(method=if(case$Noise_family=="poisson")2L else 1L,
      n=case$N,p=case$P,sparse=case$Sparse,seed=config$deconvolution_seed,
      options=list(sd_noise=1))
    return(list(X=sim$X,y=sim$y,beta=as.numeric(sim$beta_true),start=1,
      nonnegative=TRUE,family=if(case$Noise_family=="poisson")poisson(identity)else gaussian()))
  }
  set.seed(config$dense_seed)
  X <- matrix(rnorm(as.double(case$N)*case$P),case$N,case$P)
  active <- sample.int(case$P,min(config$dense_effects,case$P))
  beta <- numeric(case$P);beta[active] <- runif(length(active),-.8,.8)
  set.seed(config$dense_noise_seed)
  y <- as.numeric(X[,active,drop=FALSE] %*% beta[active])+rnorm(case$N)
  list(X=X,y=y,beta=beta,start=NULL,nonnegative=FALSE,family=gaussian())
}

solverBenchmarkFit <- function(case, solver, config) {
  problem <- solverBenchmarkProblem(case,config)
  screen <- if(case$Experiment=="dense_normal"&&case$P>case$N)config$dense_screen else "none"
  warningText <- character();errorText <- NA_character_;fit <- NULL;calls <- 0L
  mode <- if(problem$nonnegative)case$Nonnegative_method else NULL
  runFit <- function() {
    calls <<- calls+1L
    L0glm::l0glmFit(problem$X,problem$y,family=problem$family,lambda=config$criterion,
      nonnegative=problem$nonnegative,start=problem$start,screen=screen,local.search=FALSE,
      maxit=config$maxit,tol=config$tol,thresh=1e-5,ncores=config$ncores,nthreads=config$blas_threads,
      solver=solver$Solver_backend,method=solver$Solver_method,
      control.solver=L0glm::controlSolverGen(solver$Solver_backend,solver$Solver_method,
        nonnegative.method=mode,max.memory=config$max_memory))
  }
  timing <- withCallingHandlers(tryCatch({
    began <- microbenchmark::get_nanotime()
    fit <- runFit()
    microbenchmark::get_nanotime()-began
    },
    error=function(e){errorText <<- conditionMessage(e);NULL}),
    warning=function(w){warningText <<- c(warningText,conditionMessage(w));invokeRestart("muffleWarning")})
  stopifnot(calls==1L)
  values <- as.list(setNames(rep(NA_real_,14L),c("Time","GIC","adj_R2_weighted","R2",
    "RMSE_coefficients","Precision","Recall","F1","MCC","TP","FP","FN","TN","deviance")))
  if(!is.null(timing))values$Time <- timing/1e9
  if(!is.null(fit)) {
    met <- tryCatch(L0glm::l0glmMetrics(fit,beta_true=problem$beta),
      error=function(e){errorText <<- paste("Metric error:",conditionMessage(e));NULL})
    if(!is.null(met)) {
      for(pair in list(c("GIC","gic"),c("adj_R2_weighted","adj_R2_weighted"),
        c("R2","R2_coefs"),c("RMSE_coefficients","RMSE_coefs"),c("Precision","precision"),
        c("Recall","recall"),c("F1","F1_score"),c("MCC","MCC"),c("TP","TP"),
        c("FP","FP"),c("FN","FN"),c("TN","TN")))values[[pair[1L]]] <- met[[pair[2L]]]
      values$deviance <- met$deviance_true/case$N
    }
  }
  convergence <- if(is.null(fit))NULL else fit$convergence
  row <- cbind(case,solver,data.frame(Nonnegative=problem$nonnegative,
    Start=if(is.null(problem$start))"NULL"else"1",Screen=screen,Local_search=FALSE,
    Ncores=config$ncores,BLAS_threads=config$blas_threads,Timing_repetitions=1L,Warmup_fits=0L,
    Fit_calls=calls,Timing_clock="microbenchmark::get_nanotime",Maxit=config$maxit,Tolerance=config$tol,
    Iterations=if(is.null(convergence))NA_integer_ else max(convergence$iterations),
    Converged=if(is.null(convergence))NA else all(convergence$converged),
    Dispersion=if(is.null(fit))NA_real_ else fit$dispersion[1L],
    Status=if(is.null(fit))"fit_error"else if(!is.na(errorText))"metric_error"else"ok",
    Warnings=if(length(warningText))paste(unique(warningText),collapse=" | ")else NA_character_,
    Error=errorText,check.names=FALSE),as.data.frame(values,check.names=FALSE))
  list(row=row,coefficients=if(is.null(fit))NULL else as.numeric(coef(fit)),
    truth=problem$beta,convergence=convergence)
}

solverBenchmarkWorker <- function(helper,case,solver,config) {
  suppressPackageStartupMessages(library(L0glm))
  suppressPackageStartupMessages(requireNamespace("osqp",quietly=TRUE))
  suppressPackageStartupMessages(requireNamespace("glmnet",quietly=TRUE))
  L0glm::l0glmThreads(blas = config$blas_threads, openmp = config$ncores)
  env <- new.env(parent=globalenv());sys.source(helper,envir=env)
  env$solverBenchmarkFit(case,solver,config)
}
