# Repeated deconvolution experiments for conditional posterior coverage.
# Run from the repository root after installing the current L0glm package.
# "oracle" uses the true predictor support; neither workflow knows coefficients.
benchmark_posterior_coverage <- function(n=1000L,p=n,npeaks=100L,nrep=200L,channels=1L,
  output='benchmarks/6_benchmarks_posterior_coverage_results') {
  library(L0glm)
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  sim <- simul(2,n=n,p=p,n_nonzero=npeaks,seed=2)
  X <- sim$X
  B <- as.numeric(sim$beta_true)%o%if(channels==1L)1 else seq(.5,1.5,length.out=channels)
  active <- which(rowSums(B!=0)>0)
  cases <- list(poisson_high=list(family=poisson(identity),scale=1,sd=0,known=TRUE),
    poisson_moderate=list(family=poisson(identity),scale=.05,sd=0,known=TRUE),
    gaussian_known=list(family=gaussian(),scale=1,sd=5,known=TRUE),
    gaussian_estimated=list(family=gaussian(),scale=1,sd=5,known=FALSE),
    gaussian_moderate=list(family=gaussian(),scale=.05,sd=5,known=FALSE))
  if(channels>1L)cases <- cases[c('poisson_high','gaussian_known','gaussian_estimated')]
  control <- list(nsim=2000L,burnin=500L,seed=41)
  records <- failures <- list()
  for(name in names(cases)) {
    case <- cases[[name]];truth <- B*case$scale;mu <- as.matrix(X%*%truth)
    for(b in seq_len(nrep)) {
      set.seed(if(channels==1L)80000L+b else 81000L+b)
      Y <- matrix(if(case$sd==0)rpois(length(mu),mu)else rnorm(length(mu),mu,case$sd),n,channels)
      for(method in c('selected','oracle')) {
        started <- proc.time()[[3L]]
        result <- tryCatch({
          fit <- l0glmFit(if(method=='oracle')X[,active,drop=FALSE]else X,Y,
            family=case$family,dispersion=if(case$known)if(case$sd==0)1 else case$sd^2 else NULL,
            lambda=if(method=='oracle')0 else 'gic',nonnegative=TRUE,parallelize=FALSE)
          posterior <- summary(fit,inference='posterior',control.inference=control)$inference
          table <- posterior$coefficients
          beta <- if(method=='oracle')truth[active,,drop=FALSE]else truth
          target <- beta[cbind(table$row,table$column)];real <- target!=0
          covered <- target>=table$cred.low & target<=table$cred.high
          data.frame(case=name,replicate=b,method,selected_true=sum(real),true=sum(beta!=0),
            selected_false=sum(!real),covered_true=sum(covered&real),
            coverage_selected=if(any(real))mean(covered[real])else NA_real_,
            coverage_all_true=sum(covered&real)/sum(beta!=0),
            support_refit=isTRUE(fit$support.refit),prior=posterior$prior.source,sampler=posterior$sampler,
            mean_width=if(any(real))mean(table$cred.high[real]-table$cred.low[real])else NA_real_,
            dispersion=mean(fit$dispersion),seconds=proc.time()[[3L]]-started)
        },error=function(e) {
          failures[[length(failures)+1L]]<<-data.frame(case=name,replicate=b,method,error=conditionMessage(e))
          NULL
        })
        if(!is.null(result))records[[length(records)+1L]] <- result
      }
      if(b%%25L==0L)message(name,': ',b,' / ',nrep)
    }
  }
  records <- do.call(rbind,records)
  summary <- do.call(rbind,lapply(split(records,list(records$case,records$method),drop=TRUE),function(x)
    data.frame(case=x$case[1L],method=x$method[1L],replicates=nrow(x),
      coverage_selected=mean(x$coverage_selected,na.rm=TRUE),
      MC_SE=sd(x$coverage_selected,na.rm=TRUE)/sqrt(sum(!is.na(x$coverage_selected))),
      coverage_all_true=mean(x$coverage_all_true),recall=mean(x$selected_true/x$true),
      false_positives=mean(x$selected_false),refit_fraction=mean(x$support_refit),
      mean_width=mean(x$mean_width,na.rm=TRUE),dispersion=mean(x$dispersion),seconds=mean(x$seconds))))
  prefix <- file.path(output,paste0(channels,'-channel'))
  write.csv(records,paste0(prefix,'-replicates.csv'),row.names=FALSE)
  write.csv(summary,paste0(prefix,'-summary.csv'),row.names=FALSE)
  if(length(failures))write.csv(do.call(rbind,failures),paste0(prefix,'-errors.csv'),row.names=FALSE)
  saveRDS(list(n=n,p=p,npeaks=npeaks,nrep=nrep,channels=channels,cases=cases,
    control=control,session=sessionInfo()),paste0(prefix,'-settings.rds'))
  list(summary=summary,replicates=records,failures=failures)
}

if(sys.nframe()==0L) {
  single <- benchmark_posterior_coverage()
  multi <- benchmark_posterior_coverage(n=500L,p=500L,npeaks=50L,nrep=100L,channels=3L)
  print(single$summary,row.names=FALSE)
  print(multi$summary,row.names=FALSE)
}
