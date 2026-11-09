# Packed logistic GWAS with six main effects and one pairwise interaction.
# Run from the repository root. Data generation and fitting are timed separately.
benchmark_filebacked_gwas <- function(n=4000L,p=1000000L,train=floor(.8*n),
  screen.nthreads=8L,cache='codex_sandbox/filebacked_gwas',
  output='benchmarks/plots/file-backed') {
  library(L0glm)
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  stopifnot(n>train,train>2L,p>=100L)
  dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  stem <- paste0('gwas-',n,'-',p)
  prefix <- file.path(cache,stem)
  bed <- paste0(prefix,'.bed');metadata <- paste0(prefix,'.rds')
  config <- list(n=as.integer(n),p=as.integer(p),seed=101L)
  if(file.exists(bed)||file.exists(metadata)) {
    stopifnot(file.exists(bed),file.exists(metadata))
    saved <- readRDS(metadata)
    stopifnot(identical(saved$config,config))
  } else {
    Rcpp::sourceCpp('benchmarks/helpers/hwe_genotypes.cpp')
    generation <- system.time({
      maf <- benchmarkHweBed(bed,n,p,101L)
      # BED, BIM and FAM form a portable PLINK 1 dataset.
      write.table(data.frame(1,paste0('snp',seq_len(p)),0,seq_len(p),'C','A'),
        paste0(prefix,'.bim'),row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
      write.table(data.frame(paste0('sample',seq_len(n)),paste0('sample',seq_len(n)),0,0,0,-9),
        paste0(prefix,'.fam'),row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
    })[['elapsed']]
    saved <- list(config=config,maf=maf,generation_seconds=generation)
    saveRDS(saved,metadata)
  }
  X <- l0glmData(bed,n=n,p=p,control=list(block.columns=512L,max.cells=5e7))
  active <- as.integer(round(seq(1,p,length.out=6L)))
  beta <- c(1.1,1.1,-1.1,-1.1,.8,-.8)
  pair <- active[1:2];interaction <- 1.5
  signal <- as.matrix(X[,active,drop=FALSE])
  set.seed(101)
  offset <- rnorm(n,0,.1)
  eta <- -2+drop(signal%*%beta)+interaction*signal[,1]*signal[,2]+offset
  y <- rbinom(n,1,plogis(eta))
  # PLINK's binary phenotype convention is 1=control and 2=case.
  write.table(data.frame(paste0('sample',seq_len(n)),paste0('sample',seq_len(n)),0,0,0,y+1L),
    paste0(prefix,'.fam'),row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
  saved$response <- list(y=y,offset=offset,active=active,beta=beta,
    pair=pair,interaction=interaction,intercept=-2)
  saveRDS(saved,metadata)
  training <- seq_len(train);test <- seq.int(train+1L,n)
  auc <- function(y,prob) {
    positive <- sum(y);negative <- length(y)-positive
    (sum(rank(prob)[y==1])-positive*(positive+1)/2)/(positive*negative)
  }
  records <- list();product_keys <- paste(pair,collapse=':')
  for(mode in c('SIS','ISIS','ISIS + interactions')) {
    extra <- as.integer(mode=='ISIS + interactions')
    gc()
    seconds <- system.time(fit <- l0glmFit(X[training,,drop=FALSE],y[training],
      family=binomial(),intercept=TRUE,offset=offset[training],
      screen=if(mode=='SIS')'SIS'else 'ISIS',screen.nthreads=screen.nthreads,
      interactions=extra))[['elapsed']]
    terms <- fit$interactions$terms
    if(length(terms))product_keys <- union(product_keys,
      vapply(terms,function(z)paste(sort(z),collapse=':'),''))
    pair_row <- if(length(terms))which(vapply(terms,function(z)identical(sort(z),pair),TRUE))else integer()
    truth <- Matrix::sparseMatrix(i=c(1L,active+1L,p+1L+pair_row),
      j=rep(1L,7L+length(pair_row)),x=c(-2,beta,rep(interaction,length(pair_row))),
      dims=dim(coef(fit)))
    score <- l0glmMetrics(fit,beta_true=truth)
    probability <- as.numeric(predict(fit,X[test,,drop=FALSE],offset=offset[test]))
    # Main-effect recovery is always assessed against all six main effects.
    present <- as.numeric(coef(fit)[active+1L,1])!=0
    predicted_pair <- if(length(pair_row))any(coef(fit)[p+1L+pair_row,1]!=0)else FALSE
    row <- cbind(data.frame(method=mode,seconds,interaction_recovered=predicted_pair,
      main_TP=sum(present),main_FN=sum(!present),selected=Matrix::nnzero(coef(fit)[-1,,drop=FALSE]),
      test_AUC=auc(y[test],probability),
      test_logloss=-mean(dbinom(y[test],1,probability,log=TRUE))),
      score[c('gic','TP','FP','FN','TN','F1_score','MCC')])
    records[[mode]] <- row
    saveRDS(list(coefficients=coef(fit),interactions=fit$interactions,screening=fit$screening,
      active=active,beta=beta,pair=pair,interaction=interaction),
      file.path(cache,paste0(stem,'-',extra,'-',sub(' .*','',mode),'-fit.rds')))
  }
  records <- do.call(rbind,records)
  # Use a common truth/candidate universe: all original predictors plus the union
  # of generated products and the true pair, including an ungenerated true pair.
  # An omitted interaction therefore remains a false negative in every model.
  records$TP <- records$main_TP+as.integer(records$interaction_recovered)
  records$FP <- records$selected-records$TP
  records$FN <- length(active)+1L-records$TP
  records$TN <- p+length(product_keys)-records$TP-records$FP-records$FN
  records$F1_score <- 2*records$TP/(2*records$TP+records$FP+records$FN)
  records$MCC <- with(records,(as.double(TP)*TN-as.double(FP)*FN)/
    sqrt(as.double(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
  records$observations <- n;records$training <- train;records$p <- p
  records$bed_GB <- file.info(bed)$size/1e9
  records$generation_seconds <- saved$generation_seconds
  write.csv(records,file.path(output,paste0(stem,'-metrics.csv')),row.names=FALSE)
  print(records[c('method','seconds','gic','TP','FP','FN','F1_score','MCC',
    'interaction_recovered','test_AUC','test_logloss')],row.names=FALSE)
  invisible(records)
}

if(sys.nframe()==0L)benchmark_filebacked_gwas()
