test_that('source views preserve values, sparsity and screened fits', {
  set.seed(36);n <- 90L;p <- 1200L;X <- matrix(rnorm(n*p),n,p)
  y <- drop(X[,1:3]%*%c(2,-2,2))+rnorm(n)
  source <- l0glmData(X,control=list(block.columns=79L,max.cells=20000))
  expect_identical(as.matrix(source[c(9,2,9),c(4,8),drop=FALSE]),X[c(9,2,9),c(4,8)])
  expect_error(as.matrix(source),'max.cells')
  a <- sisGLM(X,y,nsis=20,nthreads=2)
  b <- sisGLM(source,y,nsis=20,nthreads=2)
  expect_identical(a$stats,b$stats)
  args <- list(y=y,intercept=TRUE,screen='ISIS',nsis=10L,screen.maxit=3L,
    local.search='auto',lambda='gic')
  a <- do.call(l0glmFit,c(list(X=X),args));b <- do.call(l0glmFit,c(list(X=source),args))
  expect_equal(unname(as.matrix(coef(a))),unname(as.matrix(coef(b))),tolerance=1e-11)
  expect_equal(unname(fitted(a)),unname(fitted(b)),tolerance=1e-11)
  expect_s3_class(b$X$data,'L0glm_data')
  sparse <- l0glmData(Matrix::Matrix(X[,1:4],sparse=TRUE))
  expect_s4_class(L0glm_data_matrix(sparse),'sparseMatrix')
})

test_that('imputation uses positive-weight training rows and follows predictor names', {
  set.seed(71);X <- matrix(rnorm(120*1500),120,1500);colnames(X) <- paste0('x',1:1500)
  y <- 2*X[,1]+rnorm(120);X[sample(length(X),2000)] <- NA
  w <- runif(120);w[1:3] <- 0;train <- 1:90;test <- 91:120
  mean <- colMeans(X[train[w[train]>0],,drop=FALSE],na.rm=TRUE)
  impute <- function(A) {na<-which(is.na(A),arr.ind=TRUE);A[na]<-mean[na[,2]];A}
  source <- l0glmData(X,missing='mean',control=list(max.cells=20000,block.columns=83))
  args <- list(y=y[train],weights=w[train],intercept=TRUE,lambda='gic',screen='SIS',
    nsis=10L,local.search=FALSE)
  a <- do.call(l0glmFit,c(list(X=impute(X[train,])),args))
  b <- do.call(l0glmFit,c(list(X=source[train,,drop=FALSE]),args))
  expect_equal(unname(as.matrix(coef(a))),unname(as.matrix(coef(b))),tolerance=1e-10)
  expect_equal(predict(b,X[test,rev(seq_len(ncol(X)))]),predict(a,impute(X[test,])),tolerance=1e-10)
  expect_equal(predict(b,as.data.frame(X[test,])),predict(a,impute(X[test,])),tolerance=1e-10)
  expect_equal(l0glmMetrics(b,X[test,rev(seq_len(ncol(X)))],y[test]),
    l0glmMetrics(a,impute(X[test,]),y[test]),tolerance=1e-9)
  null <- l0glmFit(source,y,weights=w,intercept=TRUE,start=0,screen='none',local.search=FALSE)
  expect_equal(sum(coef(null)!=0),1L)
  expect_equal(dim(coef(null)),c(1501L,1L))
})

test_that('native BED scores match resident scores for every family and threads', {
  set.seed(228);G <- matrix(sample(c(0,1,2,NA),120*19,TRUE),120,19)
  G[,18] <- NA;G[,19] <- G[,1]
  path <- tempfile(fileext='.bed');con <- file(path,'wb')
  writeBin(as.raw(c(0x6c,0x1b,1)),con)
  for(j in 1:19) {
    codes <- c(3L,2L,0L)[G[,j]+1L];codes[is.na(codes)] <- 1L
    writeBin(as.raw(colSums(matrix(codes,4)*4L^(0:3))),con)
  }
  close(con)
  source <- l0glmData(path,n=120L,p=19L)
  expect_equal(L0glm_data_raw(source,1:120,1:19),G)
  if(requireNamespace('pgenlibr',quietly=TRUE)) {
    reference <- pgenlibr::NewPgen(path,raw_sample_ct=120L)
    expect_equal(pgenlibr::ReadIntList(reference,1:19),G)
    pgenlibr::ClosePgen(reference)
  }
  for(case in sis_thread_cases()) {
    args <- case[c('y','family','weights','offset','nuisance')]
    positive <- if(is.null(dim(case$weights)))case$weights>0 else rowSums(case$weights)>0
    mean <- colMeans(G[positive,,drop=FALSE],na.rm=TRUE);mean[is.nan(mean)] <- 0
    A <- G;where <- which(is.na(A),arr.ind=TRUE);A[where] <- mean[where[,2]]
    reference <- do.call(sisGLM,c(list(X=A,nsis=8L,nthreads=1L),args))
    for(threads in c(1L,2L,4L,8L)) {
      fit <- do.call(sisGLM,c(list(X=source,nsis=8L,nthreads=threads),args))
      expect_identical(fit$indices,reference$indices)
      expect_equal(fit$stats,reference$stats,tolerance=1e-11)
    }
  }
  expect_error(l0glmData(path,n=121,p=19),'size')
})

test_that('file-backed and sparse HDF5 adapters support bounded reads', {
  skip_if_not_installed('bigmemory');skip_if_not_installed('HDF5Array')
  set.seed(24);X <- matrix(rnorm(30*8),30,8);folder <- tempfile();dir.create(folder)
  big <- bigmemory::filebacked.big.matrix(30,8,type='double',backingpath=folder,
    backingfile='matrix.bin',descriptorfile='matrix.desc');big[,] <- X
  source <- l0glmData(bigmemory::describe(big))
  expect_equal(as.matrix(source[1:4,c(3,2),drop=FALSE]),X[1:4,c(3,2)])
  h5 <- HDF5Array::writeHDF5Array(X,filepath=file.path(folder,'assay.h5'),name='matrix')
  expect_equal(as.matrix(l0glmData(h5)),X)
  expect_equal(as.matrix(l0glmData(file.path(folder,'assay.h5'))),X)
  loom <- file.path(folder,'assay.loom')
  HDF5Array::writeHDF5Array(X,filepath=loom,name='matrix')
  # Native HDF5 dimensions are genes by cells, as required by the loom schema.
  expect_equal(dim(rhdf5::h5read(loom,'matrix',native=TRUE)),rev(dim(X)))
  expect_equal(as.matrix(l0glmData(loom)),X)
  ad <- file.path(folder,'assay.h5ad');rhdf5::h5createFile(ad)
  rhdf5::h5write(t(X),ad,'X')
  expect_equal(dim(rhdf5::h5read(ad,'X',native=TRUE)),dim(X))
  expect_equal(as.matrix(suppressWarnings(l0glmData(ad))),X)
  for(layout in c('csc','csr')) {
    path <- file.path(folder,paste0(layout,'.h5ad'))
    rhdf5::h5createFile(path);rhdf5::h5createGroup(path,'X')
    A <- as(Matrix::Matrix(if(layout=='csc')X else t(X),sparse=TRUE),'dgCMatrix')
    rhdf5::h5write(A@x,path,'X/data');rhdf5::h5write(A@i,path,'X/indices')
    rhdf5::h5write(A@p,path,'X/indptr')
    file <- rhdf5::H5Fopen(path);group <- rhdf5::H5Gopen(file,'X')
    rhdf5::h5writeAttribute(as.integer(dim(X)),group,'shape')
    rhdf5::h5writeAttribute(paste0(layout,'_matrix'),group,'encoding-type')
    rhdf5::H5Gclose(group);rhdf5::H5Fclose(file)
    source <- suppressWarnings(l0glmData(path))
    expect_equal(as.matrix(source),X)
    expect_s4_class(L0glm_data_matrix(source[,1:3,drop=FALSE]),'sparseMatrix')
  }
})

test_that('Ward clustering accumulates the same cuts and aggregates in row blocks', {
  set.seed(66);n <- 257L;m <- 18L
  mean <- matrix(runif(n*3,1,8),n,3)[,rep(1:3,each=6)]
  for(family in list(poisson(),gaussian(),Gamma('log'),MASS::negative.binomial(3))) {
    Y <- if(family$family=='gaussian')mean+matrix(rnorm(n*m),n,m)else
      if(family$family=='Gamma')matrix(rgamma(n*m,4,scale=mean/4),n,m)else matrix(rpois(n*m,mean),n,m)
    args <- list(family=family,linkage='ward',k.max=4L,weights=runif(n,.5,2))
    source <- l0glmData(Y,control=list(block.cells=600,max.cells=2000))
    a <- do.call(glmCluster,c(list(Y=Y),args));b <- do.call(glmCluster,c(list(Y=source),args))
    expect_identical(a$clusterids,b$clusterids)
    for(name in c('D','history','centres','Y_agg','weights_agg'))
      expect_equal(a[[name]],b[[name]],tolerance=1e-9)
  }
})

test_that('CV predictions use fold-specific genotype means', {
  set.seed(5);n <- 90L;p <- 1000L;X <- matrix(rnorm(n*p),n,p)
  y <- 2*X[,1]+rnorm(n);X[sample(length(X),1500)] <- NA
  w <- runif(n);w[1:3] <- 0;offset <- rnorm(n,0,.1);foldid <- rep(1:3,length.out=n)
  source <- l0glmData(X,missing='mean',control=list(max.cells=8000,block.columns=83))
  tune <- l0glmTune(source,y,family=gaussian(),weights=w,offset=offset,tune='cv',
    foldid=foldid,nlambda=3,screen='SIS',nsis=8L,intercept=TRUE,
    local.search=FALSE,plot=FALSE,printbest=FALSE)
  impute <- function(A,mean) {na<-which(is.na(A),arr.ind=TRUE);A[na]<-mean[na[,2]];A}
  for(k in 1:3) {
    training <- foldid!=k;validation <- !training
    mean <- colMeans(X[training&w>0,,drop=FALSE],na.rm=TRUE)
    A <- impute(X[training,,drop=FALSE],mean);V <- impute(X[validation,,drop=FALSE],mean)
    for(i in seq_along(tune$lambdas)) {
      fit <- l0glmFit(A,y[training],weights=w[training],offset=offset[training],
        screen='SIS',nsis=8L,intercept=TRUE,local.search=FALSE,lambda=tune$cv$lambda_train[i,k])
      loss <- sum(w[validation]*(y[validation]-predict(fit,V,offset=offset[validation]))^2)/sum(w[validation])
      expect_equal(loss,unname(tune$cv$deviance[i,k]),tolerance=1e-9)
    }
  }
})

test_that('PGEN retains dosage precision, missing calls and serialized views', {
  skip_if_not_installed('pgenlibr')
  G <- matrix(c(0,.2,1.1,2,NA,0,1.5,.5,2,1,0,NA),4,3)
  path <- tempfile(fileext='.pgen')
  writer <- pgenlibr::NewPgenWriter(path,4L,3L,dosage_present=TRUE)
  for(j in 1:3)pgenlibr::AppendDosages(writer,G[,j])
  pgenlibr::ClosePgenWriter(writer)
  source <- l0glmData(path)
  expect_equal(L0glm_data_raw(source,1:4,1:3),G,tolerance=1/16384)
  expected <- as.matrix(source[c(3,2,3),c(2,1),drop=FALSE])
  saved <- tempfile();saveRDS(source,saved)
  expect_equal(as.matrix(readRDS(saved)[c(3,2,3),c(2,1),drop=FALSE]),expected)
})

test_that('VCF conversion preserves existing files at a proposed cache prefix', {
  input <- tempfile(fileext='.vcf');writeLines('',input)
  prefix <- tempfile();writeLines('keep',paste0(prefix,'.psam'))
  executable <- tempfile();writeLines('',executable)
  expect_error(L0glm_data_convert_vcf(input,list(cache=prefix,plink=executable,dosage=NULL)),
    'already contains files')
  expect_identical(readLines(paste0(prefix,'.psam')),'keep')
})
test_that("lazy source intercept checks use metadata without a full predictor probe", {
  set.seed(8611)
  X <- matrix(rnorm(10*20),10,20)
  source <- l0glmData(X,control=list(max.cells=10))
  view <- L0glm_design_view(source,intercept=FALSE,names=rep("",20))
  expect_identical(L0glm_intercepts(view),integer())
  expect_identical(L0glm_mandatory(view,3L),3L)
  added <- L0glm_design_view(source,intercept=TRUE,names=rep("",20))
  expect_identical(L0glm_intercepts(added),1L)
  # One complete source row exceeds the guard; a requested column fits.
  expect_error(view[1L,,drop=FALSE],"max.cells")
  expect_length(added[,1L],10L)
})

