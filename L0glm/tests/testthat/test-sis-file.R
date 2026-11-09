test_that('mapped backing files preserve every screening context and view', {
  skip_if_not_installed('bigmemory')
  cases <- sis_thread_cases();A <- cases[[1L]]$X
  set.seed(92);rows <- sample.int(nrow(A));columns <- sample.int(ncol(A))
  folder <- tempfile('mapped-');dir.create(folder)
  M <- bigmemory::filebacked.big.matrix(nrow(A)+2L,ncol(A)+2L,type='double',
    backingfile='design.bin',descriptorfile='design.desc',backingpath=folder)
  M[seq_len(nrow(A))+1L,seq_len(ncol(A))+1L] <- A[order(rows),order(columns)]
  view <- bigmemory::sub.big.matrix(M,firstRow=2L,lastRow=nrow(A)+1L,
    firstCol=2L,lastCol=ncol(A)+1L)
  X <- l0glmData(view)[rows,columns,drop=FALSE]
  for(case in cases) {
    context <- L0glm_sis_context(case$y,case$family,nrow(A),case$weights,
      case$offset,case$nuisance)
    for(full in c(TRUE,FALSE)) {
      expected <- sis_score_cpp(A,context,7L,full,1L)
      for(threads in c(1L,8L)) {
        actual <- L0glm_data_score(X,context,7L,full,threads)
        expect_identical(actual[names(actual)!='nthreads'],expected[names(expected)!='nthreads'])
      }
    }
  }
  context <- L0glm_sis_context(rnorm(nrow(A)),gaussian(),nrow(A),NULL,NULL,NULL)
  ids <- c(9L,2L,9L,seq_len(nrow(A)));cols <- c(4L,8L,4L)
  duplicate <- L0glm_sis_context(rnorm(length(ids)),gaussian(),length(ids),NULL,NULL,NULL)
  expected <- sis_score_cpp(A[ids,cols],duplicate,3L,TRUE,1L)
  actual <- L0glm_data_score(X[ids,cols,drop=FALSE],duplicate,3L,TRUE,8L)
  expect_identical(actual[names(actual)!='nthreads'],expected[names(expected)!='nthreads'])
  M[2L,2L] <- NA_real_
  expect_error(L0glm_data_score(X,context,7L,FALSE,8L),'finite')
  M[2L,2L] <- Inf
  expect_error(L0glm_data_score(X,context,7L,FALSE,8L),'finite')
})

test_that('mapped scores retain ties across batches and validate view bounds', {
  skip_if_not_installed('bigmemory')
  folder <- tempfile('mapped-');dir.create(folder)
  set.seed(93);A <- matrix(rep(rnorm(30),17000),30,17000)
  M <- bigmemory::filebacked.big.matrix(30,17000,type='double',
    backingfile='design.bin',descriptorfile='design.desc',backingpath=folder)
  M[,] <- A
  X <- l0glmData(M)
  for(threads in c(1L,8L))
    expect_identical(sisGLM(X,rnorm(30),nsis=11L,nthreads=threads)$indices,seq_len(11L))
  context <- L0glm_sis_context(rnorm(30),gaussian(),30,NULL,NULL,NULL)
  expect_error(sis_bigmemory_score_cpp(M@address,2:31,1:17000,context,7L,FALSE,1L),'range')
  expect_error(sis_bigmemory_score_cpp(M@address,1:30,2:17001,context,7L,FALSE,1L),'range')
})

test_that('other bigmemory types and imputation retain fallback semantics', {
  skip_if_not_installed('bigmemory')
  set.seed(94);A <- matrix(sample(0:2,40*12,TRUE),40,12);y <- rnorm(40)
  for(type in c('double','integer','short','char')) {
    M <- suppressWarnings(bigmemory::as.big.matrix(A,type=type))
    expect_equal(sisGLM(l0glmData(M),y)$stats,sisGLM(A*1,y)$stats)
  }
  folder <- tempfile('mapped-');dir.create(folder)
  M <- bigmemory::filebacked.big.matrix(40,12,type='double',
    backingfile='design.bin',descriptorfile='design.desc',backingpath=folder)
  A <- A*1;A[2,4] <- NA_real_;M[,] <- A
  source <- l0glmData(M,missing='mean')
  expect_equal(sisGLM(source,y)$stats,sisGLM(l0glmData(A,missing='mean'),y)$stats)
  path <- file.path(folder,'source.rds');saveRDS(source,path)
  expect_identical(sisGLM(readRDS(path),y)$indices,sisGLM(source,y)$indices)
})
