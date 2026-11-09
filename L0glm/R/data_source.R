#' Memory-bounded access to large predictor and outcome matrices
#' @md
#' @param x A matrix, sparse Matrix, bigmemory object/descriptor, DelayedMatrix,
#'   or path to a BED, PGEN, VCF, H5AD, h5Seurat, loom or HDF5 file.
#' @param format `"auto"` detects objects/extensions. Explicit options are
#'   `"matrix"`, `"bigmemory"`, `"delayed"`, `"bed"`, `"pgen"`, `"vcf"`,
#'   `"h5ad"`, `"h5seurat"`, `"loom"` and `"hdf5"`.
#' @param n,p Sample and variant counts for BED files. PGEN stores its dimensions.
#' @param dataset HDF5 dataset, h5Seurat assay path (default
#'   `"assays/RNA/counts"`), or H5AD layer (`NULL` selects X).
#' @param transpose Whether to transpose the backend. `NULL` uses observations
#'   as rows: H5AD, h5Seurat and loom single-cell assays become cells by genes.
#'   Generic HDF5 and resident matrices retain their supplied orientation.
#' @param missing `"error"` or `"mean"`. BED/PGEN default to mean imputation;
#'   other backends default to error. Means use nonmissing training observations
#'   with positive weight in at least one response. All-missing columns use zero.
#' @param control Named list: `block.columns=256L`, `block.cells=2e6`,
#'   `max.cells=5e7`, `nthreads=1L`, `working.columns=10000L`.
#'   `nthreads` controls selected BED decoding;
#'   `screen.nthreads` in [l0glmFit()] controls fused BED screening.
#'   VCF conversion additionally needs `cache` (new output prefix) and `plink`
#'   (PLINK 2 executable; defaults to `Sys.which("plink2")`). `dosage=NULL`
#'   imports GT hardcalls; `dosage="DS"` preserves expected ALT dosages.
#'   `working.columns` caps the resident candidate pool for large file-backed
#'   Gaussian identity ISIS fits; zero disables that cache. Its size also
#'   respects `max.cells`. Global checks preserve the ordinary ISIS path,
#'   and cached or small-support refits are accepted only when GIC improves.
#' @details
#' Subsets with `drop=FALSE` are cheap matrix views. Explicit `as.matrix()`
#' respects `max.cells`. Sparse backends retain sparse working blocks.
#' [sisGLM()] streams predictor columns, reusing its exact family-specific null
#' context. BED screening decodes two-bit calls directly inside native OpenMP
#' reductions. Ties and within-column summation order match resident matrices.
#' Double-precision bigmemory sources without missing-value imputation are
#' screened directly through their existing native mapping, using read-only access.
#' Row/column views retain their original indexing; the operating system pages
#' the backing file without a design-sized R allocation. Other source types
#' and mean-imputed data use bounded working blocks.
#'
#' Large file-backed Gaussian identity models with an unconstrained coefficient
#' vector, the GIC preset, relaxed refitting and Eigen LLT can also reuse a
#' resident pool of up to `working.columns` predictors when at least three
#' screening rounds are requested. The ordinary global
#' ISIS path is retained. Cached candidates and small-support Gaussian refits
#' are selected only when the original-p GIC improves. Exact repetition of a
#' Gaussian coefficient vector avoids a redundant final file scan.
#' Weights, offsets, compulsory terms and fixed or estimated dispersion retain
#' their meanings. This cache uses scalar or automatic starts; a supplied ridge
#' centre or validation set uses ordinary streamed ISIS. The diagnostics in
#' `fit$screening$working.set` report the pool size and full-file scan count.
#' Set `working.columns=0L` for ordinary streamed ISIS.
#'
#' [l0glmFit()] and [l0glmTune()] fit the retained SIS/ISIS columns with the
#' existing GLM/NB, joint multinomial or Cox engine. Presets and criteria retain
#' the original predictor count, and coefficients retain the original sparse
#' dimensions. Starts, nonnegativity, weights, offsets, interactions and local
#' refinement keep their fitting meanings. A small explicit starting support
#' can also be materialized directly. Intercepts in file-backed designs should
#' be requested with `intercept=TRUE`, named explicitly, or declared mandatory.
#'
#' Cross-validation learns imputation within each training fold. Predictions,
#' validation and bootstrap fits reuse or relearn the corresponding training
#' means. Source files must remain available at their saved paths.
#'
#' The working IRLS design must fit within `max.cells`; unscreened fitting of
#' larger designs requires a smaller active set or a larger explicit guard.
#' This interface streams screening and source access; it does not replace
#' every solver by an arbitrary-size row-streamed IRLS engine.
#'
#' [glmDist()] and Ward [glmCluster()] accumulate distances and cut losses over
#' row blocks. Returned centres and aggregates, and the pairwise channel-distance
#' array, obey the same guard. Bregman agglomeration uses a bounded materialized
#' matrix. Predictor groups retain their documented screening restrictions.
#'
#' BED support is PLINK 1 variant-major, with 0/1/2 copies of BIM allele 1. PGEN uses
#' the official pgenlibr dosage reader, with ALT alleles combined for multiallelic
#' records. VCF conversion creates a reusable PGEN cache with input provenance;
#' existing unrelated caches are preserved. HDF5Array handles dense and sparse
#' H5AD layers and delayed assays; h5Seurat CSC assays use its sparse seed.
#' The physical HDF5 layout determines how much I/O a column subset needs.
#' @return A serializable `L0glm_data` matrix descriptor. It retains paths or
#'   backend references and views, with bounded reads on demand.
#' @seealso [l0glmFit()], [l0glmTune()], [sisGLM()], [glmCluster()]
#'   and the [file-backed tutorial](https://tomwenseleers.github.io/L0glm/articles/file-backed.html)
#'   for reproducible million-predictor Gaussian and packed-genotype examples.
#' @examples
#' set.seed(18)
#' X <- matrix(rnorm(100*2000),100,2000)
#' y <- 2*X[,1]-2*X[,2]+rnorm(100)
#' source <- l0glmData(X,control=list(block.cells=20000,max.cells=30000))
#' fit <- l0glmFit(source,y,screen="SIS",nsis=20L,lambda="gic")
#' head(predict(fit,source[1:4,,drop=FALSE]))
#' \dontrun{
#' genotype <- l0glmData("genotypes.bed",n=4000L,p=1733480L)
#' fit <- l0glmFit(genotype,y,family=binomial(),screen="SIS",nsis=1000L)
#' assay <- l0glmData("cells.h5ad",dataset="counts")
#' fit <- l0glmFit(assay,response,family=poisson())
#' }
#' @export
l0glmData <- function(x,format='auto',n=NULL,p=NULL,dataset=NULL,
                     transpose=NULL,missing=NULL,control=list()) {
  defaults <- list(block.columns=256L,block.cells=2e6,max.cells=5e7,nthreads=1L,working.columns=10000L,
    plink=NULL,cache=NULL,dosage=NULL)
  if(!is.list(control) || length(control) && (is.null(names(control)) ||
      any(!names(control)%in%names(defaults)) || anyDuplicated(names(control))))
    stop('Use named data controls.',call.=FALSE)
  defaults[names(control)] <- control
  if(any(!vapply(defaults[1:4],function(z)is.numeric(z)&&length(z)==1L&&is.finite(z)&&z>=1,TRUE)))
    stop('Data block limits and threads must be positive.',call.=FALSE)
  value <- defaults$working.columns
  if(!is.numeric(value)||length(value)!=1L||!is.finite(value)||value<0||
     value!=floor(value)||value>.Machine$integer.max)
    stop('working.columns must be a nonnegative integer.',call.=FALSE)
  if(inherits(x,'L0glm_data'))return(x)
  if(format=='auto')format <- if(is.character(x))if(grepl('\\.vcf(\\.gz)?$',x,ignore.case=TRUE))'vcf'else
    switch(tolower(tools::file_ext(x)),bed='bed',pgen='pgen',h5ad='h5ad',loom='loom',
      h5seurat='h5seurat',h5='hdf5',stop('Specify the data format.',call.=FALSE))else
    if(inherits(x,'big.matrix')||inherits(x,'big.matrix.descriptor'))'bigmemory'else
    if(inherits(x,'DelayedMatrix'))'delayed'else 'matrix'
  object <- x
  if(format=='vcf') {x <- L0glm_data_convert_vcf(x,defaults);object <- x;format <- 'pgen'}
  if(is.null(transpose))transpose <- format%in%c('h5ad','h5seurat')
  if(!is.logical(transpose)||length(transpose)!=1L||is.na(transpose))stop('transpose must be TRUE or FALSE.',call.=FALSE)
  if(format=='bed') {
    if(length(x)!=1L||!file.exists(x)||is.null(n)||is.null(p))
      stop('BED input needs a path and sample/predictor counts n and p.',call.=FALSE)
    object <- normalizePath(x,winslash='/',mustWork=TRUE)
    dimensions <- as.integer(c(n,p));labels <- list(NULL,NULL);sparse <- FALSE
    l0glm_bed_read_cpp(object,dimensions[1],dimensions[2],integer(),integer(),1L)
  } else if(format=='pgen') {
    if(!requireNamespace('pgenlibr',quietly=TRUE))stop('Install pgenlibr for PGEN input.',call.=FALSE)
    object <- normalizePath(x,winslash='/',mustWork=TRUE)
    handle <- pgenlibr::NewPgen(object)
    dimensions <- c(pgenlibr::GetRawSampleCt(handle),pgenlibr::GetVariantCt(handle))
    pgenlibr::ClosePgen(handle);labels <- list(NULL,NULL);sparse <- FALSE
    if(!is.null(n)&&n!=dimensions[1]||!is.null(p)&&p!=dimensions[2])stop('PGEN dimensions differ from n or p.',call.=FALSE)
  } else {
    if(format=='bigmemory') {
      if(!requireNamespace('bigmemory',quietly=TRUE))stop('Install bigmemory.',call.=FALSE)
      if(inherits(x,'big.matrix.descriptor'))object <- bigmemory::attach.big.matrix(x)
    } else if(format%in%c('h5ad','hdf5','loom','h5seurat')) {
      if(!requireNamespace('HDF5Array',quietly=TRUE))stop('Install HDF5Array.',call.=FALSE)
      object <- if(format=='h5ad')HDF5Array::H5ADMatrix(x,layer=dataset)else {
        if(is.null(dataset))dataset <- if(format=='h5seurat')'assays/RNA/counts'else 'matrix'
        if(format=='h5seurat')L0glm_data_seurat(x,dataset)else HDF5Array::HDF5Array(x,dataset)
      }
    } else if(!format%in%c('matrix','delayed'))stop('Unsupported data format.',call.=FALSE)
    if(format=='matrix' && !(is.matrix(object)||inherits(object,'Matrix')))
      stop('Use a numeric matrix or a supported matrix backend.',call.=FALSE)
    if(transpose)object <- t(object)
    dimensions <- dim(object);labels <- dimnames(object)
    sparse <- inherits(object,'sparseMatrix') || inherits(object,'DelayedMatrix') &&
      DelayedArray::is_sparse(object)
  }
  if(length(dimensions)!=2L||anyNA(dimensions)||any(dimensions<1L))
    stop('Input dimensions must be positive.',call.=FALSE)
  if(is.null(missing))missing <- if(format%in%c('bed','pgen'))'mean'else 'error'
  missing <- match.arg(missing,c('error','mean'))
  cache <- new.env(parent=emptyenv());cache$means <- rep(NA_real_,dimensions[2])
  cache$rows <- seq_len(dimensions[1]);cache$prediction <- FALSE
  result <- structure(list(object=object,format=format,dimensions=as.integer(dimensions),
    labels=labels,rows=seq_len(dimensions[1]),columns=seq_len(dimensions[2]),
    sparse=sparse,missing=missing,control=defaults,cache=cache),class='L0glm_data')
  if(format=='bigmemory')result$descriptor <- bigmemory::describe(object)
  result
}

#' @export
#' @noRd
#' @method dim L0glm_data
dim.L0glm_data <- function(x)c(length(x$rows),length(x$columns))
#' @export
#' @noRd
#' @method dimnames L0glm_data
dimnames.L0glm_data <- function(x) {
  names <- list(x$labels[[1]][x$rows],x$labels[[2]][x$columns])
  if(all(vapply(names,is.null,TRUE)))NULL else names
}
#' @export
#' @noRd
#' @method print L0glm_data
print.L0glm_data <- function(x,...) {
  cat('Bounded',x$format,'design:',nrow(x),'rows and',ncol(x),'columns.\n');invisible(x)
}
#' @export
#' @noRd
#' @method [ L0glm_data
`[.L0glm_data` <- function(x,i,j,...,drop=TRUE) {
  subset <- function(index,value,names) {
    names(index) <- names;answer <- index[value]
    if(anyNA(answer))stop('Unknown or missing data indices.',call.=FALSE)
    unname(answer)
  }
  if(!missing(i))x$rows <- subset(x$rows,i,rownames(x))
  if(!missing(j))x$columns <- subset(x$columns,j,colnames(x))
  if(drop)base::drop(L0glm_data_matrix(x))else x
}
#' @export
#' @noRd
#' @method as.matrix L0glm_data
as.matrix.L0glm_data <- function(x,...)as.matrix(L0glm_data_matrix(x))

L0glm_data_raw <- function(x,rows,columns) {
  if(!length(rows)||!length(columns))return(matrix(numeric(),length(rows),length(columns)))
  if(x$format=='bed')return(l0glm_bed_read_cpp(x$object,x$dimensions[1],x$dimensions[2],
    rows,columns,as.integer(x$control$nthreads)))
  if(x$format=='pgen')return(L0glm_data_pgen(x,rows,columns))
  object <- x$object
  if(x$format=='bigmemory' && bigmemory::is.nil(object@address))
    object <- bigmemory::attach.big.matrix(x$descriptor)
  answer <- object[rows,columns,drop=FALSE]
  if(inherits(answer,'DelayedMatrix'))answer <- if(x$sparse)
    methods::as(answer,'CsparseMatrix')else as.matrix(answer)
  if(!(is.numeric(answer)||inherits(answer,'dMatrix')))
    stop('Predictor blocks must be numeric.',call.=FALSE)
  answer
}

L0glm_data_means <- function(x,columns) {
  unknown <- columns[is.na(x$cache$means[columns])]
  if(length(unknown)) {
    if(!is.null(x$cache$reference)) {
      reference <- x$cache$reference
      x$cache$means[unknown] <- L0glm_data_means(reference,
        reference$columns[match(unknown,x$columns)])
      return(x$cache$means[columns])
    }
    if(x$cache$prediction)stop('Prediction needs training imputation means.',call.=FALSE)
    if(as.double(length(x$cache$rows))*length(unknown)>x$control$max.cells)
      stop('Training imputation block exceeds max.cells.',call.=FALSE)
    block <- L0glm_data_raw(x,x$cache$rows,unknown)
    values <- if(inherits(block,'sparseMatrix'))block@x else block
    if(any(is.infinite(values)))stop('Predictors must be finite or missing.',call.=FALSE)
    mean <- Matrix::colSums(block,na.rm=TRUE)/Matrix::colSums(!is.na(block))
    mean[is.nan(mean)] <- 0
    x$cache$means[unknown] <- mean
  }
  x$cache$means[columns]
}

L0glm_data_matrix <- function(x) {
  if(!inherits(x,'L0glm_data'))return(x)
  if(as.double(nrow(x))*ncol(x)>x$control$max.cells)
    stop('Working design exceeds data max.cells; use SIS/ISIS or a smaller active set.',call.=FALSE)
  answer <- L0glm_data_raw(x,x$rows,x$columns)
  values <- if(inherits(answer,'sparseMatrix'))answer@x else answer
  if(any(is.infinite(values)))stop('Predictors must be finite or missing.',call.=FALSE)
  if(anyNA(values)) {
    if(x$missing=='error')stop('Missing predictors require missing="mean".',call.=FALSE)
    where <- which(is.na(answer),arr.ind=TRUE)
    mean <- L0glm_data_means(x,x$columns)
    answer[where] <- mean[where[,2]]
  }
  dimnames(answer) <- dimnames(x)
  answer
}

L0glm_data_training <- function(x,weights=NULL) {
  if(x$format=='bigmemory' && bigmemory::is.nil(x$object@address))
    x$object <- bigmemory::attach.big.matrix(x$descriptor)
  cache <- new.env(parent=emptyenv());cache$means <- rep(NA_real_,x$dimensions[2])
  positive <- if(is.null(weights)||length(weights)==1L)rep(TRUE,nrow(x))else
    if(is.null(dim(weights)))weights>0 else Matrix::rowSums(weights)>0
  if(length(positive)!=nrow(x))stop('Weights must match source rows.',call.=FALSE)
  cache$rows <- x$rows[positive];cache$prediction <- FALSE;x$cache <- cache;x
}

L0glm_data_restore <- function(fit,X) {
  active <- which(Matrix::rowSums(coef(fit)!=0)>0)-as.integer(fit$added_intercept)
  p <- ncol(X);original <- active[active>0 & active<=p]
  if(!is.null(fit$interactions))original <- unique(c(original,
    unlist(fit$interactions$terms[active[active>p]-p])))
  if(X$missing=='mean' && length(original))L0glm_data_means(X,X$columns[original])
  if(is.null(fit$interactions))fit$X <- L0glm_design_view(X,intercept=isTRUE(fit$added_intercept),
    names=if(is.null(colnames(X)))rep('',p)else colnames(X))else fit$X$data <- X
  fit
}

L0glm_data_newdata <- function(newdata,training) {
  if(is.data.frame(newdata)) {
    if(!all(vapply(newdata,is.numeric,TRUE)))stop('Predictors must be numeric.',call.=FALSE)
    newdata <- as.matrix(newdata)
  }
  if(is.null(dim(newdata))&&is.numeric(newdata))newdata <- if(ncol(training)==1L)
    matrix(newdata,ncol=1L)else matrix(newdata,nrow=1L,dimnames=list(NULL,names(newdata)))
  if(!inherits(newdata,'L0glm_data'))newdata <- l0glmData(newdata,missing=training$missing)
  if(ncol(newdata)!=ncol(training))stop('Prediction must supply the original predictor count.',call.=FALSE)
  expected <- colnames(training);provided <- colnames(newdata)
  if(!is.null(expected)&&!is.null(provided)&&any(!is.na(expected)&nzchar(expected))) {
    if(!anyNA(expected)&&all(nzchar(expected))&&!anyDuplicated(expected)) {
      if(anyNA(provided)||anyDuplicated(provided)||!setequal(expected,provided))
        stop('Data column names must match training predictors.',call.=FALSE)
      newdata <- newdata[,match(expected,provided),drop=FALSE]
    } else if(!identical(expected,provided))stop('Ambiguous predictor names require training order.',call.=FALSE)
  }
  cache <- new.env(parent=emptyenv());cache$means <- rep(NA_real_,newdata$dimensions[2])
  cache$means[newdata$columns] <- training$cache$means[training$columns]
  cache$rows <- integer();cache$prediction <- TRUE
  newdata$cache <- cache;newdata$missing <- training$missing;newdata
}

L0glm_data_validation <- function(X,training) {
  if(!inherits(X,'L0glm_data'))X <- l0glmData(X,missing=training$missing)
  X <- L0glm_data_training(X);X$cache$reference <- training
  X$missing <- training$missing;X
}

L0glm_data_restricted_fit <- function(args) {
  X <- args$X;p <- ncol(X)
  family <- if(is.null(args$family))gaussian()else args$family
  if(!inherits(family,'family'))family <- family()
  intercept <- if(is.null(args$intercept))identical(L0glm_ext_kind(family),'multinomial') ||
    family$family=='binomial'else isTRUE(args$intercept)
  start <- args$start
  eligible <- if(length(start)==1L)rep(start!=0,p)else {
    B <- as.matrix(start)
    if(nrow(B)==p+as.integer(intercept))B <- B[seq_len(p)+as.integer(intercept),,drop=FALSE]
    if(nrow(B)!=p)stop('start must match the source dimensions.',call.=FALSE)
    Matrix::rowSums(B!=0)>0
  }
  mandatory <- L0glm_mandatory(X,args$always.include)
  ids <- sort(unique(c(which(eligible),mandatory)))
  if(!length(ids))ids <- 1L
  if(as.double(nrow(X))*length(ids)>X$control$max.cells)
    stop('Eligible starting support exceeds data max.cells; use SIS/ISIS.',call.=FALSE)
  options <- args[setdiff(names(args),c('X','y','family','weights','offset','nonnegative',
    'screen','nsis','screen.maxit','screen.nthreads'))]
  options$intercept <- intercept
  fit <- L0glm_screen_fit(X,args$y,family,args$weights,args$offset,isTRUE(args$nonnegative),
    'SIS',FALSE,1L,1L,options,args$screen.nthreads,candidates=ids)
  fit$screening$method <- 'start';fit
}

L0glm_data_score <- function(X,context,nsis,full,nthreads) {
  if(X$format=='bed')return(sis_bed_score_cpp(X$object,X$dimensions[1],X$dimensions[2],
    X$rows,X$columns,X$cache$rows,context,nsis,full,nthreads,X$missing=='mean'))
  if(X$format=='bigmemory' && X$missing=='error') {
    d <- X$descriptor@description
    if(identical(d$type,'double') && !isTRUE(d$separated)) {
      object <- X$object
      if(bigmemory::is.nil(object@address))object <- bigmemory::attach.big.matrix(X$descriptor)
      return(sis_bigmemory_score_cpp(object@address,X$rows,X$columns,context,nsis,full,nthreads))
    }
  }
  size <- as.integer(max(1,min(X$control$block.columns,
    floor(min(X$control$block.cells,X$control$max.cells)/nrow(X)))))
  top <- integer();scores <- numeric();statistics <- if(full)numeric(ncol(X))else NULL
  signs <- if(full)integer(ncol(X))else NULL
  threads <- 1L
  for(first in seq.int(1L,ncol(X),by=size)) {
    columns <- seq.int(first,min(ncol(X),first+size-1L))
    block <- L0glm_data_matrix(X[,columns,drop=FALSE])
    part <- sis_score_cpp(L0glm_sis_design(block),context,min(nsis,length(columns)),full,nthreads)
    threads <- max(threads,part$nthreads)
    top <- c(top,columns[part$indices]);scores <- c(scores,part$top_statistics)
    order <- order(-scores,top,method='radix')[seq_len(min(nsis,length(top)))]
    top <- top[order];scores <- scores[order]
    if(full) {statistics[columns] <- part$statistics;signs[columns] <- part$sign}
  }
  list(indices=top,top_statistics=scores,statistics=statistics,sign=signs,nthreads=threads)
}

L0glm_data_marginal <- function(X,y,weights,offset,intercept,always.include) {
  n <- nrow(X);m <- NCOL(y);p <- ncol(X)+as.integer(intercept)
  if(as.double(p)*m>X$control$max.cells)stop('Marginal start exceeds max.cells.',call.=FALSE)
  design <- L0glm_design_view(X,intercept=intercept,
    names=if(is.null(colnames(X)))rep('',ncol(X))else colnames(X))
  free <- L0glm_mandatory(design,if(intercept)c(1L,always.include+1L)else always.include)
  U <- design[,free,drop=FALSE]
  w <- L0glm_ext_matrix(weights,n,m,1,'weights')
  z <- as.matrix(y)-L0glm_ext_matrix(offset,n,m,0,'offset')
  information <- cross <- response <- vector('list',m)
  for(j in seq_len(m)) {
    G <- crossprod(U,U*w[,j]);information[[j]] <- G
    response[[j]] <- if(length(free))z[,j]-as.numeric(U%*%L0glm_nuisance_solve(G,crossprod(U,w[,j]*z[,j])))else z[,j]
  }
  beta <- matrix(0,p,m)
  size <- as.integer(max(1,min(X$control$block.columns,
    floor(min(X$control$block.cells,X$control$max.cells)/n))))
  for(first in seq.int(1L,p,by=size)) {
    ids <- seq.int(first,min(p,first+size-1L));A <- design[,ids,drop=FALSE]
    for(j in seq_len(m)) {
      info <- as.numeric(Matrix::colSums(A^2*w[,j]))
      if(length(free)) {
        C <- as.matrix(Matrix::crossprod(U,A*w[,j]))
        info <- pmax(0,info-colSums(C*L0glm_nuisance_solve(information[[j]],C)))
      }
      rhs <- as.numeric(Matrix::crossprod(A,w[,j]*response[[j]]))
      beta[ids,j] <- ifelse(info>0,rhs/pmax(info,.Machine$double.xmin),0)
    }
  }
  beta[free,] <- 0
  if(length(free)) {
    eta <- matrix(0,n,m)
    for(first in seq.int(1L,p,by=size)) {
      ids <- seq.int(first,min(p,first+size-1L))
      eta <- eta+as.matrix(design[,ids,drop=FALSE]%*%beta[ids,,drop=FALSE])
    }
    for(j in seq_len(m))beta[free,j] <- L0glm_nuisance_solve(information[[j]],crossprod(U,w[,j]*(z[,j]-eta[,j])))
  }
  beta
}
