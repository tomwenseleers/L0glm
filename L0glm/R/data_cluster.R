L0glm_data_cluster <- function(Y,weights,family,method,size_fac,rowbins,nthreads,
  linkage,criterion,k,k.max,dispersion,distance.only=FALSE) {
  if(!inherits(family,'family'))family <- family()
  n <- nrow(Y);m <- ncol(Y);info <- L0glm_cluster_family(family)
  if(is.null(method))method <- if(family$family=='poisson')'poissondist'else 'deviance'
  method <- match.arg(method,c('poissondist','deviance'))
  if(method=='poissondist' && info$id!=1L)stop('poissondist needs Poisson.',call.=FALSE)
  if(is.null(linkage))linkage <- 'ward'
  if(linkage!='ward')stop('File-backed clustering supports Ward; materialize a bounded design for Bregman merges.',call.=FALSE)
  if(is.null(criterion))criterion <- if(info$likelihood)'gic'else 'dynamic'
  criterion <- match.arg(tolower(criterion),c('gic','aic','bic','ebic','mbic','dynamic'))
  if(!info$likelihood && criterion!='dynamic' && is.null(k))stop('Use a fixed or dynamic cut for this family.',call.=FALSE)
  for(value in list(rowbins,k.max,if(is.null(k))1 else k))if(length(value)!=1L ||
    !is.numeric(value)||!is.finite(value)||value<1||value!=floor(value))stop('Invalid clustering count.',call.=FALSE)
  if(!is.null(k)&&k>m)stop('k exceeds the channel count.',call.=FALSE)
  if(!is.numeric(dispersion)||!length(dispersion)%in%c(1L,m)||any(!is.finite(dispersion)|dispersion<=0))
    stop('Invalid dispersion.',call.=FALSE)
  if(info$likelihood&&info$id%in%c(1L,3L,6L)&&family$family!='quasipoisson'&&any(dispersion!=1))
    stop('This family has fixed dispersion one.',call.=FALSE)
  if(!is.null(info$scalar)&&any(dispersion!=1))stop('Specify scalar noise in the family.',call.=FALSE)
  w <- L0glm_cluster_weights(weights,n,m);threads <- L0glm_sis_threads(nthreads)
  if(rowbins>1L && (info$id!=1L||any(w!=1)||length(dispersion)!=1L))
    stop('Row binning requires Poisson, unit weights and common dispersion.',call.=FALSE)
  size <- floor(min(Y$control$block.cells,Y$control$max.cells)/m)
  size <- as.integer(max(rowbins,floor(size/rowbins)*rowbins))
  if(as.double(size)*m>Y$control$max.cells)stop('A clustering block exceeds max.cells.',call.=FALSE)
  starts <- seq.int(1L,n,by=size)
  read <- function(first,work=TRUE) {
    rows <- seq.int(first,min(n,first+size-1L));A <- as.matrix(Y[rows,,drop=FALSE])
    W <- if(length(w)==1L)w else if(length(w)==n)w[rows]else matrix(w,n,m)[rows,,drop=FALSE]
    if(work && rowbins>1L) {
      A <- rowsum(A,rep(seq_len(ceiling(length(rows)/rowbins)),each=rowbins,length.out=length(rows)),reorder=FALSE)
      W <- 1
    }
    precision <- if(length(dispersion)==1L)W/dispersion else
      matrix(W,nrow(A),m)/rep(dispersion,each=nrow(A))
    list(Y=A,W=W,precision=as.numeric(precision),rows=rows,
      work.rows=seq.int((first-1L)%/%rowbins+1L,length.out=nrow(A)))
  }
  sizes <- rep(1,m)
  if(method=='poissondist') {
    if(is.null(size_fac)) {
      sizes[] <- 0
      for(first in starts)sizes <- sizes+colSums(read(first,FALSE)$Y)
    } else sizes <- size_fac
    if(length(sizes)!=m||any(!is.finite(sizes)|sizes<0))stop('Invalid size factors.',call.=FALSE)
  }
  if(as.double(m)*(m-1)/2>Y$control$max.cells)stop('Pairwise distances exceed max.cells.',call.=FALSE)
  values <- numeric(as.double(m)*(m-1)/2)
  for(first in starts) {
    block <- read(first)
    values <- values+as.numeric(glmDist(block$Y,block$precision,family,method,sizes,threads))
  }
  D <- structure(values,Size=m,Labels=colnames(Y),Diag=FALSE,Upper=FALSE,method=method,class='dist')
  if(distance.only)return(D)
  hc <- if(m>1L)fastcluster::hclust(D,method='ward.D')else NULL
  profile <- function(groups,centres=FALSE) {
    groups <- matrix(as.integer(groups),m)
    loss <- df <- numeric(ncol(groups));C <- NULL
    if(centres) {
      if(as.double(ceiling(n/rowbins))*max(groups)>Y$control$max.cells)
        stop('Cluster centres exceed max.cells.',call.=FALSE)
      C <- matrix(0,ceiling(n/rowbins),max(groups))
    }
    for(first in starts) {
      b <- read(first)
      result <- if(is.null(info$scalar))L0glm_cluster_profile_cpp(b$Y,b$precision,sizes,groups,info$id,info$theta,threads)else {
        s <- info$scalar;L0glm_scalar_profile_cpp(b$Y,b$precision,groups,s$kind,s$parameter,s$censor,threads)
      }
      loss <- loss+result$deviance;df <- df+result$df
      if(centres)C[b$work.rows,] <- result$centres
    }
    list(deviance=loss,df=df,centres=C)
  }
  history <- NULL
  if(criterion=='dynamic'&&is.null(k)&&m>1L) {
    if(as.double(m)*m>Y$control$max.cells)stop('Dynamic-cut distances exceed max.cells.',call.=FALSE)
    group <- dynamicTreeCut::cutreeDynamic(hc,distM=as.matrix(D),deepSplit=2,minClusterSize=1,verbose=0)
    absent <- which(group==0L);group[absent] <- max(group)+seq_along(absent)
    group <- match(group,sort(unique(group)))
  } else {
    ks <- if(m==1L)1L else if(!is.null(k))as.integer(k)else seq_len(min(k.max,m))
    cuts <- if(m==1L)matrix(1L,1,1)else matrix(stats::cutree(hc,k=ks),m)
    best <- 1L
    if(info$id>0L||!is.null(info$scalar)) {
      result <- profile(cuts)
      allocation <- vapply(seq_along(ks),function(i){z<-tabulate(cuts[,i]);-2*sum(z*log(z/m))},numeric(1))
      df <- result$df+ks-1L
      score <- if(info$likelihood)get.IC(-(result$deviance+allocation)/2,max(3,n),max(2,n),df)[[
        if(criterion=='dynamic')'gic'else criterion]]else rep(NA_real_,length(ks))
      history <- data.frame(k=ks,deviance=result$deviance,allocation_loss=allocation,df=df,IC=score)
      if(info$likelihood)best <- which.min(score)
    }
    group <- cuts[,best]
  }
  centres <- if(info$id>0L||!is.null(info$scalar))profile(group,TRUE)$centres else NULL
  groups <- max(group)
  if(as.double(n)*groups>Y$control$max.cells)stop('Returned aggregates exceed max.cells.',call.=FALSE)
  aggregate <- matrix(0,n,groups);precision <- if(is.null(weights)&&info$id!=2L)NULL else aggregate
  note <- NULL
  for(first in starts) {
    b <- read(first,FALSE)
    value <- tryCatch(glmAggregate(b$Y,if(info$id==2L)b$precision else b$W,family,group),error=function(e) {
      if(grepl('Count aggregation requires shared',conditionMessage(e),fixed=TRUE))return(NULL)
      stop(e)
    })
    if(is.null(value)) {note <- 'Count weights differ within groups; aggregate original weighted responses separately.';break}
    aggregate[b$rows,] <- value$Y_agg
    if(!is.null(precision))precision[b$rows,] <- value$weights_agg
  }
  dimnames(aggregate) <- list(rownames(Y),paste0('cluster_',seq_len(groups)))
  if(!is.null(precision))dimnames(precision) <- dimnames(aggregate)
  list(D=if(m>1L)D else NULL,hc=hc,clusterids=as.integer(group),history=history,centres=centres,
    size_fac=if(method=='poissondist')sizes else NULL,linkage='ward',
    criterion=if(is.null(k))criterion else 'fixed',rowbins=as.integer(rowbins),
    Y_agg=if(is.null(note))aggregate else NULL,weights_agg=if(is.null(note))precision else NULL,
    aggregation_note=note)
}
