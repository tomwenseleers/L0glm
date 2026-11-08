# At lambda=0 ordinary/NB response channels separate. A support change in a few
# channels therefore requires no new coefficient solve in the other channels.
L0glm_local_refit_channels <- function(mask,current,d,control,channels) {
  sub <- current;sub$coefficients <- coef(current)[,channels,drop=FALSE]
  for(name in c("dispersion","dispersion.estimated","theta","theta.estimated"))
    if(length(sub[[name]])==d$m)sub[[name]]<-sub[[name]][channels]
  if(length(sub$families))sub$families<-sub$families[channels]
  for(name in c("y","weights","prior.weights","offset","z","lambdas"))
    if(!is.null(dim(sub[[name]])) && ncol(sub[[name]])==d$m)
      sub[[name]]<-sub[[name]][,channels,drop=FALSE]
  dd <- d;dd$m<-length(channels)
  dd$reference.B<-dd$reference.eta<-NULL
  for(name in c("Y","W","O","allowed"))dd[[name]]<-d[[name]][,channels,drop=FALSE]
  dd$original.y<-as.matrix(d$original.y)[,channels,drop=FALSE]
  dd$estimated<-d$estimated[channels]
  if(length(d$common.phi))dd$common.phi<-d$common.phi[channels]
  if(length(d$poisson.constant)) {
    dd$poisson.Y<-d$poisson.Y[,channels,drop=FALSE]
    dd$poisson.constant<-d$poisson.constant[channels]
  }
  ctrl<-control;ctrl$channel.refit<-FALSE
  updated<-L0glm_local_refit(mask[,channels,drop=FALSE],sub,dd,ctrl)
  answer<-current;B<-as.matrix(coef(current));B[,channels]<-as.matrix(coef(updated))
  answer$coefficients<-methods::as(Matrix::Matrix(B,sparse=TRUE),"CsparseMatrix")
  for(name in c("dispersion","dispersion.estimated","theta","theta.estimated"))
    if(length(answer[[name]])==d$m)answer[[name]][channels]<-updated[[name]]
  if(length(answer$families))answer$families[channels]<-updated$families
  for(name in c("weights","z","lambdas")) {
    a<-answer[[name]];b<-updated[[name]]
    if(is.null(a)||is.null(b))next
    if(is.null(dim(a)))a<-matrix(a,ncol=1L)
    if(is.null(dim(b)))b<-matrix(b,ncol=1L)
    if(ncol(a)==1L && d$m>1L)a<-matrix(a,nrow(a),d$m)
    if(ncol(a)==d$m && nrow(a)==nrow(b))a[,channels]<-b
    answer[[name]]<-a
  }
  answer$metrics<-answer$objective<-NULL
  answer$local.refit<-list(channels=channels,rows=updated$local.refit$rows,convergence=updated$convergence)
  answer
}

# Cheap partial-support trials with every other fitted contribution held in an
# additive LINEAR-PREDICTOR offset. This remains valid for nonidentity links,
# joint multinomial likelihoods and Cox partial likelihoods.
L0glm_local_refit_block <- function(mask,current,d,control) {
  B<-as.matrix(coef(current));previous<-B!=0;previous[d$mandatory,]<-TRUE
  changed<-which(rowSums(mask!=previous)>0)
  candidates<-which(rowSums(mask)>0)
  # Full restricted solves are cheaper for modest supports (e.g. 40-60 ISIS
  # predictors); block trials are reserved for substantially larger supports.
  if(!length(changed)||length(candidates)<=4*control$block.max)return(NULL)
  free<-union(d$mandatory,intersect(changed,candidates))
  rest<-setdiff(candidates,free)
  if(length(rest) && length(free)<control$block.max) {
    Z<-d$D[,rest,drop=FALSE];C<-d$D[,changed,drop=FALSE]
    norms<-as.numeric(Matrix::colSums(Z^2))
    if(inherits(d$D,"sparseMatrix")) {
      touch<-as.numeric(Matrix::rowSums(abs(C))>0)
      overlap<-as.numeric(Matrix::crossprod(Z^2,touch))/pmax(norms,1e-20)
    } else {
      cross<-as.matrix(Matrix::crossprod(Z,C))
      normC<-as.numeric(Matrix::colSums(C^2))
      overlap<-rowSums(cross^2/matrix(pmax(normC,1e-20),length(rest),length(changed),byrow=TRUE))/pmax(norms,1e-20)
    }
    order<-order(overlap,decreasing=TRUE)
    free<-union(free,head(rest[order[overlap[order]>0]],max(0,control$block.max-length(free))))
  }
  free<-sort(free)
  if(!length(free)||length(free)>=length(candidates))return(NULL)
  removed<-B;removed[mask]<-0;removed[free,]<-B[free,,drop=FALSE]
  removed.rows<-which(rowSums(removed!=0)>0)
  shift<-L0glm_local_eta(B,d)-as.matrix(d$D[,removed.rows,drop=FALSE]%*%removed[removed.rows,,drop=FALSE])
  if(d$nonnegative && identical(d$family$link,"identity") && all(d$O>=0) &&
      min(shift)>-1e-10*max(1,max(abs(shift))))shift<-pmax(shift,0)
  dd<-d;dd$D<-d$D[,free,drop=FALSE];dd$X<-dd$D;dd$p<-length(free)
  dd$reference.B<-dd$reference.eta<-NULL
  dd$mandatory<-match(d$mandatory,free);dd$allowed<-d$allowed[free,,drop=FALSE];dd$O<-shift
  sub<-current;sub$coefficients<-coef(current)[free,,drop=FALSE];sub$X<-dd$D
  sub$offset<-shift;sub$always.include<-dd$mandatory;sub$added_intercept<-FALSE
  for(name in c("center","scale"))sub$preprocessing[[name]]<-sub$preprocessing[[name]][free]
  sub$preprocessing$p<-length(free)-length(dd$mandatory)
  if(!is.null(sub$lambdas))sub$lambdas<-matrix(as.numeric(sub$lambdas),d$p,d$m)[free,,drop=FALSE]
  good<-seq_len(d$n)
  fixed.noise<-(is.null(d$kind) && (!any(d$estimated)||length(d$common.phi))) ||
    (identical(d$kind,"nb") && !any(current$theta.estimated))
  if(isTRUE(control$prune.rows) && inherits(dd$D,"sparseMatrix") && fixed.noise) {
    good<-which(Matrix::rowSums(abs(dd$D))>0)
    if(length(good)>0L && length(good)<.8*d$n) {
      dd$D<-dd$D[good,,drop=FALSE];dd$X<-dd$D;dd$n<-length(good)
      for(name in c("Y","W","O"))dd[[name]]<-dd[[name]][good,,drop=FALSE]
      dd$original.y<-as.matrix(dd$original.y)[good,,drop=FALSE]
      sub$X<-dd$D;sub$y<-dd$Y;sub$offset<-dd$O
    } else good<-seq_len(d$n)
  }
  ctrl<-control;ctrl$block.max<-Inf;ctrl$channel.refit<-FALSE
  updated<-L0glm_local_refit(mask[free,,drop=FALSE],sub,dd,ctrl)
  result<-current;out<-B;out[!mask]<-0;out[free,]<-as.matrix(coef(updated))
  result$coefficients<-methods::as(Matrix::Matrix(out,sparse=TRUE),"CsparseMatrix")
  for(name in c("dispersion","dispersion.estimated","theta","theta.estimated","families",
      "convergence","baseline"))if(!is.null(updated[[name]]))result[[name]]<-updated[[name]]
  for(name in c("weights","z"))if(!is.null(updated[[name]])) {
    value<-as.matrix(updated[[name]])
    if(length(good)<d$n) {
      original<-as.matrix(current[[name]])
      if(ncol(original)==1L && d$m>1L)original<-matrix(original,nrow(original),d$m)
      original[good,]<-value;value<-original
    }
    result[[name]]<-value
  }
  if(!is.null(result$lambdas) && !is.null(updated$lambdas)) {
    penalties<-matrix(as.numeric(result$lambdas),d$p,d$m)
    penalties[free,]<-as.numeric(updated$lambdas);penalties[!mask]<-Inf
    result$lambdas<-penalties
  }
  result$metrics<-result$objective<-NULL
  result$local.refit<-list(rows=free,observations=good,conditional=TRUE)
  result
}
