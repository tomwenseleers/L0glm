# Screening is shared across adapters; ISIS uses each adapter's own fitted model.
L0glm_benchmark_path <- function(method, data, grid=NULL) {
  setup <- L0glm_benchmark_setup(method, data)
  options <- data$options[[method]]
  mandatory <- if (method %in% c("L0glm", "abess")) options$always.include else NULL
  mandatory <- L0glm_mandatory(data$X, mandatory)
  start <- if (method=="L0glm") options$start else NULL
  mode <- L0glm_screen_method(data$screen, data$X, start, data$intercept, mandatory)
  if (mode=="none") return(L0glm_benchmark_path_raw(method, data, grid))
  p <- ncol(data$X); n <- nrow(data$X); added <- as.integer(data$intercept)
  allowed <- rep(TRUE, p)
  if (!is.null(start)) {
    if (any(!is.finite(start))) stop("start must contain finite coefficients.",call.=FALSE)
    if (length(start)==1L) allowed[] <- start!=0 else {
      B <- as.matrix(start)
      if (nrow(B)==p+added) B <- B[seq_len(p)+added, , drop=FALSE]
      if (nrow(B)!=p) stop("start must have the original coefficient dimensions.", call.=FALSE)
      allowed <- rowSums(B!=0)>0
    }
  }
  allowed[mandatory] <- TRUE
  eligible <- setdiff(which(allowed), mandatory)
  size <- data$nsis
  if (is.null(size)) size <- if (mode=="SIS") 2*floor(n/log(max(n,3L))) else
    max(1, round(n/(log(log(max(n,3L)))*log(max(data$p,2L)))))
  size <- min(size, length(eligible))
  niter <- if (mode=="SIS") 1L else data$screen.maxit
  criterion <- if (data$tune=="cv") "gic" else data$tune
  threads <- 1L; selected <- mandatory; previous <- NULL; best <- NULL
  offset <- data$offset; history <- list(); converged <- FALSE
  nuisance <- if (!length(mandatory)) NULL else
    cbind(if (data$kind!="cox") matrix(1,n,1L), data$X[,mandatory,drop=FALSE])
  screening_X <- if (length(eligible)==p) data$X else data$X[,eligible,drop=FALSE]
  for (iteration in seq_len(niter)) {
    proposal <- if (size) sisGLM(screening_X, data$y, family=data$family,
      weights=data$weights, offset=offset, nuisance=nuisance, nsis=size,
      return_stats=FALSE, nonnegative=data$nonnegative,
      nthreads=data$screen.nthreads) else NULL
    if (!is.null(proposal)) threads <- max(threads, proposal$nthreads)
    candidates <- sort(unique(c(selected, mandatory, eligible[proposal$indices])))
    if (identical(candidates, previous)) {converged <- TRUE; break}
    # A zero start can request a null model without calling a path backend.
    if (!length(candidates)) {
      B <- Matrix::Matrix(0, p+added, data$m, sparse=TRUE)
      if (added) {
        null <- l0glmFit(matrix(1,n,1L), data$y, family=data$family,
          weights=data$weights, offset=data$offset, lambda=0, nonnegative=data$nonnegative,
          intercept=FALSE, screen="none",local.search=FALSE)
        B[1L,] <- coef(null)
      }
      path <- list(grid=if (is.null(grid)) 0 else grid,
        coefficients=rep(list(B), if (is.null(grid)) 1L else length(grid)),
        parameter=if (method=="abess") "support.size" else "lambda",
        scope=setup$scope, approximation=setup$approximation, backend=setup$backend)
    } else {
      d <- L0glm_benchmark_subset(data, candidates)
      # Preset calibration retains the original search dimension.
      if (method=="L0glm" && is.null(grid)) {
        preset <- presetLambda(criterion, data$p, n)
        grid <- exp(seq(log(max(2,preset)*100),log(2),length.out=data$nlambda))
      }
      path <- L0glm_benchmark_path_raw(method, d, grid)
      path$metrics <- lapply(path$coefficients,function(B) if (is.null(B)) NULL else
        l0glmMetrics(L0glm_benchmark_fit(B,d),mu_floor=data$mu_floor))
      path$coefficients <- lapply(path$coefficients, L0glm_benchmark_expand,
        candidates=candidates, p=p, added=added, m=data$m)
      if (!is.null(path$adaptive)) {
        path$adaptive$coefficients <- L0glm_benchmark_expand(path$adaptive$coefficients,
          candidates,p,added,data$m)
        factors <- path$adaptive$penalty.factor
        full <- if (is.null(dim(factors))) rep(Inf,p) else matrix(Inf,p,ncol(factors))
        if (is.null(dim(factors))) full[candidates] <- factors else full[candidates,] <- factors
        path$adaptive$penalty.factor <- full
      }
    }
    if (is.null(path$metrics)) path$metrics <- lapply(path$coefficients,function(B)
      l0glmMetrics(L0glm_benchmark_fit(B,data),mu_floor=data$mu_floor))
    scores <- vapply(path$metrics,function(z) if (is.null(z)) Inf else z[[criterion]],numeric(1))
    available <- !is.na(scores) & scores<Inf
    if (!any(available)) stop("No screened candidate has an available information criterion.",call.=FALSE)
    index <- which.min(replace(scores,!available,Inf)); score <- scores[index]
    B <- path$coefficients[[index]]
    selected <- which(Matrix::rowSums(B[seq_len(p)+added,,drop=FALSE]!=0)>0)
    history[[iteration]] <- data.frame(iteration, candidates=length(candidates),
      selected=length(setdiff(selected,mandatory)), score=score)
    if (is.null(best) || score<best$score)
      best <- list(path=path, score=score, iteration=iteration, candidates=candidates)
    active <- which(Matrix::rowSums(B!=0)>0)
    offset <- as.matrix(data$design[,active,drop=FALSE]%*%B[active,,drop=FALSE])+data$offset
    previous <- candidates
    if (length(candidates)==sum(allowed)) {converged <- TRUE; break}
  }
  path <- best$path
  path$screening <- list(method=mode, screen_size=size, nthreads=threads,
    automatic=tolower(data$screen)=="auto",
    iterations=length(history), best_iteration=best$iteration, converged=converged,
    candidates=best$candidates, criterion=criterion, history=do.call(rbind,history))
  path
}

L0glm_benchmark_expand <- function(B, candidates, p, added, m) {
  if (is.null(B)) return(NULL)
  B <- methods::as(methods::as(Matrix::Matrix(B,sparse=TRUE),"generalMatrix"),"CsparseMatrix")
  entries <- summary(B)
  rows <- c(if (added) 1L, candidates+added)
  Matrix::sparseMatrix(i=rows[entries$i],j=entries$j,x=entries$x,dims=c(p+added,m))
}

L0glm_benchmark_subset <- function(data, candidates) {
  d <- data; p <- ncol(data$X); added <- as.integer(data$intercept)
  d$X <- data$X[,candidates,drop=FALSE]
  d$design <- if (added) cbind(Intercept=1,d$X) else d$X
  d$support.size <- data$support.size[data$support.size<=length(candidates)]
  d$options <- lapply(data$options,function(options) {
    for (name in c("start","prior.mean")) {
      B <- options[[name]]
      if (is.null(B) || length(B)==1L) next
      B <- as.matrix(B)
      if (added && nrow(B)==p) B <- rbind(0,B)
      rows <- if (added) c(1L,candidates+added) else candidates
      options[[name]] <- B[rows,,drop=FALSE]
    }
    for (name in c("always.include","init.active.set")) if (!is.null(options[[name]])) {
      index <- match(options[[name]],candidates)
      options[[name]] <- index[!is.na(index)]
    }
    if (!is.null(options$group.index)) options$group.index <- options$group.index[candidates]
    options
  })
  d
}
