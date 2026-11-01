#' Estimate constant or time-varying point spread functions
#'
#' Estimate a log-concave PSF from a calibration standard or an overlapping
#' sample, and construct its sparse scan-by-scan convolution design. Numeric
#' retention-time names and the bundled data's `alkanes_TIC` attributes are
#' recognised automatically.
#'
#' @param y Numeric signal vector, scan-by-channel matrix, or a list with `Y`
#'   and `time`. Sparse matrices are accepted. At least nine scans are needed.
#' @param family Identity-link family supported by [splogconNMF()], including
#'   Gaussian, Poisson, quasipoisson, negative binomial, [nb1()], ordinary
#'   [Gamma_censored()] and common-scale [Gamma_scale()]. Character names are
#'   also accepted. Gaussian observations may be negative.
#' @param time Increasing sample retention times. Defaults to numeric names
#'   or row names, or scan indices when names are absent. The standard and
#'   sample must use the same time units.
#' @param standard Optional nonoverlapping calibration signal, with its own
#'   retention-time names or `list(Y=..., time=...)`. Defaults to
#'   `attr(y, "alkanes_TIC")`, supplied with [GCMS_data1] and [GCMS_data2].
#' @param method `"auto"` calibrates from an available standard and then refines
#'   on the sample; without a standard it uses blind sample estimation.
#'   `"standard"` uses calibration alone, treating `y` as the standard when
#'   `standard=NULL`. `"blind"` estimates directly from the overlapping sample.
#'   A supplied `control$kernel_est` takes precedence over automatic calibration.
#' @param nwindows Number of PSF anchors across the sample. `1L` gives a constant
#'   PSF. Multiple anchors are smoothly interpolated on the log-profile scale.
#' @param precluster `FALSE` retains the original channels, `1L` sums all
#'   channels, and an integer requests that many outcome groups at most.
#'   `TRUE` uses at most 32 groups. The default uses the total-ion signal.
#'   Poisson distances use channel totals as size factors and native OpenMP
#'   Ward or Bregman clustering; see [glm_cluster()].
#' @param dispersion,theta Shared noise parameters, as in [splogconNMF()].
#'   Omitted parameters are estimated from temporal differences and held fixed
#'   during PSF estimation. For NB2, `theta` is the shape parameter.
#' @param weights,offset Nonnegative scalar, scan vector or matrix matching the
#'   signal used for estimation. In `method="standard"` these refer to the
#'   standard; in `"auto"` and `"blind"` they refer to the sample. An identity-link
#'   offset is an additive known signal. Grouped offsets are summed.
#' @param ncores Number of native threads; defaults to the physical core count.
#' @param control Named list. `L` fixes an odd kernel length; otherwise standard
#'   support is inferred from its fitted peak masses and blind support defaults
#'   to 41 scans. `kernel_est` supplies a vector or one kernel row per window.
#'   `maxiter=8L` and `tol=1e-3` control blind iterations; `maxiter=0L` returns
#'   the supplied kernel design directly. `fit` is a named list of additional
#'   [L0glm_blind_deco()] fitting options, such as `IC` or `relaxed`.
#'   Standard controls include `n_peaks=NULL`, `baseline_df=NULL`,
#'   `baseline_tau=0.05`, `peak_snr=5`, `peak_height=NULL`,
#'   `peak_min_distance=NULL`, `support_fraction=0.005`, `psf_mass=0.999` and
#'   `support_quantile=0.95`. `linkage="ward"` controls channel clustering.
#'   `censor=0.5` is the Gamma detection limit, taken from a supplied scalar
#'   family constructor unless explicitly overridden.
#' @param verbose Print the estimation mode, aggregation and fitted PSF width.
#'
#' @details
#' **Standard calibration.** A lower-quantile spline estimates background.
#' Isolated calibration peaks are retained using signal-to-noise and profile
#' geometry. Each peak is fitted on its original observation scale:
#' \deqn{\mu_{ij}=b_{ij}+o_{ij}+h_i s_j,\qquad h_i>0,\quad
#' \Delta^2\log h_i\leq0.}
#' Native likelihood-based profile updates minimize
#' \eqn{-\sum_{ij}w_{ij}\ell(y_{ij};\mu_{ij})} while imposing hard log concavity.
#' The logarithm parameterizes the positive profile; observations retain the
#' supplied noise distribution. Unit-area kernels are aligned by their centroid
#' and resampled to the sample scan interval. Isolated-peak fitting, automatic
#' support and interpolation share helpers with [estimate_psf_from_std()].
#'
#' **Blind refinement.** [L0glm_blind_deco()] alternates sparse nonnegative
#' spike fitting with PSF estimation. Completed candidates are compared using
#' the requested information criterion and the best checkpoint is retained.
#' The unpenalized background and any known offset enter both updates.
#' For multiple windows, local kernels are estimated independently, and
#' \eqn{\log g_t=(1-a)\log g_l+a\log g_r}, with smoothstep weights
#' \eqn{a=u^2(3-2u)}, joins adjacent anchors. The returned design has nonnegative,
#' unit-sum, log-concave columns, including truncated boundary columns.
#'
#' **Channel aggregation.** Poisson sums retain their model. Shared-dispersion
#' quasipoisson and NB1 sums retain their mean-variance relationship, Gaussian
#' variances add, and common-scale Gamma shapes add. NB2 and ordinary Gamma
#' use variance matching with observed channel fractions
#' \eqn{r_j=\sum_i y_{ij}/\sum_{i,j\in G}y_{ij}}:
#' \deqn{\theta_G=\theta/\sum_{j\in G}r_j^2,\qquad
#' \phi_G=\phi\sum_{j\in G}r_j^2.}
#' These two sum distributions are approximate when channel means have different
#' temporal shapes. Censored Gamma aggregation is also approximate. Different
#' observation weights within a group are averaged and flagged in the returned
#' diagnostics; retain the original channels for exact weighted fitting.
#'
#' Calibration and sample peaks should represent the same chromatographic
#' conditions. Broad overloaded standard peaks can require sample refinement;
#' `method="standard"` keeps the calibration unchanged. Supplied kernels are
#' useful when a particular experimental PSF is already known.
#'
#' @return A list containing the sparse design `X` (`X_spikes`), kernel vector
#'   or matrix `g_hat`, `psf_kernels`, relative-time `psf_offsets`, anchor
#'   retention times `psf_timepoints`, sample `time`, fitted `width_scans`,
#'   `method`, noise `spec`, aggregation diagnostics, optional `calibration`,
#'   blind `windows`, `iterations`, `converged` and elapsed `seconds`.
#'   The result can be passed directly as `X=psf$X` to [splogconNMF()] or
#'   [L0glm.fit()]. No observed matrix is retained.
#' @references de Rooi, J. J. and Eilers, P. H. C. (2011).
#'   Deconvolution of pulse trains with the L0 penalty.
#'   *Analytica Chimica Acta* 705, 218--226. \doi{10.1016/j.aca.2011.05.030}.
#' @seealso [GCMSdeco()], [L0glm_blind_deco()], [estimate_psf_from_std()],
#'   [splogconNMF()], [glm_cluster()]
#' @example examples/estimate_psf.R
#' @md
#' @export
estimate_psf <- function(y, family=quasipoisson("identity"), time=NULL,
                         standard=attr(y,"alkanes_TIC"),
                         method=c("auto","standard","blind"), nwindows=1L,
                         precluster=1L, dispersion=NULL, theta=NULL,
                         weights=NULL, offset=NULL, ncores=NULL,
                         control=list(), verbose=FALSE) {
  begun<-proc.time()[["elapsed"]];method<-match.arg(method)
  ncores<-L0glm_sis_threads(ncores);L0glm_deco_length(nwindows,"nwindows",1L)
  sample<-.psf_input(y,time);Y<-sample$Y;time<-sample$time;n<-nrow(Y)
  if(nwindows>n)stop("nwindows cannot exceed the number of scans.",call.=FALSE)
  defaults<-list(n_peaks=NULL,baseline_df=NULL,baseline_tau=.05,peak_snr=5,
    peak_min_distance=NULL,peak_height=NULL,support_fraction=.005,
    psf_mass=.999,support_quantile=.95,L=NULL,kernel_est=NULL,
    maxiter=8L,tol=1e-3,censor=.5,linkage="ward",fit=list())
  if(!is.list(control)||length(control)&&(is.null(names(control))||anyDuplicated(names(control))||
    any(!nzchar(names(control))))) stop("Use a named control list.",call.=FALSE)
  unknown<-setdiff(names(control),names(defaults))
  if(length(unknown)) stop("Unknown PSF control: ",paste(unknown,collapse=", "),call.=FALSE)
  if(is.list(family)&&!is.null(family$scalar)&&!"censor"%in%names(control))
    defaults$censor<-family$scalar$censor
  control<-utils::modifyList(defaults,control,keep.null=TRUE)
  L0glm_deco_length(control$maxiter,"maxiter",0L)
  if(!is.numeric(control$tol)||length(control$tol)!=1L||!is.finite(control$tol)||control$tol<0)
    stop("tol must be finite and nonnegative.",call.=FALSE)
  if(!is.numeric(control$censor)||length(control$censor)!=1L||!is.finite(control$censor)||control$censor<0)
    stop("censor must be finite and nonnegative.",call.=FALSE)
  if(!is.numeric(control$support_quantile)||length(control$support_quantile)!=1L||
    !is.finite(control$support_quantile)||control$support_quantile<=0||control$support_quantile>1)
    stop("support_quantile must lie in (0,1].",call.=FALSE)
  if(method=="auto"&&!is.null(standard)&&is.null(control$kernel_est)) {
    calibration<-estimate_psf(y,family,time,standard,"standard",nwindows,precluster,
      dispersion,theta,NULL,NULL,ncores,control,verbose)
    control$kernel_est<-calibration$g_hat
    refined<-estimate_psf(y,family,time,NULL,"blind",nwindows,precluster,
      dispersion,theta,weights,offset,ncores,control,verbose)
    refined$method<-"calibrated_blind";refined$calibration<-calibration$calibration
    refined$calibrated_psf<-calibration[c("psf_kernels","psf_offsets","psf_timepoints","width_scans")]
    refined$seconds<-proc.time()[["elapsed"]]-begun
    return(refined)
  }
  if(method=="auto") method<-"blind"
  if(method=="standard"&&is.null(standard)) standard<-y
  input<-if(method=="standard") .psf_input(standard) else sample
  spec<-slcn_family(family,dispersion,theta)
  W<-slcn_weights(weights,nrow(input$Y),ncol(input$Y))
  O<-slcn_weights(offset,nrow(input$Y),ncol(input$Y),0)
  if(any(!is.finite(W)|W<0)||!any(W>0)||any(!is.finite(O)|O<0))
    stop("Use finite nonnegative weights and offsets.",call.=FALSE)
  if(spec$kind!=0L&&any(input$Y<0)) stop("This noise family requires nonnegative observations.",call.=FALSE)
  if(spec$name%in%c("poisson","nb","nb1")&&any(abs(input$Y-round(input$Y))>1e-7))
    stop("Poisson and NB observations must be integer counts.",call.=FALSE)
  if(spec$name%in%c("Gamma","Gamma_scale")&&
    any(ifelse(input$Y==0,control$censor==0,input$Y<control$censor)))
    stop("Gamma observations must be zero (censored) or at least censor.",call.=FALSE)
  noise<-slcn_noise(input$Y,spec,W,ncores,"none",control$censor);spec<-noise$spec
  aggregated<-.psf_group(input$Y,W,O,spec,precluster,ncores,control$linkage)
  work<-aggregated$Y;W<-aggregated$W;O<-aggregated$O;local_spec<-aggregated$spec
  dt<-median(diff(time))
  if(method=="standard") {
    calibration<-.psf_standard(work,input$time,W,O,local_spec,dt,control,ncores)
    centres<-.5+(seq_len(nwindows)-.5)*n/nwindows
    points<-approx(seq_len(n),time,xout=centres,rule=2)$y
    logs<-t(apply(calibration$psf_kernels,1L,.psf_log_tail))
    kernels<-t(vapply(points,function(t) {
      z<-.psf_interp_log_kernel(t,calibration$peak_timepoints,logs,TRUE)
      h<-exp(z-max(z));h/sum(h)
    },numeric(ncol(logs))))
    offsets<-calibration$psf_offsets;windows<-NULL
  } else {
    L<-control$L
    if(is.null(L)) L<-if(is.null(control$kernel_est)) min(41L,2L*floor((n-1L)/2L)+1L) else
      if(is.matrix(control$kernel_est)) ncol(control$kernel_est) else length(control$kernel_est)
    L0glm_deco_length(L,"L",3L)
    if(is.matrix(control$kernel_est)&&nrow(control$kernel_est)!=nwindows)
      stop("kernel_est must have one row per window.",call.=FALSE)
    if(L%%2L==0L) stop("L must be odd.",call.=FALSE)
    centres<-.5+(seq_len(nwindows)-.5)*n/nwindows
    points<-approx(seq_len(n),time,xout=centres,rule=2)$y
    kernels<-matrix(0,nwindows,L);windows<-vector("list",nwindows)
    for(j in seq_len(nwindows)) {
      half<-max(2L*L,n/nwindows)
      rows<-seq.int(max(1L,ceiling(centres[j]-half)),min(n,floor(centres[j]+half)))
      initial<-if(is.matrix(control$kernel_est)) control$kernel_est[j,] else control$kernel_est
      vector_noise<-max(length(local_spec$dispersion),length(local_spec$theta))>1L
      parts<-if(vector_noise) lapply(seq_len(ncol(work)),identity) else list(seq_len(ncol(work)))
      fitted<-lapply(parts,function(channels) {
      channel_spec<-local_spec
      if(length(channel_spec$dispersion)>1L)channel_spec$dispersion<-channel_spec$dispersion[channels]
      if(length(channel_spec$theta)>1L)channel_spec$theta<-channel_spec$theta[channels]
      args<-utils::modifyList(list(y=work[rows,channels,drop=FALSE],family=.psf_model_family(channel_spec,control$censor),
        kernel_est=initial,L=L,baseline=if(length(rows)<300L) "intercept" else "bspline",
        nknots=if(is.null(control$baseline_df)) max(4L,ceiling(length(rows)/1000)+3L) else control$baseline_df,
        maxiter=control$maxiter,tol=control$tol,nwindows=1L,plot=FALSE,ncores=ncores,
        weights=if(length(W)==1L) rep(W,length(rows)) else W[rows,channels,drop=FALSE],
        offset=if(length(O)==1L) O else O[rows,channels,drop=FALSE],
        dispersion=if(channel_spec$name%in%c("gaussian","quasipoisson")) channel_spec$dispersion else NULL,
        parallelize=FALSE),control$fit)
      if(control$maxiter==0L) {
        if(is.null(initial))stop("Supply kernel_est when maxiter=0.",call.=FALSE)
        list(g_hat=L0glm_deco_kernel(initial,normalize=TRUE),iterations=0L,
          best_iteration=0L,converged=TRUE,selection_trace=numeric(),deviance_trace=numeric())
      } else do.call(L0glm_blind_deco,args)
      })
      value<-Reduce(`+`,lapply(fitted,function(z).psf_log_tail(z$g_hat)))/length(fitted)
      kernels[j,]<-exp(value-max(value));kernels[j,]<-kernels[j,]/sum(kernels[j,])
      windows[[j]]<-list(iterations=max(vapply(fitted,`[[`,0L,"iterations")),
        converged=all(vapply(fitted,`[[`,FALSE,"converged")),
        channels=lapply(fitted,function(z)z[c("best_iteration","selection_trace","deviance_trace")]))
    }
    offsets<-(seq_len(L)-(L+1L)%/%2L)*dt;calibration<-NULL
  }
  X<-.psf_build_sparse_matrix(time,points,offsets,kernels,log_floor=0,
    drop_tol=1e-12,normalize=TRUE,smooth=TRUE)
  X<-slcn_project_design_cpp(as(X,"dgCMatrix"),as.integer(ncores))
  width<-vapply(unique(round(seq(1,n,length.out=min(n,21L)))),function(j) {
    h<-X[,j];sum(h>=max(h)/2)
  },0L)
  if(method=="blind"&&control$maxiter==0L)method<-"fixed"
  if(verbose) cat("PSF:",method,",",nwindows,"anchor(s),",ncol(work),
    "calibration channel(s); median FWHM",median(width),"scans\n")
  list(X=X,X_spikes=X,g_hat=if(nwindows==1L)as.numeric(kernels) else kernels,
    psf_kernels=kernels,psf_offsets=offsets,psf_timepoints=points,time=time,
    method=method,aggregation=aggregated[c("group","original_channels","moment_matched","weight_approximation")],
    spec=spec,calibration=calibration,windows=windows,
    iterations=if(is.null(windows)) 0L else vapply(windows,`[[`,0L,"iterations"),
    converged=is.null(windows)||all(vapply(windows,`[[`,FALSE,"converged")),
    width_scans=median(width),seconds=proc.time()[["elapsed"]]-begun)
}

.psf_input<-function(y,time=NULL) {
  if(is.list(y)&&!is.null(y$Y)) {if(is.null(time))time<-y$time;y<-y$Y}
  if(is.null(time)) time<-suppressWarnings(as.numeric(if(is.null(dim(y))) names(y) else rownames(y)))
  Y<-as.matrix(y)
  if(!is.numeric(Y)||nrow(Y)<9L||!ncol(Y)||any(!is.finite(Y))) stop("Use a finite numeric signal with at least nine scans.",call.=FALSE)
  if(!length(time)) time<-seq_len(nrow(Y))
  if(length(time)!=nrow(Y)||any(!is.finite(time))||any(diff(time)<=0)) stop("Scan times must be finite and increasing.",call.=FALSE)
  storage.mode(Y)<-"double";list(Y=Y,time=time)
}

.psf_group<-function(Y,W,O,spec,k,ncores,linkage) {
  n<-nrow(Y);m<-ncol(Y)
  if(isTRUE(k)) k<-32L
  if(is.null(k)||identical(k,FALSE)) k<-m
  L0glm_deco_length(k,"precluster",1L);k<-min(m,k)
  if(k==m) return(list(Y=Y,W=W,O=O,spec=spec,group=seq_len(m),original_channels=m,moment_matched=FALSE,weight_approximation=FALSE))
  group<-if(k==1L) rep(1L,m) else {
    count<-spec$name%in%c("poisson","quasipoisson","nb","nb1")
    cl<-glm_cluster(if(count)pmax(Y,0) else Y,weights=W,
      family=if(count)poisson("identity") else .psf_model_family(spec,0),
      method=if(count)"poissondist" else "deviance",
      size_fac=if(count)colSums(pmax(Y,0)) else NULL,k=k,
      linkage=linkage,nthreads=ncores)
    cl$clusterids
  }
  P<-sparseMatrix(i=seq_len(m),j=group,x=1,dims=c(m,k));sizes<-tabulate(group,nbins=k)
  sums<-as.matrix(Y%*%P);fractions<-colSums(pmax(Y,0))/as.numeric(colSums(pmax(sums,0)))[group]
  fractions[!is.finite(fractions)]<-0
  sumsq<-as.numeric(rowsum(fractions^2,group,reorder=TRUE))
  if(spec$name=="gaussian") spec$dispersion<-spec$dispersion*sizes
  if(spec$name=="nb") spec$theta<-spec$theta/pmax(sumsq,1e-12)
  if(spec$name=="Gamma") spec$dispersion<-spec$dispersion*pmax(sumsq,1e-12)
  list(Y=sums,W=if(length(W)==1L) W else as.matrix(W%*%P)/rep(sizes,each=n),
    O=if(length(O)==1L) matrix(O,n,k)*rep(sizes,each=n) else as.matrix(O%*%P),
    spec=spec,group=group,original_channels=m,moment_matched=spec$name%in%c("nb","Gamma"),
    weight_approximation=length(W)>1L&&any(vapply(split(seq_len(m),group),
      function(j)any(W[,j,drop=FALSE]!=W[,j[1L]]),FALSE)))
}

.psf_model_family<-function(spec,censor) switch(spec$name,
  gaussian=gaussian("identity"),poisson=poisson("identity"),quasipoisson=quasipoisson("identity"),
  nb=mgcv::nb(theta=spec$theta,link="identity"),
  nb1=nb1(spec$dispersion,link="identity"),
  Gamma=Gamma_censored(spec$dispersion,link="identity",censor=censor),
  Gamma_scale=Gamma_scale(spec$dispersion,link="identity",censor=censor))

.psf_standard<-function(Y,time,W,O,spec,dt,control,ncores) {
  n<-nrow(Y);m<-ncol(Y)
  df<-if(is.null(control$baseline_df)) max(4L,min(20L,ceiling(n/1500))) else control$baseline_df
  corrected<-pmax(Y-O,0)
  baseline<-apply(corrected,2L,function(y).psf_baseline_cobs(time,y,df,control$baseline_tau))
  baseline<-matrix(baseline,n,m)
  tic<-rowSums(pmax(corrected-baseline,0))
  noise<-.psf_noise_scale(rowSums(Y),rowSums(baseline),poisson("identity"))
  height<-control$peak_height
  if(is.null(height))height<-max(control$peak_snr*noise,.01*max(tic))
  at<-.psf_detect_major_peaks(time,tic,noise,control$n_peaks,control$peak_min_distance,height,control$peak_snr)
  windows<-.psf_peak_windows(time,tic,at,noise,control$support_fraction,control$psf_mass)
  windows<-windows[vapply(windows,function(w)length(w$x)>=3L&&sum(w$signal)>0,FALSE)]
  if(!length(windows))stop("No informative calibration peaks were detected.",call.=FALSE)
  quality<-vapply(windows,function(w) {
    slack<-min(1,5*noise/max(w$peak_height,noise))
    spolo_lc(w$signal)$gap<=.25+slack &&
      max(w$signal[c(1L,length(w$signal))])<=max(w$signal)*(.1+slack)
  },FALSE)
  if(!any(quality)) {
    strongest<-which.max(vapply(windows,function(w)w$peak_height,0))
    quality[strongest]<-TRUE
  }
  discarded<-vapply(windows[!quality],`[[`,0,"peak_time")
  windows<-windows[quality]
  fitted<-lapply(windows,function(w) {
    rows<-match(w$x,time);bg<-baseline[rows,,drop=FALSE]+if(length(O)==1L)O else O[rows,,drop=FALSE]
    signal<-pmax(Y[rows,,drop=FALSE]-bg,0);s<-colSums(signal)
    initial<-pmax(rowSums(signal),max(rowSums(signal))*1e-10);initial<-spolo_lc(initial)$h
    mean<-bg+initial%o%s
    z<-slcn_log_profile(Y[rows,,drop=FALSE],mean,initial,s,spec$kind,spec$dispersion,
      spec$theta,if(length(W)==1L)W else as.numeric(W[rows,,drop=FALSE]),
      length(rows),0L,20L,control$censor,12L,as.integer(ncores),mode=0L)
    h<-z$profile/sum(z$profile);centre<-sum(w$x*h)
    q<-.psf_weighted_quantile(w$x,h,c(.0005,.9995))
    list(time=w$peak_time,x=w$x-centre,h=h,half=max(abs(q-centre)),
      height=w$peak_height,width=diff(range(w$x[h>=max(h)/2])))
  })
  half<-as.numeric(quantile(vapply(fitted,`[[`,0,"half"),control$support_quantile,names=FALSE))
  L<-if(is.null(control$L)) 2L*ceiling(half/dt)+1L else control$L
  L0glm_deco_length(L,"L",3L);if(L%%2L==0L)stop("L must be odd.",call.=FALSE)
  offsets<-(seq_len(L)-(L+1L)%/%2L)*dt
  kernels<-t(vapply(fitted,function(z) {
    h<-approx(z$x,z$h,xout=offsets,yleft=0,yright=0)$y
    positive<-which(h>0);core<-seq.int(min(positive),max(positive))
    h[core]<-spolo_lc(h[core])$h;h/sum(h)
  },numeric(L)))
  list(psf_kernels=kernels,psf_offsets=offsets,
    peak_timepoints=vapply(fitted,`[[`,0,"time"),peak_indices=vapply(windows,`[[`,0L,"peak_index"),
    peak_widths=vapply(fitted,`[[`,0,"width"),baseline_df=df,excluded_timepoints=discarded)
}
