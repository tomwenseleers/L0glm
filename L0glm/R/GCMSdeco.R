GCMS_range <- function(value, name) {
  if (!is.null(value) && (!is.numeric(value) || length(value)!=2L || anyNA(value) || value[1L]>=value[2L]))
    stop(name, " must contain increasing lower and upper retention times in minutes.",call.=FALSE)
  value
}

GCMS_read_cdf <- function(cdf_file, time_range=NULL, mass=NULL) {
  if (!requireNamespace("ncdf4",quietly=TRUE)) stop("Install ncdf4 to read ANDI-MS files.",call.=FALSE)
  if (!is.character(cdf_file)||length(cdf_file)!=1L||is.na(cdf_file)||!nzchar(cdf_file))
    stop("cdf_file must be a local file path or HTTP(S) URL.",call.=FALSE)
  local <- cdf_file
  if (grepl("^https?://",cdf_file,ignore.case=TRUE)) {
    local <- tempfile(fileext=".cdf");on.exit(unlink(local),add=TRUE)
    utils::download.file(cdf_file,local,mode="wb",quiet=TRUE)
  }
  nc <- ncdf4::nc_open(local)
  on.exit(ncdf4::nc_close(nc),add=TRUE,after=FALSE)
  required <- c("point_count","scan_index","scan_acquisition_time","mass_values","intensity_values")
  if (!all(required %in% c(names(nc$var),names(nc$dim))))
    stop("The CDF needs ANDI-MS scan, mass and intensity variables.",call.=FALSE)
  get <- function(name,...) ncdf4::ncvar_get(nc,name,...)
  counts <- as.numeric(get("point_count"));starts <- as.numeric(get("scan_index"))
  time <- as.numeric(get("scan_acquisition_time"));n <- length(time)
  points <- prod(nc$var[["mass_values"]]$varsize)
  if (!n || length(counts)!=n || length(starts)!=n ||
      any(!is.finite(counts)|counts<0|counts!=floor(counts)) ||
      any(!is.finite(starts)|starts<0|starts!=floor(starts)) ||
      any(starts+counts>points) || any(diff(starts)<head(counts,-1L)) ||
      prod(nc$var[["intensity_values"]]$varsize)!=points)
    stop("CDF scan indices and point counts are inconsistent.",call.=FALSE)
  if (any(!is.finite(time)|time<0)||any(diff(time)<=0)) stop("Scan times must be finite and increasing.",call.=FALSE)
  att <- ncdf4::ncatt_get(nc,"scan_acquisition_time","units")
  units <- if(isTRUE(att$hasatt)) tolower(trimws(as.character(att$value))) else ""
  factor <- if (!nzchar(units)||grepl("^s(ec(ond)?s?)?$",units)) 1/60 else
    if(grepl("^min(ute)?s?$",units)) 1 else if(grepl("^ms$|^milliseconds?$",units)) 1/60000 else
      stop("Unrecognized scan time units: ",units,call.=FALSE)
  minutes <- time*factor;time_range <- GCMS_range(time_range,"time_range")
  rows <- if(is.null(time_range)) seq_len(n) else which(minutes>=time_range[1L]&minutes<=time_range[2L])
  if(!length(rows)) stop("No scans in the requested retention-time range.",call.=FALSE)
  first <- starts[rows[1L]];span <- max(starts[rows]+counts[rows])-first
  if(span<=0||sum(counts[rows])==0) {
    if(is.null(mass)) stop("No mass/intensity points in the requested scans.",call.=FALSE)
    return(list(Y=matrix(0,length(rows),length(mass)),time=minutes[rows],mass=mass,
      source=cdf_file,scans=rows,scan_time_units=units))
  }
  take <- sequence(as.integer(counts[rows]))+rep(starts[rows]-first,counts[rows])
  mz <- as.numeric(get("mass_values",start=first+1,count=span))[take]
  intensity <- as.numeric(get("intensity_values",start=first+1,count=span))[take]
  if(any(!is.finite(mz)|mz<0)||any(!is.finite(intensity)|intensity<0))
    stop("CDF masses and intensities must be finite and nonnegative.",call.=FALSE)
  nominal <- floor(mz+.5)
  if(is.null(mass)) mass <- seq.int(min(nominal),max(nominal))
  channel <- match(nominal,mass);keep <- !is.na(channel)
  Y <- as.matrix(Matrix::sparseMatrix(i=rep.int(seq_along(rows),counts[rows])[keep],
    j=channel[keep],x=intensity[keep],dims=c(length(rows),length(mass))))
  dimnames(Y) <- list(format(minutes[rows],digits=12,trim=TRUE),as.character(mass))
  list(Y=Y,time=minutes[rows],mass=mass,source=cdf_file,scans=rows,
    scan_time_units=if(nzchar(units)) units else "ANDI-MS seconds (units omitted)")
}

GCMS_compounds <- function(model,time,sample_name="L0glm") {
  apex <- vapply(seq_len(ncol(model$H)),function(k) which.max(model$H[,k]),0L)
  order <- order(time[apex]);area <- rowSums(model$S);delta <- diff(time)
  quadrature <- if(length(time)>1L) c(delta[1L]/2,
    if(length(time)>2L) (head(delta,-1L)+tail(delta,-1L))/2,tail(delta,1L)/2) else 0
  integrated <- as.numeric(Matrix::crossprod(model$H,quadrature))*area
  names <- sprintf("%s_%04d_RT%.5f",sample_name,seq_along(order),time[apex[order]])
  data.frame(sample=rep(sample_name,length(order)),compound=order,name=names,
    retention_time=time[apex[order]],apex_scan=apex[order],area=area[order],
    area_min=integrated[order],n_ions=rowSums(model$S>0)[order],row.names=NULL)
}

GCMS_write_msp <- function(model,time,mass,msp_file,compounds=NULL) {
  if(nrow(model$H)!=length(time)||ncol(model$S)!=length(mass))
    stop("Profile, spectrum and scan/mass dimensions differ.",call.=FALSE)
  if(!is.character(msp_file)||length(msp_file)!=1L||is.na(msp_file)||!nzchar(msp_file)||
      !grepl("\\.msp$",msp_file,ignore.case=TRUE)) stop("msp_file must end in .msp.",call.=FALSE)
  if(is.null(compounds)) compounds <- GCMS_compounds(model,time)
  order <- compounds$compound;names <- compounds$name
  lines <- unlist(lapply(seq_along(order),function(index) {
    k<-order[index];s<-model$S[k,];keep<-which(s>0)
    if(!length(keep)) stop("Cannot export an empty compound spectrum.",call.=FALSE)
    abundance<-999*s[keep]/max(s[keep])
    c(paste0("Name: ",names[index]),
      sprintf("Comments: Sample %s; retention time %.7g min; integrated ion intensity %.9g; splogconNMF",
        compounds$sample[index],compounds$retention_time[index],sum(s)),
      sprintf("Num Peaks: %d",length(keep)),sprintf("%d %.8g;",mass[keep],abundance),"")
  }),use.names=FALSE)
  dir.create(dirname(msp_file),recursive=TRUE,showWarnings=FALSE)
  con<-file(msp_file,"wb");on.exit(close(con))
  writeLines(as.character(lines),con,sep="\r\n",useBytes=TRUE)
  compounds
}

GCMS_calibration_start <- function(calibration,time,nwindows,L=NULL) {
  dt<-stats::median(diff(time))
  if(is.null(L)) L<-2L*ceiling(max(abs(calibration$psf_offsets))/dt)+1L
  if(L%%2L==0L)L<-L+1L
  centres<-.5+(seq_len(nwindows)-.5)*length(time)/nwindows
  times<-stats::approx(seq_along(time),time,xout=centres,rule=2)$y
  log_kernels<-t(apply(calibration$psf_kernels,1L,.psf_log_tail))
  offsets<-(seq_len(L)-(L+1L)%/%2L)*dt
  kernels<-t(vapply(times,function(t) {
    z<-.psf_interp_log_kernel(t,calibration$peak_timepoints,log_kernels,smooth=TRUE)
    h<-exp(stats::approx(calibration$psf_offsets,z-max(z),xout=offsets,
      yleft=-Inf,yright=-Inf)$y)
    positive<-which(h>0)
    if(!length(positive)) stop("The calibrated PSF has no positive sample-grid support.",call.=FALSE)
    core<-seq.int(min(positive),max(positive));h[core]<-spolo_lc(h[core])$h
    h/sum(h)
  },numeric(L)))
  if(nwindows==1L) as.numeric(kernels) else kernels
}

#' Deconvolve a GC/MS chromatogram and export its compound spectra
#'
#' Read an ANDI-MS netCDF chromatogram, estimate a smoothly time-varying point
#' spread function, recover log-concave compounds with [splogconNMF()], and
#' write a NIST-compatible MSP spectral library. All original mass channels
#' are retained in the final spectra.
#'
#' @param cdf_file Local ANDI-MS \code{.cdf} path or HTTP(S) URL. The optional
#'   \pkg{ncdf4} package is required. Masses are binned to the nearest integer
#'   m/z; traces between the minimum and maximum recorded mass are retained.
#' @param msp_file Output \code{.msp} path. Parent directories are created if
#'   necessary; an existing output file is replaced.
#' @param xlsx_file Optional Excel output path. The workbook contains compound
#'   areas, unit-L1 spectra and run metadata, including the sample name. Requires
#'   the optional \pkg{openxlsx} package. \code{NULL} omits Excel export.
#' @param sample_name Sample identifier, defaulting to the CDF filename without
#'   its extension. Included in the compound table, MSP entries and workbook.
#' @param normalize Spectral coefficient normalization: \code{"none"} returns
#'   integrated channel intensities; \code{"l1"} gives a unit-sum spectrum per
#'   compound, with its magnitude carried by the compound area.
#' @param family Identity-link noise model accepted by [splogconNMF()].
#'   The default \code{quasipoisson("identity")} estimates one pooled dispersion
#'   across channels, suitable for overdispersed instrumental intensities.
#' @param nwindows Number of overlapping blind-PSF windows, default \code{10L}.
#'   \code{1L} uses a shared PSF. The returned sparse design interpolates the
#'   PSFs smoothly across the chromatogram.
#' @param alkanes_file Optional ANDI-MS alkane-standard chromatogram. Isolated
#'   peaks calibrate the initial PSF through [estimate_psf_from_std()], with
#'   automatic support and resampling to the sample's scan interval.
#' @param time_range Optional increasing pair of sample retention-time bounds
#'   in minutes. Infinite bounds are accepted. This can exclude the solvent
#'   interval or select a chromatographic window.
#' @param calibration_range Optional retention-time bounds in minutes for the
#'   standard. The whole standard is used when omitted.
#' @param baseline_df Dimension of the nonnegative cubic B-spline background.
#'   \code{NULL} chooses a low-dimensional basis from the scan count and typical
#'   PSF width; an integer of at least four fixes its dimension.
#' @param psf_control Named list passed to [L0glm_blind_deco()], for example
#'   \code{list(L=61L, maxiter=8L, tol=1e-3)}. \code{kernel_est} supplies a known
#'   initial PSF. Data and \code{nwindows} are set by this function. The scalar
#'   Poisson working weight is matched to the final noise model by default;
#'   \code{weights} can override it. Multiple windows require a scalar weight.
#' @inheritParams splogconNMF
#' @param object,x A fitted \code{GCMSdeco} object.
#'
#' @details
#' **File formats.** ANDI-MS variables \code{scan_index}, \code{point_count},
#' \code{scan_acquisition_time}, \code{mass_values} and \code{intensity_values}
#' define the scan spectra. Scan indices are zero-based. Padding between scans
#' is accepted. Intensities must be finite and nonnegative. Times are converted
#' to minutes; absent units use the ANDI-MS seconds convention. The selected
#' retention-time interval is read directly from the point arrays.
#'
#' MSP entries are ordered by profile apex time, with a unique compound name,
#' retention time and integrated ion intensity in the comments, followed by
#' mass/intensity pairs. Exported spectra have base peak 999 and retain every
#' positive fitted ion. \code{coef()} returns unnormalized fitted spectra.
#' With \code{normalize="l1"}, \code{coef()} returns spectra summing to one.
#' Their compound profiles are \eqn{P_{ir}=H_{ir}\sum_j S_{rj}}, preserving
#' \eqn{HS=PT}. The workbook reports \code{area}, equal to \eqn{\sum_j S_{rj}}, in
#' intensity-scans, and \code{area_min}, the trapezoidal integral of
#' \eqn{P_{\cdot r}} over retention time, in intensity-minutes.
#'
#' The returned object stores the factors and source CDF path. Observed data
#' are loaded from the CDF only when a plot needs them. Fitted signals are
#' reconstructed from the factors; the observed matrix is released after fitting.
#'
#' **PSF and background.** [L0glm_blind_deco()] fits the total ion chromatogram
#' in overlapping windows, optionally starting from the alkane standard.
#' Interpolated columns are projected to log concavity on their local support
#' and normalized to unit sum. Their typical full width at half maximum,
#' \eqn{w}, sets the background dimension to
#' \deqn{q=\max\{4,\min[40,3+\lceil n/(200w)\rceil]\}.}
#' Background coefficients are nonnegative and unpenalized. This keeps the
#' baseline much smoother than an individual elution peak.
#'
#' **Coherent noise scaling.** Noise parameters are estimated once from the
#' original channels and shared by the PSF pilot and final factorization.
#' If \eqn{\widehat\mu_{i+}} and \eqn{\widehat V_{i+}} are pilot total means
#' and variances, the Poisson working weight for PSF estimation is
#' \deqn{a=\frac{\sum_i\widehat\mu_{i+}}
#' {\kappa\sum_i\widehat V_{i+}}.}
#' For quasipoisson and NB1 this becomes \eqn{1/(\kappa\phi)}, and for
#' common-scale Gamma it is \eqn{1/(\kappa b)}. The lag-one adjustment is
#' \eqn{\kappa=1+2\widehat\rho_1}. Gaussian, NB2 and standard Gamma use their
#' summed pilot variances in the same formula. This gives the PSF pilot an
#' appropriate average precision; final factorization uses the specified
#' family's full original-channel likelihood.
#'
#' \code{precluster=32L} performs selection and fusion in at most 32 size-adjusted
#' groups, followed by refinement on every original channel. Rank selection,
#' hard log concavity, censored Gamma intensities and compound-count information
#' criteria follow [splogconNMF()].
#'
#' @return An object of class \code{GCMSdeco} containing \code{model} (the full
#'   factorization), a \code{compounds} table with names, retention times,
#'   original-file apex scan indices, integrated intensities and ion counts,
#'   the blind fit and sparse design in \code{psf}, optional standard calibration
#'   in \code{psf_standard}, \code{background_df}, \code{psf_width_scans},
#'   the actual PSF working weight, retention-time and mass axes, original scan
#'   indices, sample identity, output file paths and stage timings. Coefficient, fitted-value,
#'   log-likelihood and information-criterion methods delegate to \code{model}.
#' @references NIST MS Search user manual, section on NIST text format:
#'   \url{https://chemdata.nist.gov/mass-spc/ms-search/docs/Ver20Man.pdf}.
#' @seealso [plot.GCMSdeco()], [splogconNMF()], [L0glm_blind_deco()], [estimate_psf_from_std()]
#' @examples
#' \dontrun{
#' fit <- GCMSdeco("beer.cdf", "beer.msp", alkanes_file="alkanes.cdf",
#'                 time_range=c(6, Inf), xlsx_file="beer.xlsx")
#' fit
#' fit$compounds
#' coef(fit)                 # unnormalized spectra, all channels
#' GIC(fit)
#' plot(fit, type="lines", file="beer_traces.pdf")
#' plot(fit, type="heatmap", file="beer_heatmaps.pdf")
#' matplot(fit$time, as.matrix(fit$model$H), type="l", lty=1,
#'         xlab="Retention time (min)", ylab="Elution profile")
#' }
#' @export
GCMSdeco <- function(cdf_file,msp_file,family=quasipoisson("identity"),lambda="gic",nwindows=10L,
                     alkanes_file=NULL,time_range=NULL,calibration_range=NULL,
                     baseline_df=NULL,precluster=32L,dispersion=NULL,theta=NULL,
                     correlation=c("auto","none","lag1"),
                     ncores=NULL,control=list(),psf_control=list(),verbose=FALSE,
                     xlsx_file=NULL,sample_name=NULL) {
  begun<-proc.time()[["elapsed"]];seconds<-numeric()
  ncores<-L0glm_sis_threads(ncores);L0glm_deco_length(nwindows,"nwindows",1L)
  correlation<-slcn_correlation(correlation)
  if(is.null(sample_name)) sample_name<-tools::file_path_sans_ext(basename(sub("[?#].*$","",cdf_file)))
  if(!is.character(sample_name)||length(sample_name)!=1L||is.na(sample_name)||!nzchar(trimws(sample_name)))
    stop("sample_name must be a nonempty string.",call.=FALSE)
  if(!is.null(xlsx_file)) GCMS_xlsx_path(xlsx_file)
  GCMS_range(calibration_range,"calibration_range")
  if(!is.character(msp_file)||length(msp_file)!=1L||is.na(msp_file)||
      !grepl("\\.msp$",msp_file,ignore.case=TRUE)) stop("msp_file must end in .msp.",call.=FALSE)
  if(!is.list(psf_control)||length(psf_control)&&(is.null(names(psf_control))||
      any(!nzchar(names(psf_control)))||anyDuplicated(names(psf_control)))) stop("Use a named psf_control list.",call.=FALSE)
  if(any(c("y","nwindows") %in% names(psf_control))) stop("Pass data and nwindows as GCMSdeco arguments.",call.=FALSE)
  timed<-function(label,expression) {
    seconds[label]<<-system.time(value<-force(expression))[["elapsed"]]
    if(verbose){cat(label,": ",round(seconds[label],2)," s\n",sep="");flush.console()}
    value
  }
  raw<-timed("read_cdf",GCMS_read_cdf(cdf_file,time_range));Y<-raw$Y;n<-nrow(Y)
  if(n<9L||nwindows>n) stop("Use at least nine scans and at most one PSF window per scan.",call.=FALSE)
  spec<-slcn_family(family,dispersion,theta)
  censor<-if(!is.null(control$censor)) control$censor else
    if(is.list(family)&&!is.null(family$scalar)) family$scalar$censor else .5
  noise<-timed("noise",slcn_noise(Y,spec,1,ncores,correlation,censor));spec<-noise$spec
  factor<-spec$information_factor
  # A shared Poisson working weight gives the TIC pilot the same average
  # variance scale as the final likelihood, including its covariance adjustment.
  if(spec$name %in% c("poisson","quasipoisson","nb1","Gamma_scale")) {
    working_weight<-1/(spec$dispersion*factor)
  } else {
    level<-if(is.null(noise$pilot)) slcn_pilot(Y) else noise$pilot
    total<-mean(rowSums(level))
    variance<-switch(spec$name,gaussian=ncol(Y)*spec$dispersion,
      nb=mean(rowSums(level+level^2/spec$theta)),Gamma=spec$dispersion*mean(rowSums(level^2)))
    working_weight<-total/max(variance,1e-15)/factor
  }
  calibration<-NULL;initial<-psf_control$kernel_est
  if(!is.null(alkanes_file)) {
    std<-timed("read_standard",if(identical(alkanes_file,cdf_file)&&identical(calibration_range,time_range))
      raw else GCMS_read_cdf(alkanes_file,calibration_range))
    calibration<-timed("calibrate_psf",estimate_psf_from_std(rowSums(std$Y),time=std$time,
      output_timepoints=raw$time,n_peaks=16L,return_details=TRUE))
    if(is.null(initial)) initial<-GCMS_calibration_start(calibration,raw$time,nwindows,psf_control$L)
  }
  pilot_df<-max(4L,min(40L,3L+ceiling(n/(200*15))))
  args<-utils::modifyList(list(y=rowSums(Y),family=poisson("identity"),kernel_est=initial,
    L=if(is.matrix(initial)) ncol(initial) else if(!is.null(initial)) length(initial) else 41L,
    baseline="bspline",nknots=pilot_df,nwindows=nwindows,plot=FALSE,ncores=ncores,
    IC=if(is.character(lambda)) lambda else "gic",
    weights=working_weight),psf_control)
  if(nwindows>1L&&length(args$weights)>1L)
    stop("Use a scalar PSF working weight with multiple windows.",call.=FALSE)
  working_weight<-args$weights
  blind<-timed("blind_psf",do.call(L0glm_blind_deco,args))
  X<-methods::as(methods::as(Matrix::Matrix(blind$X,sparse=TRUE),"generalMatrix"),"dgCMatrix")
  # Interpolated PSFs are projected on their local support, retaining structural zeros.
  X<-timed("project_psf",slcn_project_design_cpp(X,ncores))
  blind$X<-blind$X_spikes<-X
  sample<-unique(round(seq(1,n,length.out=min(n,21L))))
  widths<-vapply(sample,function(j){h<-X[,j];sum(h>=max(h)/2)},0)
  width<-max(1,stats::median(widths))
  if(is.null(baseline_df)) baseline_df<-max(4L,min(40L,3L+ceiling(n/(200*width))))
  L0glm_deco_length(baseline_df,"baseline_df",4L)
  if(baseline_df>n) stop("baseline_df exceeds the number of scans.",call.=FALSE)
  B<-as.matrix(splines::bs(raw$time,df=baseline_df,intercept=TRUE))
  model<-timed("factorization",splogconNMF(X,Y,family=family,lambda=lambda,background=B,precluster=precluster,
    dispersion=spec$dispersion,theta=if(spec$name=="nb") spec$theta*(if(spec$theta_update) -1 else 1) else NULL,
    correlation=spec$correlation,ncores=ncores,control=control,verbose=verbose))
  model$spec$dispersion_estimated<-spec$dispersion_estimated
  model$spec$theta_estimated<-spec$theta_estimated
  model$spec$correlation_estimate<-spec$correlation_estimate
  model$spec$noise_diagnostics<-spec$noise_diagnostics
  model$nuisance_df<-as.integer(spec$dispersion_estimated)+as.integer(spec$theta_estimated)
  model$ic_df<-ncol(Y)*model$rank+model$nuisance_df
  compounds<-GCMS_compounds(model,raw$time,sample_name)
  compounds<-timed("write_msp",GCMS_write_msp(model,raw$time,raw$mass,msp_file,compounds))
  compounds$apex_scan<-raw$scans[compounds$apex_scan]
  if(!is.null(xlsx_file)) timed("write_xlsx",GCMS_write_xlsx(model,raw$time,raw$mass,
    compounds,xlsx_file,cdf_file,alkanes_file,seconds))
  model$y<-model$mu<-NULL
  if(!is.null(model$aggregation$model)) model$aggregation$model$mu<-NULL
  blind$fitted_values<-blind$fitted_baseline<-NULL
  seconds["total"]<-proc.time()[["elapsed"]]-begun
  structure(list(model=model,compounds=compounds,psf=blind,background_df=baseline_df,
    psf_standard=calibration,psf_width_scans=width,time=raw$time,mass=raw$mass,
    scans=raw$scans,psf_working_weight=working_weight,
    cdf_file=if(grepl("^https?://",cdf_file)) cdf_file else normalizePath(cdf_file,winslash="/"),
    sample_name=sample_name,msp_file=normalizePath(msp_file,winslash="/"),
    xlsx_file=if(is.null(xlsx_file)) NULL else normalizePath(xlsx_file,winslash="/"),
    seconds=seconds),class="GCMSdeco")
}

#' @rdname GCMSdeco
#' @export
coef.GCMSdeco <- function(object,...,normalize=c("none","l1")) {
  normalize<-match.arg(normalize);S<-stats::coef(object$model,...)
  if(normalize=="l1"&&nrow(S)) S<-S/rowSums(S)
  S
}
#' @rdname GCMSdeco
#' @export
fitted.GCMSdeco <- function(object,...) stats::fitted(object$model,...)
#' @rdname GCMSdeco
#' @export
nobs.GCMSdeco <- function(object,...) stats::nobs(object$model,...)
#' @rdname GCMSdeco
#' @export
logLik.GCMSdeco <- function(object,...) stats::logLik(object$model,...)
#' @rdname GCMSdeco
#' @export
AIC.GCMSdeco <- function(object,...,k=2) stats::AIC(object$model,...,k=k)
#' @rdname GCMSdeco
#' @export
BIC.GCMSdeco <- function(object,...) stats::BIC(object$model,...)
#' @rdname GCMSdeco
#' @export
GIC.GCMSdeco <- function(object,...) GIC(object$model,...)
#' @rdname GCMSdeco
#' @export
EBIC.GCMSdeco <- function(object,...) EBIC(object$model,...)
#' @rdname GCMSdeco
#' @export
MBIC.GCMSdeco <- function(object,...) MBIC(object$model,...)
#' @rdname GCMSdeco
#' @method print GCMSdeco
#' @export
print.GCMSdeco <- function(x,...) {
  cat("GC/MS deconvolution: ",x$model$rank," compounds; ",length(x$mass),
    " mass channels; ",length(x$time)," scans\n",
    "GIC = ",format(GIC(x),digits=9),"; ",round(x$seconds[["total"]],2)," s\n",
    "Spectra: ",x$msp_file,"\n",sep="")
  invisible(x)
}
