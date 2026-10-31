GCMS_plot_selection <- function(x,time_range,mass) {
  time_range<-GCMS_range(time_range,"time_range")
  rows<-if(is.null(time_range)) seq_along(x$time) else which(x$time>=time_range[1L]&x$time<=time_range[2L])
  cols<-if(is.null(mass)) seq_along(x$mass) else match(mass,x$mass)
  if(!length(rows)||!length(cols)||anyNA(cols)) stop("Select recorded retention times and masses.",call.=FALSE)
  list(rows=rows,cols=sort(unique(cols)))
}

GCMS_plot_observed <- function(x,rows,cols) {
  raw<-GCMS_read_cdf(x$cdf_file,range(x$time[rows])+c(-1,1)*1e-10,mass=x$mass[cols])
  target<-x$time[rows];left<-pmax(1L,findInterval(target,raw$time))
  right<-pmin(left+1L,length(raw$time))
  ix<-ifelse(abs(raw$time[left]-target)<=abs(raw$time[right]-target),left,right)
  if(any(abs(raw$time[ix]-target)>1e-9)) stop("The source CDF scan grid differs from the fitted sample.",call.=FALSE)
  raw$Y[ix,,drop=FALSE]
}

GCMS_plot_offset <- function(x,rows,cols) {
  O<-x$model$offset
  if(length(O)==1L) O else O[rows,cols,drop=FALSE]
}

GCMS_pearson <- function(x,Y,mu,rows,cols) {
  spec<-x$model$spec;phi<-spec$dispersion;theta<-spec$theta
  variance<-switch(spec$name,gaussian=matrix(phi,nrow(mu),ncol(mu)),poisson=mu,
    quasipoisson=phi*mu,nb=mu+mu^2/theta,nb1=phi*mu,
    Gamma=phi*mu^2,Gamma_scale=phi*mu)
  expected<-mu
  detection<-x$model$control$censor
  if(is.null(detection)) detection<-0
  if(spec$name%in%c("Gamma","Gamma_scale")&&detection>0) {
    positive<-which(mu>0)
    a<-if(spec$name=="Gamma") rep(1/phi,length(positive)) else mu[positive]/phi
    b<-if(spec$name=="Gamma") phi*mu[positive] else rep(phi,length(positive))
    expected[]<-variance[]<-0
    expected[positive]<-mu[positive]*stats::pgamma(detection,a+1,scale=b,lower.tail=FALSE)
    second<-a*(a+1)*b^2*stats::pgamma(detection,a+2,scale=b,lower.tail=FALSE)
    variance[positive]<-pmax(0,second-expected[positive]^2)
  }
  W<-x$model$weights
  if(is.null(W)) W<-1
  if(length(W)>1L) W<-W[rows,cols,drop=FALSE]
  delta<-Y-expected
  residual<-sqrt(W)*delta/sqrt(variance)
  residual[variance==0&delta==0]<-0
  residual[!is.finite(W)|W==0]<-NA_real_
  residual
}

#' @rdname GCMSdeco
#' @param type Residual type: \code{"pearson"} uses the fitted marginal variance;
#'   \code{"response"} returns observed minus fitted latent intensity.
#' @param mass Optional nominal m/z values for residual extraction.
#' @export
residuals.GCMSdeco <- function(object,type=c("pearson","response"),
                               time_range=NULL,mass=NULL,...) {
  type<-match.arg(type);selection<-GCMS_plot_selection(object,time_range,mass)
  rows<-selection$rows;cols<-selection$cols;model<-object$model
  Y<-GCMS_plot_observed(object,rows,cols)
  mu<-as.matrix(model$H[rows,,drop=FALSE]%*%model$S[,cols,drop=FALSE])+
    model$background[rows,,drop=FALSE]%*%model$alpha[,cols,drop=FALSE]+GCMS_plot_offset(object,rows,cols)
  if(type=="response") Y-mu else GCMS_pearson(object,Y,mu,rows,cols)
}

#' Plot GC/MS signals, compound profiles and recovered spectra
#'
#' Draw four signal heatmaps, background-corrected observed mass traces with
#' fitted compound profiles, or recovered mass spectra. Observations are read
#' from the source CDF for the requested interval and released after plotting.
#'
#' @param x A fitted [GCMSdeco()] object.
#' @param type \code{"lines"}, \code{"heatmap"} or \code{"spectra"}.
#' @param file Optional PDF output path. \code{NULL} uses the current device.
#'   PDF reports can contain multiple pages.
#' @param time_range Retention-time bounds in minutes; \code{NULL} uses the full
#'   fitted interval.
#' @param mass Optional vector of recorded nominal m/z values.
#' @param which Compound indices, as returned in \code{x$compounds$compound}.
#'   \code{NULL} shows all compounds.
#' @param order Spectral panel order: \code{"rt"} (default) uses increasing apex
#'   retention time; \code{"area"} or \code{"height"} uses decreasing compound
#'   area or total-ion peak height; \code{"index"} preserves the selected indices.
#' @param nrows,ncols Rows per PDF page, and columns for spectral panels.
#' @param minutes_per_panel Retention-time span of each line panel, default 2.5
#'   minutes. Consecutive panels cover the whole selected interval.
#' @param transform Line-display transform: signed square root (default),
#'   identity or signed \code{log1p}. Negative background-corrected traces are
#'   retained. Tick labels use the original intensity scale.
#' @param peak_scale \code{"basepeak"} scales each profile by its strongest ion
#'   in the displayed mass range, matching the individual grey mass traces.
#'   \code{"total"} uses the compound's summed ion intensity and shows the
#'   observed total ion trace in grey.
#' @param labels Label line-profile apices by compound index.
#' @param max_time Maximum heatmap time bins. Data and fitted signals are averaged
#'   over the same bins, retaining all selected mass channels.
#' @param clip Upper quantile used for heatmap colour saturation, default .999.
#'   Set to one to display the full transformed range.
#' @param gamma Positive heatmap contrast exponent, default .13 as in [plotmat()].
#' @param residual_limit Symmetric Pearson-residual colour limit, default 3.
#' @param ... Graphical arguments passed to the line or spectral panel.
#' @return Invisibly, \code{x}. Plotting does not attach observed data to it.
#' @details Heatmaps show the observed signal, observed signal minus modeled
#'   background, fitted chromatogram and Pearson residuals. Spline background and
#'   any known offset are subtracted in the second heatmap and grey line traces.
#'   Signal heatmaps use the [plotmat()] black-blue-green-yellow-orange-red
#'   palette and truncate negative corrected intensities at zero. Pearson
#'   residuals use blue-black-red for negative-zero-positive deviations, without
#'   an intensity transform. Display time bins contain their average residual.
#'   Spectral panels use base-peak percentages for readability; \code{coef(x,
#'   normalize="l1")} and the Excel export provide unit-L1 spectra.
#' @seealso [GCMSdeco()], [splogconNMF()]
#' @examples
#' \dontrun{
#' fit <- GCMSdeco("beer.cdf", "beer.msp", xlsx_file="beer.xlsx")
#' plot(fit, type="heatmap", file="beer_heatmaps.pdf")
#' plot(fit, type="lines", file="beer_traces.pdf", minutes_per_panel=2.5)
#' plot(fit, type="spectra", file="beer_spectra.pdf",
#'      which=head(fit$compounds$compound[order(-fit$compounds$area)],24))
#' }
#' @export
plot.GCMSdeco <- function(x,type=c("lines","heatmap","spectra"),file=NULL,
                         time_range=NULL,mass=NULL,which=NULL,nrows=4L,ncols=2L,
                         minutes_per_panel=2.5,transform=c("sqrt","identity","log1p"),
                         peak_scale=c("basepeak","total"),labels=FALSE,
                         max_time=1600L,clip=.999,gamma=.13,residual_limit=3,
                         order=c("rt","area","height","index"),...) {
  type<-match.arg(type)
  transform<-match.arg(transform);peak_scale<-match.arg(peak_scale)
  order<-match.arg(order)
  for(name in c("nrows","ncols","max_time")) L0glm_deco_length(get(name),name,1L)
  if(length(minutes_per_panel)!=1L||!is.finite(minutes_per_panel)||minutes_per_panel<=0||
      length(clip)!=1L||!is.finite(clip)||clip<=0||clip>1||
      length(gamma)!=1L||!is.finite(gamma)||gamma<=0||
      length(residual_limit)!=1L||!is.finite(residual_limit)||residual_limit<=0)
    stop("Invalid plotting span, contrast or colour limits.",call.=FALSE)
  K<-ncol(x$model$H)
  ids<-if(is.null(which)) seq_len(K) else which
  if(!is.numeric(ids)||anyNA(ids)||any(ids!=floor(ids)|ids<1|ids>K)) stop("which must index fitted compounds.",call.=FALSE)
  selection<-GCMS_plot_selection(x,time_range,mass);rows<-selection$rows;cols<-selection$cols
  if(!is.null(file)) {
    if(!is.character(file)||length(file)!=1L||!grepl("[.]pdf$",file,ignore.case=TRUE)) stop("file must end in .pdf.",call.=FALSE)
    dir.create(dirname(file),recursive=TRUE,showWarnings=FALSE)
    grDevices::pdf(file,width=11.7,height=8.3,onefile=TRUE,useDingbats=FALSE)
    on.exit(grDevices::dev.off(),add=TRUE)
  }
  old<-graphics::par(no.readonly=TRUE);on.exit(graphics::par(old),add=TRUE,after=FALSE)
  S<-x$model$S;H<-x$model$H;B<-x$model$background;A<-x$model$alpha
  colour<-if(K) grDevices::hcl(h=(seq_len(K)*137.508)%%360,c=75,l=48) else character()
  trans<-switch(transform,sqrt=function(z) sign(z)*sqrt(abs(z)),identity=identity,
    log1p=function(z) sign(z)*log1p(abs(z)))
  inverse<-switch(transform,sqrt=function(z) sign(z)*z^2,identity=identity,
    log1p=function(z) sign(z)*expm1(abs(z)))
  sample<-if(is.null(x$sample_name)) tools::file_path_sans_ext(basename(x$cdf_file)) else x$sample_name
  if(type=="spectra") {
    if(length(ids)&&order!="index") {
      key<-switch(order,rt=vapply(ids,function(k) x$time[which.max(H[,k])],0),
        area=-rowSums(S[ids,,drop=FALSE]),
        height=-vapply(ids,function(k) max(H[,k])*sum(S[k,]),0))
      ids<-ids[base::order(key)]
    }
    graphics::par(mfrow=c(nrows,ncols),mar=c(3,3,2.5,1),mgp=c(1.8,.6,0))
    if(!length(ids)) {graphics::plot.new();graphics::title("No recovered compounds");return(invisible(x))}
    for(k in ids) {
      s<-S[k,cols];scale<-max(s)
      graphics::plot(x$mass[cols],if(scale>0) 100*s/scale else s,type="h",col=colour[k],
        xlab="m/z",ylab="Relative intensity (%)",ylim=c(0,105),
        main=sprintf("%s | compound %d | RT %.3f min",sample,k,x$time[which.max(H[,k])]),...)
    }
  } else if(type=="heatmap") {
    Y<-GCMS_plot_observed(x,rows,cols)
    bins<-pmin(max_time,ceiling(seq_along(rows)*min(max_time,length(rows))/length(rows)))
    size<-tabulate(bins);Q<-Matrix::sparseMatrix(i=bins,j=seq_along(rows),x=1/size[bins])
    t<-as.numeric(Q%*%x$time[rows])
    pearson<-matrix(0,length(size),length(cols))
    for(first in seq.int(1L,length(rows),by=512L)) {
      at<-seq.int(first,min(first+511L,length(rows)));ix<-rows[at]
      mean<-as.matrix(H[ix,,drop=FALSE]%*%S[,cols,drop=FALSE])+
        B[ix,,drop=FALSE]%*%A[,cols,drop=FALSE]+GCMS_plot_offset(x,ix,cols)
      r<-GCMS_pearson(x,Y[at,,drop=FALSE],mean,ix,cols)
      r[is.na(r)]<-0
      r[is.infinite(r)]<-sign(r[is.infinite(r)])*1e100
      pearson<-pearson+as.matrix(Q[,at,drop=FALSE]%*%r)
    }
    Y<-as.matrix(Q%*%Y)
    background<-as.matrix((Q%*%B[rows,,drop=FALSE])%*%A[,cols,drop=FALSE])
    O<-GCMS_plot_offset(x,rows,cols)
    background<-background+if(length(O)==1L) O else as.matrix(Q%*%O)
    mu<-as.matrix((Q%*%H[rows,,drop=FALSE])%*%S[,cols,drop=FALSE])+background
    graphics::par(mfrow=c(4L,1L),mar=c(3,4,2.2,2),mgp=c(1.8,.6,0))
    titles<-c("Observed chromatogram","Observed minus modeled background","Fitted chromatogram","Mean Pearson residual")
    signal_limit<-max(1e-12,as.numeric(stats::quantile(pmax(Y,0)^gamma,clip,names=FALSE)))
    for(panel in 1:4) {
      z<-if(panel==4L) pearson else pmax(switch(panel,Y,Y-background,mu),0)^gamma
      bound<-if(panel==4L) residual_limit else signal_limit
      signed<-panel==4L
      palette<-if(signed) grDevices::colorRampPalette(c("blue","black","red"))(1000L) else
        grDevices::colorRampPalette(c("black","black","blue","green","yellow","orange","red"))(1000L)
      z[]<-pmax(if(signed) -bound else 0,pmin(bound,z))
      graphics::image(seq(min(t),max(t),length.out=length(t)),x$mass[cols],z,
        col=palette,zlim=if(signed) c(-bound,bound) else c(0,bound),
        useRaster=length(cols)<3L||length(unique(diff(x$mass[cols])))==1L,
        xlab="Retention time (min)",ylab="m/z",main=paste(sample,titles[panel],sep=" | "))
    }
  } else {
    graphics::par(mfrow=c(nrows,1L),mar=c(3,6,2.2,1),mgp=c(1.8,.6,0))
    starts<-seq(x$time[min(rows)],x$time[max(rows)],by=minutes_per_panel)
    for(start in starts) {
      ix<-rows[x$time[rows]>=start&x$time[rows]<start+minutes_per_panel]
      if(!length(ix)) next
      Y<-GCMS_plot_observed(x,ix,cols)
      bg<-as.matrix(B[ix,,drop=FALSE]%*%A[,cols,drop=FALSE])+GCMS_plot_offset(x,ix,cols)
      observed<-if(peak_scale=="total") matrix(rowSums(Y-bg),ncol=1L) else Y-bg
      active<-ids[Matrix::colSums(H[ix,ids,drop=FALSE])>0]
      scale<-if(peak_scale=="total") rowSums(S[active,cols,drop=FALSE]) else
        if(length(active)) apply(S[active,cols,drop=FALSE],1L,max) else numeric()
      peaks<-sweep(as.matrix(H[ix,active,drop=FALSE]),2L,scale,"*")
      value<-trans(observed);ylim<-range(value,trans(peaks),0,finite=TRUE)
      if(diff(ylim)==0) ylim<-ylim+c(-1,1)
      graphics::matplot(x$time[ix],value,type="l",lty=1,col=grDevices::adjustcolor("grey65",.4),
        xlab="Retention time (min)",ylab="",yaxt="n",ylim=ylim,
        main=sprintf("%s | %s display | %s-scaled compound profiles",sample,transform,peak_scale),...)
      ticks<-pretty(ylim,4L)
      graphics::axis(2,at=ticks,labels=formatC(inverse(ticks),format="g",digits=2),las=1)
      graphics::mtext("Intensity",side=2,line=4.5)
      if(length(active)) graphics::matlines(x$time[ix],trans(peaks),lty=1,lwd=1.4,col=colour[active])
      if(isTRUE(labels)&&length(active)) for(j in seq_along(active)) {
        apex<-which.max(peaks[,j]);graphics::text(x$time[ix[apex]],trans(peaks[apex,j]),
          active[j],col=colour[active[j]],pos=3,cex=.6,xpd=NA)
      }
      graphics::abline(h=0,col="grey80",lty=3)
    }
  }
  invisible(x)
}
