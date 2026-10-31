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
#' @param nrows,ncols Rows per PDF page, and columns for spectral panels.
#' @param minutes_per_panel Retention-time span of each line panel, default 2.5
#'   minutes. Consecutive panels cover the whole selected interval.
#' @param transform Display transform: signed square root (default), identity
#'   or signed \code{log1p}. Negative background-corrected values and residuals
#'   are retained. Line-axis tick labels use the original intensity scale.
#' @param peak_scale \code{"basepeak"} scales each profile by its strongest ion
#'   in the displayed mass range, matching the individual grey mass traces.
#'   \code{"total"} uses the compound's summed ion intensity and shows the
#'   observed total ion trace in grey.
#' @param labels Label line-profile apices by compound index.
#' @param max_time Maximum heatmap time bins. Data and fitted signals are averaged
#'   over the same bins, retaining all selected mass channels.
#' @param clip Upper quantile used for heatmap colour saturation, default .999.
#'   Set to one to display the full transformed range.
#' @param ... Graphical arguments passed to the line or spectral panel.
#' @return Invisibly, \code{x}. Plotting does not attach observed data to it.
#' @details Heatmaps show the observed signal, observed signal minus modeled
#'   background, fitted chromatogram and signed residuals. Spline background and
#'   any known offset are subtracted in the second heatmap and grey line traces.
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
                         max_time=1600L,clip=.999,...) {
  type<-match.arg(type);transform<-match.arg(transform);peak_scale<-match.arg(peak_scale)
  for(name in c("nrows","ncols","max_time")) L0glm_deco_length(get(name),name,1L)
  if(length(minutes_per_panel)!=1L||!is.finite(minutes_per_panel)||minutes_per_panel<=0||
      length(clip)!=1L||!is.finite(clip)||clip<=0||clip>1) stop("Invalid plotting span or clip quantile.",call.=FALSE)
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
    t<-as.numeric(Q%*%x$time[rows]);Y<-as.matrix(Q%*%Y)
    background<-as.matrix((Q%*%B[rows,,drop=FALSE])%*%A[,cols,drop=FALSE])
    O<-GCMS_plot_offset(x,rows,cols)
    background<-background+if(length(O)==1L) O else as.matrix(Q%*%O)
    mu<-as.matrix((Q%*%H[rows,,drop=FALSE])%*%S[,cols,drop=FALSE])+background
    graphics::par(mfrow=c(4L,1L),mar=c(3,4,2.2,2),mgp=c(1.8,.6,0))
    titles<-c("Observed chromatogram","Observed minus modeled background","Fitted chromatogram","Residual: observed minus fitted")
    for(panel in 1:4) {
      z<-trans(switch(panel,Y,Y-background,mu,Y-mu))
      bound<-as.numeric(stats::quantile(abs(z),clip,names=FALSE));bound<-max(bound,1e-12)
      signed<-panel%in%c(2L,4L)&&any(z<0)
      palette<-if(signed) grDevices::colorRampPalette(c("#ae267e","#f9fafb","#123e68"))(201L) else
        grDevices::colorRampPalette(c("#ffffff","#9de2db","#087e83","#102d4d"))(201L)
      z[]<-pmax(if(signed) -bound else 0,pmin(bound,z))
      graphics::image(seq(min(t),max(t),length.out=length(t)),x$mass[cols],z,
        col=palette,zlim=if(signed) c(-bound,bound) else c(0,bound),
        useRaster=length(cols)<3L||length(unique(diff(x$mass[cols])))==1L,
        xlab="Retention time (min)",ylab="m/z",main=paste(sample,titles[panel],sep=" | "))
    }
  } else {
    graphics::par(mfrow=c(nrows,1L),mar=c(3,4,2.2,1),mgp=c(1.8,.6,0))
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
        xlab="Retention time (min)",ylab="Intensity",yaxt="n",ylim=ylim,
        main=sprintf("%s | %s display | %s-scaled compound profiles",sample,transform,peak_scale),...)
      ticks<-pretty(ylim,4L)
      graphics::axis(2,at=ticks,labels=format(signif(inverse(ticks),3),trim=TRUE),las=1)
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
