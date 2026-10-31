gcms_test_cdf <- function(Y,path,units="seconds",time=seq_len(nrow(Y)),padding=0L) {
  n <- nrow(Y);m <- ncol(Y);counts <- rowSums(Y!=0)
  starts <- c(0,head(cumsum(counts+padding),-1L))
  points <- sum(counts+padding)
  mass <- intensity <- numeric(points)
  for(i in seq_len(n)) {
    js <- which(Y[i,]!=0)
    at <- starts[i]+seq_along(js)
    mass[at] <- js+39;intensity[at] <- Y[i,js]
  }
  scan <- ncdf4::ncdim_def("scan_number","",seq_len(n),create_dimvar=FALSE)
  point <- ncdf4::ncdim_def("point_number","",seq_len(points),create_dimvar=FALSE)
  vars <- list(ncdf4::ncvar_def("scan_index","",list(scan),prec="integer"),
    ncdf4::ncvar_def("point_count","",list(scan),prec="integer"),
    ncdf4::ncvar_def("scan_acquisition_time",units,list(scan),prec="double"),
    ncdf4::ncvar_def("mass_values","m/z",list(point),prec="double"),
    ncdf4::ncvar_def("intensity_values","counts",list(point),prec="double"))
  nc <- ncdf4::nc_create(path,vars);on.exit(ncdf4::nc_close(nc))
  for(name in c("scan_index","point_count","scan_acquisition_time","mass_values","intensity_values"))
    ncdf4::ncvar_put(nc,name,switch(name,scan_index=starts,point_count=counts,
      scan_acquisition_time=time,mass_values=mass,intensity_values=intensity))
  invisible(path)
}

test_that("ANDI import preserves silent ions, scan padding and time units", {
  skip_if_not_installed("ncdf4")
  read <- getFromNamespace("GCMS_read_cdf","L0glm")
  set.seed(81);Y <- matrix(rpois(12*6,3),12,6);Y[,3] <- 0;Y[4,] <- 0
  for(units in c("seconds","minutes","milliseconds","")) {
    divisor <- switch(units,seconds=60,minutes=1,milliseconds=60000,60)
    file <- tempfile(fileext=".cdf")
    gcms_test_cdf(Y,file,units=units,padding=2L)
    z <- read(file)
    expect_equal(unname(z$Y),unname(Y),tolerance=0)
    expect_identical(z$mass,40:45)
    expect_equal(z$time,seq_len(12)/divisor,tolerance=1e-14)
    a <- read(file,c(2.9,8.1)/divisor)
    expect_equal(unname(a$Y),unname(Y[3:8,,drop=FALSE]),tolerance=0)
    expect_identical(a$scans,3:8)
    unlink(file)
  }
})

test_that("bulk native projection preserves storage, inputs and thread parity", {
  H <- Matrix::sparseMatrix(i=c(2:6,3:8),j=c(rep(1,5),rep(2,6)),
    x=c(1,5,1,4,1,2,1,3,8,1,2),dims=c(9L,2L))
  values <- H@x
  project <- getFromNamespace("slcn_project_design_cpp","L0glm")
  local <- getFromNamespace("spolo_lc","L0glm")
  a <- project(H,1L);b <- project(H,2L)
  expect_identical(a,b);expect_identical(H@x,values)
  expect_identical(a@i,H@i);expect_identical(a@p,H@p)
  for(j in 1:2) {
    at <- seq.int(H@p[j]+1L,H@p[j+1L])
    expect_identical(a@x[at],local(H@x[at])$h)
  }
})

test_that("CDF-to-MSP workflows preserve channels, scores and profile constraints", {
  skip_if_not_installed("ncdf4")
  set.seed(54);n <- 150L;m <- 6L
  kernel <- convkernel("gaussian",15L,list(sigma=1.5));X <- convmatrix(kernel,n)
  S <- matrix(rexp(3*m,rate=1/1500),3,m)
  mu <- as.matrix(X[,c(25,70,120)]%*%S)+10
  Y <- matrix(rpois(n*m,mu),n,m);Y[,3] <- 0
  file <- tempfile(fileext=".cdf");gcms_test_cdf(Y,file)
  families <- list(gaussian("identity"),poisson("identity"),quasipoisson("identity"),
    Gamma_censored(.3,"identity"),nb1(4,"identity"),Gamma_scale(4,"identity",censor=.5))
  if(requireNamespace("mgcv",quietly=TRUE)) families <- c(families,list(mgcv::nb(link="identity")))
  for(i in seq_along(families)) {
    f <- families[[i]];out <- tempfile(fileext=".msp")
    excel <- if(i==2L&&requireNamespace("openxlsx",quietly=TRUE)) tempfile(fileext=".xlsx") else NULL
    z <- GCMSdeco(file,out,family=f,nwindows=if(i==2L) 2L else 1L,
      dispersion=if(i==1L) 16 else if(i==3L) 5 else NULL,
      theta=if(grepl("negative binomial",f$family,ignore.case=TRUE)) -10 else NULL,
      correlation="none",precluster=3L,ncores=1L,baseline_df=4L,
      psf_control=list(kernel_est=kernel,maxiter=1L),
      control=list(broad_repeats=0L,full_sweeps=0L,final_polishes=0L),
      xlsx_file=excel,sample_name="Test sample")
    expect_s3_class(z,"GCMSdeco")
    expect_null(z$model$y);expect_null(z$model$mu)
    expect_equal(z$sample_name,"Test sample")
    expect_true(all(z$compounds$sample=="Test sample"))
    if(!is.null(excel)) {expect_true(file.exists(z$xlsx_file));unlink(excel)}
    expect_equal(dim(fitted(z)),c(n,m));expect_equal(ncol(coef(z)),m)
    expect_true(all(coef(z)[,3]==0))
    expect_equal(GIC(z),GIC(z$model));expect_equal(AIC(z),AIC(z$model))
    expect_equal(as.numeric(logLik(z)),z$model$loglik)
    expect_equal(nobs(z),n)
    expect_true(is.finite(GIC(z))&&file.exists(out))
    expect_output(print(z),"GC/MS deconvolution")
    expect_true(is.finite(z$psf_working_weight)&&z$psf_working_weight>0)
    if(i==2L) expect_equal(z$psf_working_weight,1)
    if(i==3L) expect_equal(z$psf_working_weight,.2)
    lines <- readLines(out,warn=FALSE)
    names <- which(startsWith(lines,"Name:"))
    expect_length(names,z$model$rank)
    counts <- as.integer(sub("Num Peaks: ","",grep("^Num Peaks:",lines,value=TRUE)))
    expect_identical(counts,as.integer(z$compounds$n_ions))
    for(k in seq_along(names)) {
      begin <- which(startsWith(lines,"Num Peaks:"))[k]+1L
      if(!counts[k]) next
      pairs <- scan(text=gsub(";"," ",paste(lines[seq.int(begin,length.out=counts[k])],collapse=" ")),
        what=numeric(),sep="",quiet=TRUE)
      expect_length(pairs,2L*counts[k])
      pairs <- matrix(pairs,ncol=2L,byrow=TRUE)
      s <- coef(z)[z$compounds$compound[k],];keep <- which(s>0)
      expect_equal(pairs[,1],z$mass[keep])
      expect_equal(pairs[,2],unname(999*s[keep]/max(s[keep])),tolerance=1e-7)
    }
    unlink(out)
  }
  expect_error(GCMSdeco(file,file),"msp_file")
  expect_error(GCMSdeco(file,tempfile(fileext=".msp"),nwindows=0),"nwindows")
  unlink(file)
})

gcms_test_object <- function() {
  set.seed(712);n<-70L;m<-6L;t<-seq_len(n)/60
  H<-outer(seq_len(n),c(15,35,55),function(i,j) exp(-.5*((i-j)/2)^2))
  H<-Matrix::Matrix(sweep(H,2L,colSums(H),"/"),sparse=TRUE)
  S<-matrix(rexp(3*m,rate=1/500),3,m);colnames(S)<-as.character(40:45)
  B<-matrix(1,n,1L);A<-matrix(3,1,m);mu<-as.matrix(H%*%S)+B%*%A
  Y<-matrix(rpois(n*m,mu),n,m)
  file<-tempfile(fileext=".cdf");gcms_test_cdf(Y,file)
  model<-list(H=H,S=S,background=B,alpha=A,offset=0,n=n,p=n,rank=3L,
    loglik=sum(dpois(Y,mu,log=TRUE)),spec=list(name="poisson",dispersion=1,
      theta=1,dispersion_estimated=FALSE,theta_estimated=FALSE,information_factor=1))
  class(model)<-"splogconNMF"
  compounds<-GCMS_compounds(model,t,"beer_test")
  structure(list(model=model,time=t,mass=40:45,cdf_file=file,sample_name="beer_test",
    compounds=compounds,seconds=c(total=1)),class="GCMSdeco")
}

test_that("unit-L1 spectra and areas reconstruct the original fitted signal", {
  skip_if_not_installed("ncdf4")
  x<-gcms_test_object();on.exit(unlink(x$cdf_file))
  S<-coef(x);T<-coef(x,normalize="l1");area<-rowSums(S)
  expect_equal(rowSums(T),rep(1,nrow(T)))
  expect_equal(as.matrix((x$model$H%*%Matrix::Diagonal(x=area))%*%T),
    as.matrix(x$model$H%*%S),tolerance=1e-12)
  expect_equal(x$compounds$area,area[x$compounds$compound])
  q<-rep(1/60,length(x$time));q[c(1L,length(q))]<-1/120
  expected<-as.numeric(Matrix::crossprod(x$model$H,q))*area
  expect_equal(x$compounds$area_min,expected[x$compounds$compound])
  expect_equal(fitted(x),as.matrix(x$model$H%*%S)+3)
  expect_true(all(x$compounds$sample=="beer_test"));expect_null(x$model$y)
})

test_that("Excel exports preserve sample names, areas and every normalized ion", {
  skip_if_not_installed("ncdf4");skip_if_not_installed("openxlsx")
  x<-gcms_test_object();file<-tempfile(fileext=".xlsx")
  on.exit(unlink(c(x$cdf_file,file)))
  GCMS_write_xlsx(x$model,x$time,x$mass,x$compounds,file,x$cdf_file,NULL,c(fit=1))
  expect_identical(openxlsx::getSheetNames(file),c("Compounds","Spectra","Run"))
  tab<-openxlsx::read.xlsx(file,"Compounds")
  expect_equal(tab$sample,rep("beer_test",3L));expect_equal(tab$area,x$compounds$area)
  expect_equal(tab$area_min,x$compounds$area_min)
  spectra<-openxlsx::read.xlsx(file,"Spectra")
  expect_equal(unname(as.matrix(spectra[,-(1:4)])),
    unname(coef(x,normalize="l1")[x$compounds$compound,,drop=FALSE]),tolerance=1e-12)
  expect_error(GCMS_xlsx_path("bad.csv"),"xlsx_file")
})

test_that("GCMS plots reload observed intervals and create multi-panel PDFs", {
  skip_if_not_installed("ncdf4")
  x<-gcms_test_object();before<-serialize(x,NULL)
  paths<-paste0(tempfile(),c("_lines.pdf","_heatmap.pdf","_spectra.pdf"))
  on.exit(unlink(c(x$cdf_file,paths)))
  for(i in seq_along(paths)) {
    expect_invisible(plot(x,type=c("lines","heatmap","spectra")[i],file=paths[i],
      nrows=2L,ncols=2L,minutes_per_panel=.25,max_time=19L))
    expect_gt(file.info(paths[i])$size,1000)
  }
  expect_identical(serialize(x,NULL),before)
  expect_error(plot(x,which=99),"which")
  expect_error(plot(x,mass=1000),"recorded")
  z<-GCMS_plot_observed(x,10:20,c(1L,3L))
  expect_identical(dim(z),c(11L,2L));expect_null(x$model$y)
})
