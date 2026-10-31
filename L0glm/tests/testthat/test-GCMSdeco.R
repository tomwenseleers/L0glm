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
    z <- GCMSdeco(file,out,family=f,nwindows=if(i==2L) 2L else 1L,
      dispersion=if(i==1L) 16 else if(i==3L) 5 else NULL,
      theta=if(grepl("negative binomial",f$family,ignore.case=TRUE)) -10 else NULL,
      correlation="none",precluster=3L,ncores=1L,baseline_df=4L,
      psf_control=list(kernel_est=kernel,maxiter=1L),
      control=list(broad_repeats=0L,full_sweeps=0L,final_polishes=0L))
    expect_s3_class(z,"GCMSdeco")
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
