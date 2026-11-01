test_that("ANDI-MS matrix exports preserve times, intensities and silent mass channels", {
  skip_if_not_installed("ncdf4")
  set.seed(27);Y <- matrix(rpois(15L*5L,3),15L,5L);Y[,3L] <- 0
  time <- cumsum(c(1,seq(.01,.02,length.out=14)))
  mass <- 40:44;dimnames(Y) <- list(as.character(time),as.character(mass))
  for(scale in c(1,.125,1e9)) {
    file <- tempfile(fileext=".cdf");on.exit(unlink(file),add=TRUE)
    write_GCMS_cdf(Y*scale,file,sample_name="Matrix example")
    a <- GCMS_read_cdf(file)
    expect_equal(unname(a$Y),unname(Y*scale),tolerance=0)
    expect_equal(a$time,time,tolerance=1e-12);expect_identical(a$mass,mass)
    nc <- ncdf4::nc_open(file)
    expect_equal(ncdf4::ncatt_get(nc,0,"dataset_completeness")$value,"C1+C2")
    expect_equal(ncdf4::ncatt_get(nc,0,"experiment_title")$value,"Matrix example")
    expect_equal(ncdf4::ncatt_get(nc,"scan_acquisition_time","units")$value,"seconds")
    expect_equal(as.numeric(ncdf4::ncvar_get(nc,"total_intensity")),unname(rowSums(Y*scale)),tolerance=0)
    expect_true(nc$dim$point_number$unlim);expect_length(nc$var,20L)
    expect_equal(as.numeric(ncdf4::ncvar_get(nc,"point_count")),rep(5L,15L))
    ncdf4::nc_close(nc)
  }
  file <- tempfile(fileext=".cdf")
  expect_error(write_GCMS_cdf(Y,file,time=rev(time)),"scan times")
  expect_error(write_GCMS_cdf(Y,file,mass=rep(1,5)),"integer m/z")
  expect_error(write_GCMS_cdf(-Y,file),"nonnegative")
  expect_error(write_GCMS_cdf(Y,"file.txt"),"end in")
  expect_error(write_GCMS_cdf(Y,file,sample_name=""),"sample_name")
})

test_that("short CDF windows can use a constant unpenalized background", {
  skip_if_not_installed("ncdf4")
  set.seed(278);n <- 70L
  kernel <- convkernel("gaussian",11L,list(sigma=1.5))
  X <- convmatrix(kernel,n)
  mean <- as.matrix(X[,c(20L,48L)]%*%matrix(c(80,40,20,100,50,10),2L,3L))+2
  Y <- matrix(rpois(length(mean),mean),n,3L)
  cdf <- tempfile(fileext=".cdf");msp <- tempfile(fileext=".msp")
  on.exit(unlink(c(cdf,msp)))
  write_GCMS_cdf(Y,cdf,time=seq_len(n)/60,mass=40:42)
  f <- GCMSdeco(cdf,msp,nwindows=1L,baseline_df=1L,precluster=2L,
    family=poisson("identity"),correlation="none",ncores=1L,
    psf_control=list(kernel_est=kernel,maxiter=0L),
    control=list(profile_cycles=0L,profile_fusion=FALSE))
  expect_equal(f$background_df,1L)
  expect_equal(f$model$background,matrix(1,n,1L))
  expect_equal(ncol(coef(f)),3L)
  expect_equal(dim(fitted(f)),dim(Y));expect_true(file.exists(msp))
  expect_lt(abs(f$model$ledger_error),1e-4)
  expect_identical(f$psf$iterations,0L)
  expect_equal(f$psf$g_hat,kernel/sum(kernel),tolerance=1e-12)
  expect_identical(f$psf$method,"fixed")
  expect_equal(dim(f$psf$X),c(n,n))
  expect_error(L0glm_blind_deco(rowSums(Y),maxiter=0L,plot=FALSE),"kernel_est")
  windows <- L0glm_blind_deco(rowSums(Y),kernel_est=kernel,nwindows=2L,
    maxiter=0L,plot=FALSE,baseline="bspline",nknots=4L)
  expect_identical(windows$iterations,c(0L,0L))
  expect_equal(unname(windows$g_hat),rbind(kernel/sum(kernel),kernel/sum(kernel)),tolerance=1e-12)
  expect_true(windows$converged)
})
