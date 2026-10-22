sparseglm_plot_fixture <- function() {
  truth <- Matrix::Matrix(rbind(c(2,0),c(0,3),c(0,0),c(4,1),c(0,0)),sparse=TRUE)
  beta <- Matrix::Matrix(rbind(c(2,0),c(0,0),c(1,0),c(5,1),c(0,0)),sparse=TRUE)
  model <- function(B) structure(list(coefficients=B),class="L0glm_fit")
  structure(list(fits=list(first=model(beta),second=model(truth),missing=NULL),
    beta_true=truth,timings=data.frame(seconds=c(1.2,2.3,NA_real_)),
    data=list(intercept=FALSE,y=as.matrix(truth)+.25)),class="sparseglm")
}

with_plot_device <- function(code) {
  file <- tempfile(fileext=".pdf")
  grDevices::pdf(file)
  on.exit({grDevices::dev.off();unlink(file)})
  force(code)
}

test_that("coefficient plots avoid repeated joint zeros and retain sparse indices", {
  object <- sparseglm_plot_fixture()
  beta <- coef(object)$first
  points <- L0glm_plot_coefficients(beta,object$beta_true,"all")
  expect_equal(sum(points$estimated == 0 & points$truth == 0),1L)
  expect_equal(length(points$index),6L)
  expect_equal(points$estimated,as.numeric(beta[points$index]))
  expect_equal(points$truth,as.numeric(object$beta_true[points$index]))
  expect_equal(sort(points$color),sort(c("navy","navy","magenta","red","navy","navy")))
  huge <- Matrix::sparseMatrix(i=c(2L,1000000L),j=c(1L,3L),x=c(3,4),dims=c(1000000L,3L))
  sparse <- L0glm_plot_coefficients(huge,huge,"all")
  expect_equal(length(sparse$estimated),3L)
  expect_equal(sum(sparse$estimated == 0 & sparse$truth == 0),1L)
  zero <- Matrix::Matrix(0,5,2,sparse=TRUE)
  expect_equal(length(L0glm_plot_coefficients(zero,zero,"all")$index),1L)
  expect_length(L0glm_plot_coefficients(zero,NULL,"all")$index,0L)
})

test_that("row mean plots retain support classification before signed cancellation", {
  beta <- Matrix::Matrix(rbind(c(1,-1),c(1,2),c(0,0)),sparse=TRUE)
  truth <- Matrix::Matrix(rbind(c(0,0),c(2,2),c(3,-1)),sparse=TRUE)
  points <- L0glm_plot_coefficients(beta,truth,"mean")
  expect_equal(points$estimated,as.numeric(Matrix::rowMeans(beta))[points$index])
  expect_equal(points$truth,as.numeric(Matrix::rowMeans(truth))[points$index])
  expect_equal(points$color[points$index == 1L],"red")
  expect_equal(points$color[points$index == 3L],"magenta")
})

test_that("method selection preserves requested order and graphics state", {
  object <- sparseglm_plot_fixture()
  titles <- character()
  layouts <- list()
  testthat::local_mocked_bindings(plot.default=function(x,y,...) {
    titles <<- c(titles,list(...)$main)
    layouts[[length(layouts)+1L]] <<- graphics::par("mfrow")
    graphics::plot.new()
  }, .package="graphics")
  with_plot_device({
    original <- graphics::par("mfrow")
    expect_identical(plot(object,which=c(2,1),channels="mean"),object)
    expect_equal(titles,c("second (2.3 s)","first (1.2 s)"))
    expect_equal(layouts[[1L]],c(2L,1L))
    expect_equal(graphics::par("mfrow"),original)
    expect_warning(plot(object,which=c(3,1)),"missing")
    for (which in list(0,4,1.5,NA_real_,Inf,TRUE))
      expect_error(plot(object,which=which),"whole-number")
    expect_error(plot(object,type="unknown"),"arg")
    expect_error(plot(object,channels="unknown"),"arg")
    expect_error(plot(object,type="heatmap",channels="mean"),"coefficient plots")
    absent <- object; absent$beta_true <- NULL
    expect_error(plot(absent,type="heatmap"),"beta_true")
    single <- object
    single$fits <- lapply(single$fits,function(fit) if(!is.null(fit)) {
      fit$coefficients <- fit$coefficients[,1,drop=FALSE];fit
    })
    expect_error(plot(single,type="heatmap"),"multichannel")
    plot(object,which=c(1,2,1))
    expect_equal(tail(layouts,1L)[[1L]],c(3L,1L))
  })
})

test_that("million-predictor paths draw only varying paths and one zero", {
  B <- Matrix::sparseMatrix(i=c(2L,999999L),j=c(1L,1L),x=c(2,-3),
    dims=c(1000000L,1L))
  object <- list(coefficients=list(B,B/2),always.include=integer(),
    metrics=data.frame(lambda=c(20,2),support_size_estimated=c(2,2),gic=c(5,6)),
    tune="gic",lambda_best=20,lambda_preset=20,lambda_max=20)
  plotted <- NULL
  testthat::local_mocked_bindings(matplot=function(x,y,...) {
    plotted <<- y
    graphics::plot(x,range(y),type="n")
  }, .package="graphics")
  with_plot_device(L0glm_plot_tuning(object))
  expect_equal(dim(plotted),c(2L,3L))
  expect_equal(sort(colSums(abs(plotted))),c(0,3,4.5))
})

test_that("heatmaps stack responses, truth and methods with correct support lines", {
  object <- sparseglm_plot_fixture()
  calls <- lines <- list()
  testthat::local_mocked_bindings(plotmat=function(M,...) {
    calls[[length(calls)+1L]] <<- c(list(M=M,layout=graphics::par("mfrow")),list(...))
    graphics::plot.new();graphics::plot.window(xlim=c(1,nrow(M)),ylim=c(1,ncol(M)))
  }, .package="L0glm")
  testthat::local_mocked_bindings(abline=function(...) {
    lines[[length(lines)+1L]] <<- list(...)
  }, .package="graphics")
  with_plot_device({
    plot(object,type="heatmap")
    expect_equal(length(calls),3L)
    expect_equal(vapply(calls,`[[`,character(1),"main"),c("Observed data","Ground truth","first (1.2 s)"))
    expect_equal(calls[[1L]]$M,object$data$y,ignore_attr=TRUE)
    expect_equal(calls[[1L]]$layout,c(3L,1L))
    expect_equal(calls[[2L]]$zlim,calls[[3L]]$zlim)
    expect_equal(lines[[1L]]$v,c(1,2,4))
    expect_equal(lines[[2L]]$v,c(1,4))
    expect_equal(lines[[3L]]$v,3)
    expect_equal(lines[[4L]]$v,2)
    expect_equal(vapply(lines,`[[`,character(1),"col"),
      grDevices::adjustcolor(c("cyan","cyan","red","magenta"),alpha.f=.8))
    calls <- lines <- list()
    plot(object,type="heatmap",which=c(2,1),maxresolution=2)
    expect_equal(vapply(calls,`[[`,character(1),"main"),
      c("Observed data","Ground truth","second (2.3 s)","first (1.2 s)"))
    expect_equal(calls[[1L]]$layout,c(4L,1L))
    expect_equal(calls[[2L]]$zlim,calls[[4L]]$zlim)
    expect_equal(lines[[1L]]$v,1+(c(1,2,4)-1)/4)
    # An added intercept is not a time-index row of the images.
    object$data$intercept <- TRUE
    calls <- lines <- list()
    plot(object,type="heatmap",which=1)
    expect_equal(nrow(calls[[1L]]$M),5L)
    expect_equal(nrow(calls[[2L]]$M),4L)
    expect_equal(lines[[1L]]$v,c(1,3))
    # Dense support overlays remain visible without hiding the image itself.
    object$data$intercept <- FALSE
    object$data$y <- NULL
    object$beta_true <- Matrix::Matrix(1,1000,2,sparse=TRUE)
    object$fits <- list(first=structure(list(coefficients=object$beta_true),class="L0glm_fit"))
    calls <- lines <- list()
    plot(object,type="heatmap")
    expect_equal(lines[[1L]]$col,grDevices::adjustcolor("cyan",alpha.f=.25))
    expect_length(lines[[1L]]$v,1000L)
  })
})

test_that("signed heatmaps preserve negatives and one-dimensional image axes work", {
  captured <- NULL
  testthat::local_mocked_bindings(image=function(...) {captured <<- list(...)}, .package="graphics")
  M <- matrix(c(-4,0,1,9),2)
  plotmat(M,gamma=.5,signed=TRUE,axes=FALSE)
  expect_equal(unname(captured$z),matrix(c(-2,0,1,3),2))
  plotmat(M,gamma=.5,axes=FALSE)
  expect_equal(unname(captured$z),matrix(c(0,0,1,3),2))
  testthat::local_mocked_bindings(axis=function(...) NULL,box=function(...) NULL,.package="graphics")
  expect_no_error(plotmat(matrix(0,1,3)))
})
