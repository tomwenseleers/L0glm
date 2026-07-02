# FUNCTIONS TO DECOMPOSE FUNCTION AS A SUM OF GAUSSIAN OR LOG-CONCAVE PEAKS
# USING BACKFITTING

# 1. COUPLE OF UTILITY FUNCTIONS ####

# GAUSSIAN PEAK FUNCTION WITH MODE AT u, WIDTH AT INFLECTION POINT w (FWHM=2.355*w) AND PEAK HEIGHT h ####
gauspeak=function(x, u, w, h=1) h*exp(((x-u)^2)/(-2*(w^2)))
# EXAMPLES
plot(x=seq(0,1,length.out=1000), y=gauspeak(x=seq(0,1,length.out=1000), u=0.5, w=0.1, h=1), ylim=c(0,1), type="l")


# FUNCTION TO CALCULATE FULL WIDTH AT HALF MAXIMUM OF A FITTED SIGNAL ####
# plus corresponding width of a Gaussian peak function if signal were Gaussian
fwhm = function(x, y, interpol=FALSE) { 
  halfheight = max(y)/2
  id.maxy = which.max(y)
  y1 = y[1:id.maxy]
  y2 = y[id.maxy:length(y)]
  x1 = x[1:id.maxy]
  x2 = x[id.maxy:length(y)]
  if (interpol) {
  x.start = approx(x=y1,y=x1, xout=halfheight, method="linear")$y # use spline() if you would like spline interpolation
  x.stop = approx(x=y2,y=x2, xout=halfheight, method="linear")$y # use spline() if you would like spline interpolation
  } else { 
    x.start = x[which.min(abs(y1-halfheight))]
    x.stop = x[which.min(abs(y2-halfheight))+id.maxy]
    }
  fwhm = x.stop-x.start
  width = fwhm/(2*sqrt(2*log(2)))
  return(list(fwhm=fwhm, width=width))
}

# EXAMPLE
x = seq(0,200,length.out=1000)
y = gauspeak(x=x, u=100, w=10, h=100)
y.fwhm <- fwhm(x=x, y=y)
# $`fwhm`
# [1] 23.54934
# 
# $width
# [1] 10.00048

# FUNCTION TO CALCULATE INFORMATION CRITERIA AIC & BIC GIVING GOODNESS OF FIT ####
# SSs are calculated on given transformed scale (e.g. sqrt() variance stabilizing function for Poisson data)
IC = function (y, yhat, npars, transf = function (y) sqrt(y)) { 
  nobs = length(y)
  RSS = sum((transf(y)-transf(yhat))^2) # residual SS on sqrt() scale
  min2LL = nobs + nobs*log(2*pi) + nobs*log(RSS/nobs) # -2*logL
  AIC = min2LL + 2*npars # 2143.1
  BIC = min2LL + log(nobs)*npars # 2202.0
  return(list(AIC=AIC, BIC=BIC))
}

# FUNCTION TO SIMULATE SUM OF GAUSSIAN PEAKS WITH POISSON NOISE ####
# at given peak modes u, peak widths w (FWHM=2.355*w) and heights h
sumgauspeaks = function(x=1:1000, u=c(200, 320, 600, 500), w=c(50,50,50,200), h=c(1000, 500, 300, 1000), plot=TRUE) {
  npeaks = length(u)
  n = length(x)
  Y_nonoise = do.call(cbind, lapply(1:npeaks, function (peak) gauspeak(x, u=u[peak], w=w[peak], h=h[peak]) ))
  y_nonoise = rowSums(Y_nonoise)
  set.seed(101) # set random seed to make examples reproducible
  Y = apply(Y_nonoise, 2, function(col) rpois(n,col)) # peaks with Poisson noise
  y = rowSums(Y)
  if (plot==TRUE) {
    plot(x, y, col="blue", ylab="Signal (points=with noise, lines=no noise)")
    lines(x, y_nonoise, lwd=3, col="red")
    matplot(Y, type="p", pch=1, col="grey", add=TRUE)
    matplot(Y_nonoise, type="l", lty=1, lwd=3, add=TRUE)
  }
  return(list(Y_nonoise=Y_nonoise, y_nonoise=y_nonoise, Y=Y, y=y))
}

# FUNCTION THAT RETURNS PEAK LOCATIONS IN P SPLINE SMOOTHED SIGNAL ####
# (peak mode & nearby local minima + peak height) 
getpeaks = function(x, y, sortby="height", decreasing=TRUE, plot=TRUE) {
  # first we determine the location of all peaks based on local maxima in the signal & nearby local minima
  offs = 1E-20
  y_filt = pmax(exp(smooth.spline(x=x, y=log(y+offs))$y)-offs, 0) # smoothing spline with optimal smoothing determined using generalized cross validation (sqrt variance normalizing transformation is performed)
  # PS for peak detection we could also look for minima in 2nd derivative, cf. plot(diff(y_filt, order=2),type="l")
  # as in https://www.originlab.com/doc/Origin-Help/FindPeak-PA
  x_maxima = which(diff(sign(diff(y_filt)))==-2)
  # x_maxima = x_maxima[c(which(diff(x_maxima)>5),tail(x_maxima,1))] # join nearby peaks
  x_maxima2 = c(1,x_maxima,length(x))
  x_minima = sapply(1:(length(x_maxima2)-1), function (i) which.min(y_filt[x_maxima2[i]:x_maxima2[i+1]])+x_maxima2[i] )-1
  if (plot) {
    plot(x,y, col="blue", ylab="Signal (points=signal, line=P spline smooth)", main="Signal with P spline smooth & local maxima (red) & minima (blue)")
    lines(x,y_filt,col="red",lwd=3)
    abline(v=x_maxima,col="red")
    abline(v=x_minima,col="blue")
  }
  peakM = do.call(rbind,lapply(1:length(x_maxima), function (peaknr) c(
    mode=x_maxima[peaknr], # peak maxima + start & stop (nearby local minima)
    start=x_minima[peaknr],
    stop=x_minima[peaknr+1],
    range=x_minima[peaknr+1]-x_minima[peaknr],
    absheight=y_filt[x_maxima[peaknr]], # abs peak height
    height=y_filt[x_maxima[peaknr]]-(y_filt[x_minima[peaknr]] + 
                                       (y_filt[x_minima[peaknr+1]]-
                                          y_filt[x_minima[peaknr]])/(x_minima[peaknr+1]-x_minima[peaknr]) * 
                                       (x_maxima[peaknr]-x_minima[peaknr])), # estimated relative peak height = abs peak height - linearly interpolated baseline at peak mode = "peak prominence" 
    width=NA # for estimated Gaussian peak width w - for later use
  ) ) ) 
  sortorder = order(peakM[,sortby], decreasing=decreasing) # sort peaks from large to small or in elution order
  peakM = peakM[sortorder,]
  return(peakM)
}

# EXAMPLE
x = 1:1000
y = sumgauspeaks(x=x, u=c(200, 320, 600, 500), w=c(50,50,50,200), h=c(1000, 500, 300, 1000), plot=TRUE)$y
peakM = getpeaks(x, y, sortby="height", plot=TRUE) # peak matrix


# FAST FUNCTION TO COMPUTE B-SPLINE BASIS ####
# of degree "deg" (3=cubic), at positions "x",
# on a uniform grid with "ndx" knots between "x0" and "x1"
# Similar but not identical to bs(), and 100x faster than bs() 
# (bs() also has boundary knots, code below produces identical interior knots)
fastbs = function(x, nknots=20, degree=3) {
  x = as.vector(x) 
  x0 = min(x)
  x1 = max(x)
  x = as.matrix(x, ncol=1)
  dx = (x1 - x0) / nknots
  t = matrix(x0 + dx * ((-degree):(nknots-1)), nrow=1)
  T = matrix(1,nrow(x), ncol(x)) %*% t
  X = x %*% matrix(1, nrow(t), ncol(t))
  D = (X - T) / dx
  B = t(apply(t(cbind(matrix(0, nrow(x), ncol(x)), (D <= 1)*1)), 2, function(row) diff(row)))
  r = matrix(c(2:length(t),1),nrow=1)
  for (k in 1:degree) B = (D * B + (k + 1 - D) * B[, r]) / k
  return(B)
}
# EXAMPLES
matplot(fastbs(x=seq(0,1,length.out=100)), type="l")
library(microbenchmark)
microbenchmark(fastbs(x=seq(0,1,length.out=100)), times=100) # 1.8 ms, >100x faster than spline::bs() function
library(splines) # splines::bs() is ca 100x slower
microbenchmark(bs(x = seq(0,1,length.out=100),
                  knots = seq(0, 1, length.out=20), # 20 knots
                  degree = 3,
                  intercept = FALSE), times=100) # 168 ms


# FUNCTION TO CARRY OUT FAST GAUSSIAN PEAK FIT (USING TER BRAAK / OKSANEN PARAMETERISATION) #### 
# i.e. fit  a quadratic model on a log scale using a least square model fit on a log transformed Y scale (if control$fast==TRUE) or using a generalized linear model with Poisson noise fit on a log link scale (if control$fast==FALSE)
# see https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/0012-9658(2001)082%5B1191:CIFTOI%5D2.0.CO%3B2
gausfit = function(x, y, xout=x, control=list(start=NULL, fast=TRUE)) {
  # control_default = list(w=0.05, u=x[which.max(y)], h=max(y), wtolfac=5)
  # if (length(control)==0) { control = control_default } else { control = modifyList(control_default, control) }
  # wlim <- c(control$w/control$wtolfac, control$w*control$wtolfac)
  maxy <- max(y)
  if (maxy==0) { return(list(x=x, y=y, fitted=rep(0,length(y)))) }
  if (maxy<=1) { ymax <- 1E5 } else { ymax <- maxy }
  yscaling <- ymax/maxy
  y <- pmax(round(y*yscaling),0)
  # control$h <- control$h*yscaling
  ymax_idx <- which.max(y)
  subs_start <- suppressWarnings(max(which(y[1:ymax_idx]==0)))
  if (is.infinite(subs_start)) subs_start <- 1
  subs_stop <- (ymax_idx+suppressWarnings(min(which(y[ymax_idx:length(y)]==0))))
  if (is.infinite(subs_stop)) subs_stop <- length(y)
  subset <- subs_start:subs_stop # y!=0 #  # remove zeros on the side
  # start <- c((-control$u^2+2*control$w^2*log(control$h))/(2*control$w^2), control$u/control$w^2, -1/(2*control$w^2)) # starting values for fitted coefficients based on initial guesses for u, h & w (assuming gaussian shape, i.e. quadratic model on log link scale)
  x_subs <- x[subset]
  y_subs <- y[subset]
  # we initialize coefficients with weighted least squares model fit on log transformed Y scale
  if (is.null(control$start)) { offs <- 1E-10
                                weights <- y_subs^2/(y_subs+offs)^2 # sqrt(y_subs+offs)
                                c <- control$start <- .lm.fit(x=cbind(1,x_subs,x_subs^2)*sqrt(weights), 
                                                              y=log(y_subs+offs)*sqrt(weights))$coefficients 
  } 
  
  # now do actual quadratic GLM fit (if fast==FALSE)
  if (control$fast==FALSE) { 
  c <- suppressWarnings(glm.fit(cbind(1,x_subs,x_subs^2), y_subs, 
                                family=poisson(link=log), 
                                # or we could use gaussian(link=log) or Gamma(link=log) if we work with y_subs+1E-100
                                start=control$start, 
                                control=glm.control(epsilon = 1e-8, maxit = 50, trace = FALSE)
                                ))$coefficients # we fit quadratic model on log link scale - note: speedglm, glmnet or MatrixModels::glm4 was not faster, for some minimal GLM fitting code see https://bwlewis.github.io/GLM/, but all were slower than R's inbuilt glm.fit (with the possibl exception of glm.svd, https://github.com/bwlewis/GLM/issues/4, but not for this small scale problem)
  }
  
  if (c[[3]]<0) {
  u_fitted <- -c[[2]]/2/c[[3]] # inferred mode
  w_fitted <- sqrt(-(1/2)/c[[3]]) # inferred width
  h_fitted <- exp(c[1] - c[2]^2/4/c[3]) # inferred peak height (* yscaling)
  fitted <- gauspeak(x=x, u=u_fitted, w=w_fitted, h=h_fitted/yscaling)
  preds <- gauspeak(x=xout, u=u_fitted, w=w_fitted, h=h_fitted/yscaling)
  } else { 
    u_fitted <- mean(x) # inferred mode
    w_fitted <- diff(range(x))/2 # inferred width
    h_fitted <- 0 # inferred peak height (* yscaling)
    fitted <- rep(0,length(x))
    preds <- rep(0,length(xout))
    }
  out <- list(x=x, y=y, 
              fitted=fitted,
              xout=xout,
              preds=preds,
              fitfunc="gausfit",
              IC=IC(y=y, yhat=fitted, npars = 3), # AIC & BIC measured on sqrt() transformed scale (to take into account Poisson noise)
              # fitted.func=function (x) with(data.frame(u=u_fitted, w=w_fitted, h=h_fitted), gauspeak(x=x, u=u, w=w, h=h)),
              fitted.pars=c(u=u_fitted, w=w_fitted, h=h_fitted),
              npars=3)
  return(out)
}

# EXAMPLE
library(microbenchmark)
x = seq(0,200,length.out=1000)
pars.exp = c(u=100, w=10, h=100) # expected parameter values
y_nonoise = gauspeak(x=x, u=pars.exp["u"], w=pars.exp["w"], h=pars.exp["h"])
y = rpois(length(y_nonoise),y_nonoise)
microbenchmark(fit <- gausfit(x=x, y=y, control=list(fast=TRUE)), times=100) # 282 ?s using least squares fit on log Y scale
microbenchmark(fit <- gausfit(x=x, y=y, control=list(fast=FALSE)), times=100) # 1.1 ms using Poisson GLM fit with log link scale

xout = 1:200
plot(x,y)
#plot(xout,gausfit(x=x, y=y, xout=xout, control=list(fast=TRUE))$preds, type="l")
lines(xout,gausfit(x=x, y=y, xout=xout, control=list(fast=TRUE))$preds, type="l", col = 2)

# check bias if fit is done using LS fit on log transformed scale vs using GLM
x = seq(0,200,length.out=1000)
pars.exp = c(u=100, w=10, h=1000) # expected parameter values
nseeds = 1000
relbias.wls = matrix(nrow=nseeds, ncol=3, dimnames=list(NULL, c("u","w","h")))
relbias.glm = relbias.wls
y_nonoise = gauspeak(x=x, u=pars.exp["u"], w=pars.exp["w"], h=pars.exp["h"])
# plot(x, y_nonoise, type="l")
for (seed in 1:nseeds) {
set.seed(seed)
y = rpois(length(y_nonoise),y_nonoise)
# plot(x,y,type="l")
fit.wls = gausfit(x=x, y=y, control=list(fast=TRUE))
fit.glm = gausfit(x=x, y=y, control=list(fast=FALSE))
pars.wls = fit.wls$fitted.pars
pars.glm = fit.glm$fitted.pars
relbias.wls[seed,] = (pars.wls-(pars.exp+1E-100))/(pars.exp+1E-100)
relbias.glm[seed,] = (pars.glm-(pars.exp+1E-100))/(pars.exp+1E-100)
}
par(mfrow=c(2,3))
col="steelblue"
hist(relbias.wls[,"u"], col=col)
hist(relbias.wls[,"w"], col=col)
hist(relbias.wls[,"h"], col=col)
hist(relbias.glm[,"u"], col=col)
hist(relbias.glm[,"w"], col=col)
hist(relbias.glm[,"h"], col=col)

dev.off()
plot(x, y, type="l", col="lightgrey")
lines(x,y_nonoise,col=adjustcolor("green", alpha.f=0.5), lwd=4)
lines(x,fit.wls$fitted,col=adjustcolor("blue", alpha.f=0.5), lwd=4)
lines(x,fit.glm$fitted,col=adjustcolor("red", alpha.f=0.5), lwd=4)


# FUNCTION TO DO LOG-CONCAVE SPLINE FIT USING QUADRATIC OR LEAST ABS DEVIATION LOSS ####

logconcfit = function(x, y, xout=x, control=list(method="cobs", tau=0.5, nknots=20)) { 
  # methods: conspline / cobs / conreg / splconreg / (eilers method to be implemented, but right now not stable enough)
  
  # control_default = list(w=0.05, u=x[which.max(y)], h=max(y), wtolfac=5)
  # if (length(control)==0) { control = control_default } else { control = modifyList(control_default, control) }
  y <- pmax(y,0)
  ymax_idx <- which.max(y)
  subs_start <- suppressWarnings(max(which(y[1:ymax_idx]==0)))
  if (is.infinite(subs_start)) subs_start <- 1
  subs_stop <- (ymax_idx+suppressWarnings(min(which(y[ymax_idx:length(y)]==0))))
  if (is.infinite(subs_stop)) subs_stop <- length(y)
  subset <- subs_start:subs_stop # y!=0 #  # remove zeros on the side
  x_subs <- x[subset]
  y_subs <- y[subset]
  
  offs <- 1E-10 # 1E2
  weights <- y_subs^2/(y_subs+offs)^2 # sqrt(y_subs+offs)
  xxoutsame=FALSE
  if (length(x)==length(xout)) { if (all.equal(x,xout)) xxoutsame=TRUE }
  
  switch(control$method,
         
      conspline={
      require(ConSpline) # concave spline with quadratic loss
      fit <- conspline(y=log(y_subs+offs), x=x_subs, type=4, 
                       wt=weights+1E-100, 
                       knots=seq(min(x_subs),max(x_subs),length.out=control$nknots))
      extrapolfun = splinefun(x=x_subs, y=pmax(exp(fit$muhat)-offs,0), method = "natural") # extrapolate using natural cubic spline TO DO: could be done on log+eps scale
      fitted <- pmax(extrapolfun(x=x),0) # CHECK
      npars <- control$nknots
      if (xxoutsame) { preds <- fitted } else { preds <- pmax(extrapolfun(x=xout),0) }
      },
      
      cobs={
      require(cobs)
      fit <- cobs(x=x_subs,y=log(y_subs+offs), constraint= "concave", lambda=0, # concave spline with quantile loss
                  knots=seq(min(x_subs),max(x_subs),length.out=control$nknots),  # nr of knots
                  nknots=control$nknots, knots.add = FALSE, repeat.delete.add = FALSE,
                  w=weights, # obs weights
                  tau=control$tau) # do median or quantile regression 
      fitted <- pmax(exp(predict(fit, z=x)[,2])-offs,0)
      npars <- control$nknots
      if (xxoutsame) { preds <- fitted } else { preds <- pmax(exp(predict(fit, z=xout)[,2])-offs,0) }
      },
      
      conreg={
      require(cobs)
      fit <- conreg(x=x_subs, y=log(y_subs+offs), # concave regression with quadratic loss but not smoothly interpolated
                              w=weights, 
                              convex=FALSE, 
                              method="Duembgen06_R") 
      extrapolfun = splinefun(x=x_subs, y=pmax(exp(fit$yf)-offs,0), method = "natural") # extrapolate using natural cubic spline TO DO: could be done on log+eps scale
      fitted <- pmax(extrapolfun(x=x),0) # CHECK
      npars <- length(fit$iKnots)
      if (xxoutsame) { preds <- fitted } else { preds <- pmax(extrapolfun(x=xout),0) }
      },
      
      splconreg={
      require(cobs)
      fit <- conreg(x=x_subs, y=log(y_subs+offs), # concave regression with quadratic loss but not smoothly interpolated
                      w=weights, 
                      convex=FALSE, 
                      method="Duembgen06_R") 
      npars <- length(fit$iKnots)
      fit <- interpSplineCon(fit) # concave regression with quadratic loss & then interpolated with cubic spline between inferred kinks
      fitted <- pmax(exp(predict(fit, x)$y)-offs,0)
      if (xxoutsame) { preds <- fitted } else { fitted <- pmax(exp(predict(fit, xout)$y)-offs,0) }
      }
  )
  
  out <- list(x=x, y=y, 
              fitted=fitted,
              xout=xout,
              preds=preds,
              fitfunc="logconcfit",
              fitted.pars=c(u=x[which.max(fitted)], w=fwhm(x=x, y=fitted)$width, h=max(fitted)),
              IC=IC(y=y, yhat=fitted, npars = npars), # AIC & BIC measured on sqrt() transformed scale (to take into account Poisson noise)
              npars=npars)
  return(out)
}

# EXAMPLE
library(microbenchmark)
x = seq(0,200,length.out=1000)
pars.exp = c(u=100, w=10, h=100) # expected parameter values
y_nonoise = gauspeak(x=x, u=pars.exp["u"], w=pars.exp["w"], h=pars.exp["h"])
set.seed(123)
y = rpois(length(y_nonoise),y_nonoise)
plot(x,y)
microbenchmark(fit <- logconcfit(x=x, y=y, control=list(yhatstart=NULL, method="conspline", nknots=20)), times=100) # 4.3ms
lines(x, fit$fitted, col="red") 
cor(y_nonoise,fit$fitted)^2 # 99.969 -> best
logconcfit(x=x, y=y, xout=0:100, control=list(yhatstart=NULL, method="conspline", nknots=20))
microbenchmark(fit <- logconcfit(x=x, y=y, control=list(yhatstart=NULL, method="cobs", tau=0.5, nknots=20)), times=100) # 8ms
lines(x, fit$fitted, col="green3") # OK
cor(y_nonoise,fit$fitted)^2 # 99.932 -> 2nd best
microbenchmark(fit <- logconcfit(x=x, y=y, control=list(yhatstart=NULL, method="conreg", nknots=20)), times=100) # 17ms
lines(x, fit$fitted, col="blue") # this one seems too pointy, but maybe good for triangular peaks??
cor(y_nonoise,fit$fitted)^2 # 99.947
microbenchmark(fit <- logconcfit(x=x, y=y, control=list(yhatstart=NULL, method="splconreg")), times=100) # 20ms
lines(x, fit$fitted, col="orange") # OK
cor(y_nonoise,fit$fitted)^2 # 99.912


# FUNCTION TO CARRY OUT UNIVARIATE LOG-CONCAVE SPLINE FITS ACCORDING TO METHOD OF EILERS (2005) ####
# original code by Paul Eilers (1998)
# see paper by Eilers (2005), with modifications by Tom Wenseleers (2018)
# currently not used since function has some numerical stability problems if domain over which fit is done is chosen too large
# plus it's quite slow & there are too many parameters to tune
logconcfit.eilers = function (x, y, 
                              lambda = 1, # smoothness penalty regularization penalty, set to 0 to disable
                              ord = 4, # order of differences to use in smoothness penalty
                              nknots = round(length(x[x!=0])/3),  # number of spline knots
                              degree = 3, # B spline degree
                              lambda_ridge = 1E-8, # small ridge regularization parameter, smaller value results in curve better following the peaks
                              gamma = 1E7, # regularization parameter to enforce log-concavity, higher values enforce stronger log-concavity, e.g. 1E7, smaller, e.g. 5, less
                              delta = -1E-9, # added to help solve numerical issues
                              yshift = 0.0001, # added to help solve numerical issues
                              yinit = NULL, # vector with initialization for fit
                              w = 0.001, # peak width of Gaussian peak shape used for initialization if yinit=NULL
                              tol = 0.001, # convergence criterion in terms of proportional change in MSE
                              niter = 50, # number of desired iterations (set to large value if you specify tol)
                              verbose = TRUE # print fit progress? 
                              ) { # plot fit?
  
  if (!is.null(yinit)) yinit = matrix(yinit, ncol=1)
  
  # rescale x and y to go between 0 and 1 
  # yshift = 0
  x = as.vector(x)
  y = as.vector(y)
  xmin = min(x)
  xmax = max(x)
  x = (x - xmin) / (xmax - xmin)
  ymax = max(y)
  y = y / ymax + yshift
  x = matrix(x,ncol=1)
  y = matrix(y,ncol=1)
  m = nrow(y)
  # plot(x,y)
  
  # set up penalty matrix P for smoothness
  B = fastbs(x, nknots = nknots, degree = degree)
  # matplot(B, type="l", lty=1)
  nb = ncol(B)
  D = diff(diag(1, nb), differences = ord)
  P = lambda * t(D) %*% D # P spline bases
  # matplot(P, type="l", lty=1)
  
  # set up ridge penalty regularization matrix R
  R = lambda_ridge * diag(1, nb)
  # matplot(R, type="l", lty=1)
  
  # set up  penalty matrix for unimodality (log-concavity) D2
  D2 = diff(diag(1, nb), differences = 2)
  # matplot(D2, type="l", lty=1)
  
  # do log-concave spline fit
  require(MASS) # for ginv function
  w2 = matrix(0, nb - 2, 1)
  it = 0
  # cm = log(mean(y)) * matrix(1, nb, 1) # check for overflows?
  if (is.null(yinit)) yinit = matrix(gausfit(x=x, y=(y-yshift)*ymax, control=list(fast=TRUE))$fitted/ymax,
                                     ncol=1) # we initialize with a gaussian peak fit if no initialization vector yinit is given
  # lines(x,yinit,col="red")
  cm = ginv(B, tol = 0) %*% yinit
  C = c()
  oldMSE = 1e10
  newMSE = 1e9
  propchangeMSE = 1
  
  while ((propchangeMSE > tol) & (it < niter)) {
    eta = B %*% cm
    eta[nrow(eta), 1] = eta[nrow(eta)-1, 1] # ADDED AS A FIX FOR LAST OUTLIER POINT
    eta[eta<(-570)] = -570 # added as fix for underflows
    eta[eta>10] = 10 # added as fix for overflows
    # plot(eta)
    mu = exp(eta) # this one can easily give overflows
    # plot(mu)
    Mu = mu %*% matrix(1, 1, nb)
    G = B * Mu
    c0 = cm
    Q = t(D2) %*% diag(as.vector(w2)) %*% D2
    A = (t(G) %*% G + P + gamma * Q + R)
    r = t(G) %*% (y - mu + mu * eta)
    cm = ginv(A, tol = 0) %*% r
    C = cbind(C, cm)
    # delta = -yshift # y shift 1E-8
    w2 = (D2 %*% cm > delta) * 1
    it = it + 1
    # print(norm(cm - c0, "2"))
    oldMSE = newMSE
    newMSE = sum((y - mu) ^ 2) / m # mean squared error
    # if (newMSE>oldMSE) { newMSE = oldMSE
    #                      cm=c0 }
    if (verbose) print(paste0("MSE iteration ", as.character(it), ": ", newMSE))
    propchangeMSE = abs((newMSE - oldMSE) / oldMSE)
    # plot(x, y, main=paste0("Log-concave spline fit (iteration ",as.character(it),")") )
    # lines(x, mu, col="red") # log-concave fit
  }
  
  # scale x and y and mu back to original scale
  x = x * (xmax - xmin) + xmin
  y = (y - yshift) * ymax
  mu = (mu - yshift) * ymax
  
  return(list(x = x, y = y, fitted = mu))
  
}

# EXAMPLE

library(microbenchmark)
x = seq(0,200,length.out=1000)
pars.exp = c(u=100, w=10, h=100) # expected parameter values
y_nonoise = gauspeak(x=x, u=pars.exp["u"], w=pars.exp["w"], h=pars.exp["h"])
y = rpois(length(y_nonoise),y_nonoise)
microbenchmark(fit <- logconcfit.eilers(x=x, y=y, verbose=FALSE), times=5) # ca 1s
plot(x, y, main = "Log-concave spline fit (red)")
lines(x, fit$fitted, col = "red") # log-concave fit



# FUNCTION TO FIT PEAK USING GAUSSIAN OR LOG-CONCAVE PEAK MODEL AND RETURN THE FIT WITH BEST AIC or BIC ####

fitpeak = function(x, y, xout=x, fitfuncs=list("gausfit", "logconcfit"), crit="BIC", 
                   control=list(start=NULL, fast=TRUE, method="cobs", tau=0.5, nknots=20),
                   envir = parent.frame() ) { 
  fits = Map(function(fitfunc) do.call(fitfunc, list(x=x, y=y, xout=xout, control=control), envir = envir), setNames(fitfuncs, fitfuncs))
  infcrit = sapply(fits, function (fit) fit[["IC"]][[crit]]) # vector of AICs or BICs of fits
  # print(infcrit)
  best = which.min(infcrit)
  bestfit = fits[[best]]
  return(bestfit)
}

# EXAMPLE
library(microbenchmark)
x = seq(0,200,length.out=1000)
pars.exp = c(u=100, w=10, h=100) # expected parameter values
y_nonoise = gauspeak(x=x, u=pars.exp["u"], w=pars.exp["w"], h=pars.exp["h"])
y = rpois(length(y_nonoise),y_nonoise)
microbenchmark(bestfit <- fitpeak(x=x, y=y, fitfuncs=list("gausfit", "logconcfit"),
                                  crit="BIC", 
                                  control=list(start=NULL, fast=TRUE, method="cobs", tau=0.5, nknots=20) ), times=100) # 8.8ms
plot(x, y, main = paste0("Best fit = ",bestfit$fitfunc," (red)"))
lines(x, bestfit$fitted, col = "red") # best fit





# 2. FUNCTION TO FIT SUM OF GAUSSIAN OR LOG-CONCAVE FUNCTIONS USING BACKFITTING APPROACH ####

# function to decompose signal as sum of peaks located at peak maxima with given shape (gaussian, log-concave spline,...) 
# using backfitting
fitpeakmaxima = function(x, y, 
                         peakM = getpeaks(x, y, sortby="height", plot=TRUE),
                         fitfuncs = list("gausfit", "logconcfit"),
                         crit="BIC", 
                         control=list(start=NULL, fast=FALSE, method="conspline", tau=0.5, nknots=20),
                         Cinit=NULL, maxiter=10, plot=TRUE, verbose=TRUE) {
  
  # TO DO: chage fitgaus to fitpeak
  # maybe after a couple of iterations sort peakM matric by elution order?
  # after couple of iterations no longer use all fitfuncs, but only type inferred in previous iteration
  # implement active set algo & check for convergence per peak, update peakM at end with fitted peak params or width inferred from fwhm
  require(nnls)
  npeaks = nrow(peakM)
  nobs = length(x)
  npars = vector(mode="numeric", length=npeaks) # nr of parameters for each peak that is fitted
  fits = vector("list", length=npeaks) # list with fits per peak

  if (is.null(Cinit)) { C = matrix(0, ncol=npeaks, nrow=nobs) } else { C = Cinit }

for (iter in 1:maxiter) { # backfitting iterations
  Cold = C
  for (peaknr in 1:npeaks) { # TO DO: only iterate over peaks that haven't reached stationary point yet, ie use active set
    skippeak = FALSE
    resid = y-rowSums(C[,-peaknr]) # partial residuals for peak peaknr
    # if (iter==-1) resid = resid-min(resid) # to avoid negative residuals in 1st iteration
    # plot(x,resid,type="l") 
    resid = pmax(resid,0) # cutt off residuals at 0
    
    # During the first iteration we first try to fit on local range 
    # corresponding to local peak & nearby local minima.
    # If that extrapolates well outside this range we use fit on extended range (e.g. for triangular peaks). 
    # In later iterations we always fit on a local range, determined after first determining a suitable range after the first iteration and updating peaks accordingly.
    
    # fit on local range
    x_start = peakM[peaknr,"start"]
    x_stop = peakM[peaknr,"stop"]
    x_range = x_start:x_stop # subsetted range corresponding to given peaknr
    x_local = x[x_range] # maybe extend range somwhat?? (by factor of 2?)
    y_local = resid[x_range]
    l = x_stop-x_start+1
    if ((iter==-1)&(is.null(Cinit))) {
      # we subtract a linear baseline connecting two endpoints of each peak during the 1st iteration
      # currently disabled (results in very narrow peaks after 1st iteration)
      baseline = y_local[1] + (y_local[l]-y_local[1])/(x_local[l]-x_local[1]) * (x_local-x_local[1]) 
      y_local = pmax(y_local-baseline, 0)
    } else { baseline = rep(0,length(y_local)) }
    if (sum(y_local)==0) skippeak=TRUE 
    if (!skippeak ) {
      # if (iter==1) {yinit=NULL} else {yinit=C[x_range,peaknr]} # use warm starts if possible
      fit_local <- fitpeak(x=x_local, y=y_local, xout=x, fitfuncs=fitfuncs, crit=crit, control=control)
      fits[[peaknr]] <- fit_local

      # if predictions of fit_local extrapolate well outside of fitted domain we carry out fit on extended range (e.g. for triangular peaks)
      if (iter==1) {
        y_extrapol = rep(0,nobs)
        y_fitted = fit_local$fitted+baseline
        xmax = which.max(y_fitted)
        y_extrapol[x_range] = y_fitted
        if (x_range[1]>1) y_extrapol[1:(x_range[1]-1)] = y_fitted[1] + (y_fitted[xmax]-y_fitted[1])/(x_local[xmax]-x_local[1]) * ((1:(x_range[1]-1))-x_local[1]) # linear extrapolation on left-hand side
        if (x_range[l]<nobs) y_extrapol[x_range[l]:nobs] = y_fitted[xmax] + (y_fitted[l]-y_fitted[xmax])/(x_local[l]-x_local[xmax]) * ((x_range[l]:nobs)-x_local[xmax]) # linear extrapolation on right-hand side
        y_extrapol = pmax(y_extrapol,0)
        # plot(x,resid)
        # lines(x,y_extrapol, col="red",lwd=2)
        # cor(y_extrapol,resid) # 0.94
        y_extrapol[is.na(y_extrapol)] = 0 
      }
      
      if ((iter==2)&(cor(y_extrapol,resid)>0.85)) { 
        # fit on full range (or modify this to use an extended range, but not the full range??)
        fit_full <- fitpeak(x=x, y=resid, fitfuncs=fitfuncs, crit=crit, control=control)
        fits[[peaknr]] <- fit_full
        y_fitted = fit_full$fitted
        npars[peaknr] = fit_full$npars
        } else { 
          y_fitted = fit_local$preds
          npars[peaknr] = fit_local$npars
          }
      
      C[,peaknr] = y_fitted
      if (iter==2) {
         # update start & stop of peak in peaks list
         peakM[peaknr,"start"] = max(head(rev(which(y_fitted[1:which.max(y_fitted)]<1E-10)),1)-10,1)
         peakM[peaknr,"stop"] = min(which.max(y_fitted)+head(which(y_fitted[which.max(y_fitted):nobs]<1E-10),1)+10,nobs)
       }  
      # matplot(C, type="l") 
    }
    # peaknr = peaknr+1
  }
  # after each cycle of updates we optimally rescale the obtained elution profiles using nnls (enhances convergence speed)
  nnlsfit = nnls(C, y)
  C = C %*% diag(coef(nnlsfit))
  Cs = cbind(Cold, C)
  # cossim = sapply(1:ncol(C), function(col) tcrossprod(C[,col,drop=F]/sqrt(rowSums(C[,col,drop=F])), Cold[,col,drop=F]/
  #                                                      sqrt(rowSums(Cold[,col,drop=F]))) )
  # print(cossim)
  
  if (plot|iter==maxiter) { # (iter %% 5)==1
    plot(x,y, col="blue", lwd=1, ylab="Signal (blue) + fit (red)", main=paste0("Explained variance after ", as.character(iter), " backfitting iteration(s) (R2=", as.character(round(100*cor(sqrt(rowSums(C)), sqrt(y))^2,2)), "%)"))
    lines(x,rowSums(C), type="l", col="red", lwd=5, lty=1)
    matplot(C, type="l", add=TRUE, lwd=2, lty=1) 
  }
  if (verbose&(ncol(C)!=0)&sum(rowSums(C))!=0) message(cor(sqrt(rowSums(C)), sqrt(y))^2) }
  if (plot) plot(x,y-rowSums(C), type="l", ylab="Residuals", main=paste0("Residuals after ", as.character(iter), " backfitting iterations"))

  # update peak start, stop, range, height & width in peaks list
peakM[,"start"] = apply(C, 2, function (y) max(head(rev(which(y[1:which.max(y)]<1E-10)),1),1))
peakM[,"stop"] = apply(C, 2, function (y) min(which.max(y)+head(which(y[which.max(y):nobs]<1E-10),1),nobs))
peakM[,"height"] = sapply(1:npeaks, function (peaknr) fits[[peaknr]]$fitted.pars[["h"]])
peakM[,"width"] = sapply(1:npeaks, function (peaknr) fits[[peaknr]]$fitted.pars[["w"]])
peakM[,"mode"] = sapply(1:npeaks, function (peaknr) fits[[peaknr]]$fitted.pars[["u"]])
peakM[,"range"] = peakM[,"stop"] - peakM[,"start"]

return(list(x=x, y=y, peakM=peakM, C=C, fits=fits, IC=IC(y=y, yhat=rowSums(C), npars = sum(npars)), maxiter=maxiter))
}

# EXAMPLE
x = 1:1000
gpeaks = sumgauspeaks(x=x, u=c(200, 320, 600, 500), w=c(50,50,50,200), h=c(1000, 500, 300, 1000), plot=TRUE)
y = gpeaks$y
microbenchmark(peakM <- getpeaks(x, y, sortby="height", plot=TRUE), times=1) # 3.7 ms with plot=FALSE, 43 ms with plot=TRUE
microbenchmark(fit1 <- fitpeakmaxima(x=x, y=y, 
                       peakM = peakM,
                       fitfuncs = list("gausfit","logconcfit"), # specify just gausfit if you would like to force Gaussian peaks 
                       crit="BIC",
                       control=list(start=NULL, fast=FALSE, method="conspline", tau=0.5, nknots=20),
                       Cinit=NULL, maxiter=10, plot=T), times=1) # 12ms with plot=FALSE

matplot(gpeaks$Y, add=TRUE, type="l", lty=1, lwd=0.1)
fit1$IC$AIC # 2991.96
fit1$IC$BIC # 3036.13
cor(sqrt(y),sqrt(rowSums(fit1$C)))^2 # 98.9


### now do forward stepwise addition of peaks based on residual (1 by one) as long as this improves AIC or BIC
# code from here on is not finished yet...
crit = "BIC"

C = fit1$C
matplot(C, type="l")

resid = y-rowSums(C)
# resid = resid-min(resid) # to avoid negative residuals
resid = pmax(resid,0) # cutt off residuals at 0
plot(x, resid, type="l") 
# now identify biggest peak in residual

microbenchmark(peaks_new <- getpeaks(x=x, y=resid, sortby="height", plot=TRUE), times=5) # 45ms
peakM2 = rbind(peakM, peaks_new[1,]) # add biggest peak to existing peak matrix
C2 = cbind(C,0)
# plot(x,y)
# lines(x,rowSums(C2),col="blue")
# now redo some backfitting iterations and check if AIC or BIC got better
fit1$IC$AIC # 2991.96
fit1$IC$BIC # 3036.13

microbenchmark(fit2 <- fitpeakmaxima(x=x, y=y, 
                                     peakM = peakM2,
                                     fitfuncs = list("gausfit"), # ,"logconcfit"
                                     crit="BIC",
                                     control=list(start=NULL, fast=FALSE, method="conspline", tau=0.5, nknots=20),
                                     Cinit=C2, maxiter=10, plot=TRUE, verbose=FALSE), times=1) # 456ms with plot&verbose=FALSE
matplot(gpeaks$Y, add=TRUE, type="l", lty=2, lwd=0.1)
fit2$IC$AIC # 1466.674
fit2$IC$BIC # 1525.567
peakM3 <- fit2$peakM

# comparison with nls (Levenberg-Marquardt) 
library(minpack.lm)
microbenchmark(nlsfit2 <- nlsLM(y ~ gauspeak(x, u1, w1, h1) + gauspeak(x, u2, w2, h2) 
                                  + gauspeak(x, u3, w3, h3) + gauspeak(x, u4, w4, h4),
                 data = data.frame(x=x, y=y),
                 control = nls.lm.control(maxiter = 100),
                 start = list(u1 = peakM3[1, "mode"],
                              u2 = peakM3[2, "mode"],
                              u3 = peakM3[3, "mode"],
                              u4 = peakM3[4, "mode"],
                              w1 = 75,
                              w2 = 75,
                              w3 = 75,
                              w4 = 75,
                              h1 = peakM3[1, "height"],
                              h2 = peakM3[2, "height"],
                              h3 = peakM3[3, "height"],
                              h4 = peakM3[4, "height"]
                 )), times=1) # 334ms
# solution now looks perfect :
Yhat <- cbind(gauspeak(x, coef(nlsfit2)["u1.mode"], coef(nlsfit2)["w1"], coef(nlsfit2)["h1.height"]),
              gauspeak(x, coef(nlsfit2)["u2.mode"], coef(nlsfit2)["w2"], coef(nlsfit2)["h2.height"]),
              gauspeak(x, coef(nlsfit2)["u3.mode"], coef(nlsfit2)["w3"], coef(nlsfit2)["h3.height"]),
              gauspeak(x, coef(nlsfit2)["u4.mode"], coef(nlsfit2)["w4"], coef(nlsfit2)["h4.height"]))
yhat <- rowSums(Yhat)
plot(x,y, col="blue")
lines(x,yhat,col="red",lwd=2)
matplot(Yhat, add=TRUE, type="l", lty=1, lwd=4)
matplot(gpeaks$Y, add=TRUE, type="l", lty=2, lwd=0.1)
cor(sqrt(y),sqrt(yhat))^2 # 99.75702
IC(y=y, yhat=yhat, npars=3*4) # AIC=1462.172, BIC=1521.065


# example of fit on NNDSVD results
library(R.matlab)
y=as.vector(readMat("exampleNNSVDH.mat")$H) # example factor from a NNDSVD
x=1:length(y)
plot(x,y,type="l")
write.csv(data.frame(x=x,y=y),file="difficult.peaks.csv")

microbenchmark(peakM <- getpeaks(x, y, sortby="height", plot=TRUE), times=1) # TO DO: maybe also cluster peaks that are too close together??
microbenchmark(fit <- fitpeakmaxima(x=x, y=y, # still some bugs here - check out where they come from...
                                     peakM = peakM,
                                     fitfuncs = list("gausfit","logconcfit"), # specify just gausfit if you would like to force Gaussian peaks 
                                     crit="BIC",
                                     control=list(start=NULL, fast=TRUE, method="cobs", tau=0.5, nknots=10),
                                     Cinit=NULL, maxiter=10, plot=TRUE), times=1) 

fit$IC$AIC # 2991.96
fit$IC$BIC # 3036.13
cor(sqrt(y),sqrt(rowSums(fit1$C)))^2 # 98.9

