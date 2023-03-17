# Simulate data
# GC-MS data is simulated as a linear combination of spectrum profiles S and
# elution profiles C.

n <- p <- 100 # take as many covariate p as the number of observations n
x <- 1:n
npeaks <- 10 # nr of peaks = nr of spectra
set.seed(1)
peakhrange.log <- c(2, 5) # peak height range
w1 = 1.25 # simulated peak widths bigaussian peak shape
w2 = 2.5

# Function to generate an asymmetrical Gaussian peak
bigauspeak = function(x, u, w1, w2, h=1) {
  out <- rep(0, length(x))
  out[x<=u] <- h*exp(-(x[x<=u]-u)^2/(2*w1^2))
  out[x>u] <- h*exp(-(x[x>u]-u)^2/(2*w2^2))
  return(out)
}
# Banded matrix with theoretical peak shape function used, this will serve as
# the covariate matrix from which we sample the true set of variables for selection
# TODO you can still add background covariates in this one as well (e.g. B spline basis)
X <- sapply(1:p, function (u) bigauspeak(x, u = u, w1 = w1, w2 = w2, h = 1))
# Generate matrix C containing the true set of covariates
u <- sample(x, npeaks, replace=FALSE) # simulated peak locations
# h <- 10^runif(npeaks, min = min(peakhrange), max = max(peakhrange)) # simulated peak heights (in total ion current)
h <- 10^seq(2, 4, length.out = npeaks)
C <- sweep(X[,u], 2, h, "*")
par(mfrow = c(2,1))
matplot(C+1, log = "y", type = "l")
image(C^0.2, col = topo.colors(255))

# Generate the spectrum profiles
# Load real-life profiles
file <- "../benchmarks/multichannel/simulated_spectra.rds"
simulated_spectra <- readRDS(file)
m <- ncol(simulated_spectra) # 460
# Randomly select npeaks spectra
S <- simulated_spectra[sample(1:nrow(simulated_spectra), npeaks),]

# Simulate response variable
Y0 <- C %*% S # noiseless simulated signal = linear convolution of spike train with peak shape function
Y <- apply(Y0 + 1, 2, function (col) {
  out <- rpois(n, col) # simulated signal with random poisson noise on it
  out[is.na(out)] <- rnorm(n = sum(is.na(out)), mean = col[is.na(out)],
                           sd = sqrt(col[is.na(out)])) # to deal with overflows
  return(out)
})

# Plot the data
const <- 1
par(mfrow = c(1,1))
matplot(C+const, type="l", log="y", lwd=3, lty=1, col="red")
matlines(Y+const, lty = 1, col = rgb(0,0,0,0.2))
abline(v=u, col="blue", lwd=log10(h)/3)

# Case I: Fit a multichannel glm (no penalty)
system.time(fit1 <- group.pen.nnglm(X = X, Y = Y, family = poisson(identity),
                                   # Weights = matrix(1, nrow = nrow(Y), ncol = ncol(Y)),
                                   intercept = TRUE,
                                   group.pen.fun = function(mat) sqrt(rowSums(mat^2)), # default = L2 norm of the estimated coefficients per component
                                   lambdas = 0,
                                   maxit = 10))
par(mfrow = c(3,1))
image(Y0^0.2, col = topo.colors(255), main = "Ground truth")
image(Y^0.2, col = topo.colors(255), main = "Input data")
image(fit1$fitted.values^0.2, col = topo.colors(255), main = "Fitted data")

# Case II: Fit a multichannel L0 penalized glm
system.time(fit2 <- group.pen.nnglm(X = X, Y = Y, family = poisson(identity),
                                    # Weights = matrix(1, nrow = nrow(Y), ncol = ncol(Y)),
                                   intercept = TRUE,
                                   lambdas = rep(1E-2, m),
                                   # lambdas = 1E-5*(1+apply(Y, 2, max))^0.001, # lambda should be function of scaling of each mass channel, we specify small lambda because we initialise with nnls solution

                                   group.pen.fun = function(mat){
                                     out <- sqrt(rowSums(mat^2)) # row L2 norms are used for adaptive ridge penalty
                                     out <- out/max(out) # L infinity norm normalized L2 row norms to make them range between 0 and 1 - @CHRIS - maybe figure out why this normalisation works best - I also tried e.g. m*rownorms/sum(rownorms) but that didn't work so well
                                     return(out)
                                   },
                                   maxit = 30))

par(mfrow = c(3,1))
image(y = 0:m, x = 0:n, Y0^0.2, col = topo.colors(255), main = "Ground truth")
abline(v=u, col="green")
image(y = 0:m, x = 0:n, Y^0.2, col = topo.colors(255), main = "Input data")
abline(v=u, col="green")
image(y = 0:m, x = 0:n, fit2$fitted.values^0.2, col = topo.colors(255), main = "Fitted data")
abline(v=which(rowSums(fit2$coefficients)>1), col="green")


image(fit2$coefficients^0.2, col = topo.colors(255))



image(x=x, y=1:ncol(betaM_theor), betaM_theor^0.2,col = topo.colors(255), main="Ground truth", xaxt='n', xlab="", ylab="mass channel")
abline(v=u, col="green")
image(x=1:101, y=1:ncol(fit2$coefficients), fit2$coefficients^0.2,col = topo.colors(255), main="Per channel WNNLS sol (green=correct, red=false pos, orange=false neg)", xaxt='n', xlab="", ylab="mass channel")
# abline(v=u, col="red")
abline(v=intersect(which(rowSums(fit2$coefficients)>1E2),u+1), col="green")
abline(v=setdiff(which(rowSums(fit2$coefficients)>1E2),u+1), col="red")
abline(v=setdiff(u,which(rowSums(fit2$coefficients)>1E2)), col="orange")
image(x=x, y=1:ncol(betaM), betaM^0.2,col = topo.colors(255), main="Group L0-pen sol (green=correct, red=false pos, orange=false neg)", xaxt='n', xlab="", ylab="mass channel")
# abline(v=u, col="red")
abline(v=intersect(which(rowSums(betaM)>1E2),u), col="green")
abline(v=setdiff(which(rowSums(betaM)>1E2),u), col="red")










# load some spectra from adams library
library(matrixStats)
spectomat = function(lib,mzout,norm=TRUE) {
  minmz=min(mzout)
  nmz=length(mzout)
  outmat=do.call(cbind,lapply(1:length(lib),function(i) { ms=lib[[i]]$pspectrum
  out <- numeric(nmz)
  out[ms[,1]-minmz+1] <- ms[, 2]
  out } ))
  if (norm) outmat <- sweep(outmat, 2, colSums(outmat), FUN="/")
  # storage.mode(outmat) <- "integer"
  outmat
}
GCMSPrimedir="D:/Documents/GCMS/MyRPackages/GCMSPrime/"
#GCMSPrimedir = "Z:/Chris/GCMSPrime/"
load(file=paste0(GCMSPrimedir,"adamsDB.rda"))
mzout=41:500
spectralM=t(spectomat(adamsDB,mzout=mzout))
dim(spectralM) # 2204 rows (spectra), columns are m/z values
dev.off()
image(spectralM^0.2, col = topo.colors(255))
sum(rowSums(spectralM)==1) # 2202

# simulated data
require(Matrix)
n = 100 # nr of scanlines
x = 1:n # time
npeaks = 10 # nr of peaks = nr of spectra
set.seed(1234)
u = sample(x, npeaks, replace=FALSE) # simulated peak locations
peakhrange = c(100,1E5) # peak height range
h = 10^runif(npeaks, min=log10(min(peakhrange)), max=log10(max(peakhrange))) # simulated peak heights (in total ion current)
w1 = 1.25 # simulated peak widths bigaussian peak shape
w2 = 2.5
a_theor = rep(0, n) # simulated spike train
a_theor[u] = h


bM = do.call(cbind, lapply(1:n, function (u) bigauspeak(x, u=u, w1=w1, w2=w2, h=1) )) # banded matrix with theoretical peak shape function used, you can still add background covariates in this one as well (e.g. B spline basis)
startfrom = round(nrow(spectralM)/2)-npeaks # spectrum nr to start from in spectral library
S = spectralM[startfrom:(startfrom+npeaks-1),,drop=FALSE] # matrix with selected spectra in rows, npeaks x 460 m/z
image(S^0.2, col = topo.colors(255))
C = do.call(cbind, lapply(1:npeaks, function (i) bigauspeak(x, u=u[i], w1=w1, w2=w2, h=h[i])) ) # elution profiles
image(C^0.1, col = topo.colors(255))
Y_nonoise = C %*% S # noiseless simulated signal = linear convolution of spike train with peak shape function
image(x=x, y=mzout, Y_nonoise^0.2, col = topo.colors(255))
abline(v=u, col="red")
Y = apply(Y_nonoise, 2, function (col) { out = rpois(n, col) # simulated signal with random poisson noise on it
out[is.na(out)] = rnorm(n=sum(is.na(out)), mean=col[is.na(out)], sd=sqrt(col[is.na(out)])) # to deal with overflows
return(out) } )


const=10
matplot(Y+const, type="l", lty=1, log="y")
matplot((C/4+const), type="l", lwd=3, lty=1, col="red", add=TRUE)
abline(v=u, col="blue", lwd=log10(h)/3)
X=bM
m = ncol(X) # nr of covariates (e.g. peak shape models)
k = ncol(Y) # nr of linked Y channels (e.g. mass channels)
betaM_theor = matrix(0,nrow=m,ncol=k)
for (i in 1:npeaks) { betaM_theor[u[i],]=S[i,]*h[i] }
image(x=x, y=1:ncol(betaM_theor), betaM_theor^0.2, col = topo.colors(255)) # theoretical coef matrix as we would like to reconstruct it
abline(v=u,col="red")

