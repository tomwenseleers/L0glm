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

npeaks <- 10

GCMSPrimedir = "Z:/Chris/GCMSPrime/"
load(file=paste0(GCMSPrimedir,"adamsDB.rda"))
mzout=41:500

ord <- order(sapply(adamsDB, "[[", "RI"))
adamsDB <- adamsDB[ord]

spectralM=t(spectomat(adamsDB,mzout=mzout))
dim(spectralM) # 2204 rows (spectra), columns are m/z values
sum(rowSums(spectralM)==1) # 2202

start <- round(nrow(spectralM)/2)
nspecs <- 100
S <- 1E3 * spectralM[start+(1:nspecs),,drop=FALSE] # matrix with selected spectra in rows, npeaks x 460 m/z
S <- t(apply(S, 1, function(x) pmax(x - min(x[x!=0]), 0)))
S <- S + matrix(rpois(length(S), lambda = as.vector(S)), ncol = ncol(S), nrow = nrow(S))
S <- t(apply(S, 1, function(x) x/sqrt(sum(x^2))))

saveRDS(object = S, file = "./inst/Simulated spectrum.rds")
