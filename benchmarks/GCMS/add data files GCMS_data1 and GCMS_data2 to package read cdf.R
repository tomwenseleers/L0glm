# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("xcms")
library(xcms)
readcdf = function(file, mz.ignore = NULL, profmethod = "bin", profparam = list(offset = -0.5)) {
  f=xcmsRaw(file, profstep=1, profmethod=, profparam = profparam) # set profstep to 0 to read in the raw data
  mzrange=f@mzrange
  mzvals=mzrange[[1]]:mzrange[[2]]
  scantimes=f@scantime/60 # scantimes in mins
  M=t(f@env$profile) # data matrix
  M[M<=0] = 0
  colnames(M)=mzvals
  rownames(M)=scantimes
  if(!is.null(mz.ignore)){
    M = M[,-which(colnames(M) %in% mz.ignore)] # remove these ion traces
  }
  M[is.na(M)]=0
  M[M<0]=0
  return(M)
}



# function to plot matrix as image, also taking sparse matrix as input
plotmat <- function(M, 
                    row_subsamp=round(nrow(M)/min(nrow(M),2000)), 
                    col_subsamp=round(ncol(M)/min(ncol(M),2000)),
                    gamma=0.12, ...) { 
  cols = colorRampPalette(c("black", "black","blue","green","yellow","orange", "red"))(1000)
  M = M[seq(1,nrow(M),by=row_subsamp),
        seq(1,ncol(M),by=col_subsamp), drop=FALSE]
  M=as.matrix(M) # convert sparse to dense matrix
  M[is.na(M)]=0
  # M[M<0]=0
  if (!is.null(as.numeric(rownames(M))[[1]])) { xrange = range(as.numeric(rownames(M)))
  xvals = seq(min(xrange), max(xrange), length.out = nrow(M)) } else { xvals = 1:nrow(M) }
  if (!is.null(as.numeric(colnames(M))[[1]])) { yrange = range(as.numeric(colnames(M)))
  yvals = seq(min(yrange), max(yrange), length.out = ncol(M)) } else { yvals = 1:ncol(M) }
  image(x=xvals,
        y=yvals,
        z=M^gamma, 
        col = cols,
        useRaster=TRUE,
        ...)
}

## SAMPLE 1 ####
# sample 1: queen hydrocarbon profiles of the ant Cataglyphis iberica (small section of chromatogram), cf. Van Oystaeyen et al. Science 2014
# here all peaks have a well behaved near-Gaussian shape
samplefile1 = "C:/Users/bherr/Dropbox/gcms_deconvolution_project/christophe/chromatograms/cataglyphis_iberica/queens/CI23_Q_120828_01.CDF"
library(ncdf4)
cdf_info <- nc_open(samplefile1) # see https://pjbartlein.github.io/REarthSysSci/netCDF.html how to parse info
cdf_info
system.time(raw1 <- readcdf(file=samplefile1)) # 3.9s
plotmat(raw1, subsamp=10) # whole chromatogram
# small section containing at 
# RT 26.84 15-methylnonacosane (diagn ions 224+225)
# RT 26.91 13-Methylnonacosane (diagn ions 196+225)
# RT 26.94 11-Methylnonacosane (diagn ions 168+280)
# RT 26.95 Eicosamethylcyclodecasiloxane (diagn ions 73+147+221+281+355+429+503)
# RT 27.01 9-Methylnonacosane (diagn ions 140+308)
# RT 27.11 7-Methylnonacosane (diagn ions 112+336)
section <- as.numeric(rownames(raw1))>26.8&as.numeric(rownames(raw1))<27.2
raw1_sect <- raw1[section,,drop=F] # small section
plotmat(raw1_sect, subsamp=1)
rownames(raw1_sect)
colnames(raw1_sect) # m/z 30 to 700
dim(raw1_sect) # 119 671
saveRDS(raw1_sect, file="../benchmarks/GCMS/GCMS_data1.rds")
saveRDS(raw1, file="../benchmarks/GCMS/GCMS_data1_full.rds")
# C9 to C40 linear alkane ladder from same sample batch
samplefile1_alk = "C:/Users/bherr/Dropbox/gcms_deconvolution_project/christophe/chromatograms/cataglyphis_iberica/alkstandards/STD_AL0001_120716.CDF"
system.time(raw1_alk <- readcdf(file=samplefile1_alk)) # 3.9s
saveRDS(raw1_alk, file="../benchmarks/GCMS/GCMS_data1_alk_full.rds")
plotmat(raw1_alk, subsamp=10) # whole chromatogram
samplefile1_alk_TIC = rowSums(raw1_alk) # total ion chromatogram
names(samplefile1_alk_TIC) = rownames(raw1_alk) # RTs
plot(x=as.numeric(names(samplefile1_alk_TIC)), y=samplefile1_alk_TIC, type="l", 
     xlab="Retention time (min)", ylab="Total ion intensity", main="C9 to C40 alkane ladder") 
GCMS_data1 <- raw1_sect
attr(GCMS_data1, "alkanes_TIC") <- samplefile1_alk_TIC # we add the alkane ladder TIC as an attribute to data(GCMS_data1)
save(GCMS_data1, file="./data/GCMS_data1.rda")


## SAMPLE 2 ####
# sample 2: more challenging example: SPME-HS of beer volatiles, from Scheurs et al. Nature Comm. 2024
# here peaks with a a well behaved near-Gaussian shape co-occur with a large
# acid peak with a triangular sharkfin shape (data below is a section of the whole chromatogram)
samplefile2 = "C:/Users/bherr/Dropbox/gcms_deconvolution_project/christophe/chromatograms/beers/20170124/54040014_a.cdf"
library(ncdf4)
cdf_info <- nc_open(samplefile2) # see https://pjbartlein.github.io/REarthSysSci/netCDF.html how to parse info
cdf_info
system.time(raw2 <- readcdf(file=samplefile2)) # 3.9s
plotmat(raw2, subsamp=10) # whole chromatogram
# raw2 <- readRDS("C:/Users/bherr/Documents/Github/L0glm/benchmarks/GCMS/GCMS_beersample.rds")
# colnames(raw2) <- as.character(33:550)
# rownames(raw2) <- seq(20, 30, length.out=nrow(raw2))
# dim(raw2) # 947 x 518 (time x mass)
section <- as.numeric(rownames(raw2))>19&as.numeric(rownames(raw2))<25
raw2_sect <- raw2[section,,drop=F] # small section
plotmat(raw2_sect, subsamp=1)
rownames(raw2_sect) # 19.00 to 25.00 mins
colnames(raw2_sect) # 33 to 550
dim(raw2_sect) # 1764  518
saveRDS(raw2_sect, file="../benchmarks/GCMS/GCMS_data2.rds")
saveRDS(raw2, file="../benchmarks/GCMS/GCMS_data2_full.rds")
# C7 to C22 linear alkane ladder from same sample batch
samplefile2_alk = "C:/Users/bherr/Dropbox/gcms_deconvolution_project/christophe/chromatograms/beers/20170124/alkanes.cdf"
system.time(raw2_alk <- readcdf(file=samplefile2_alk)) # 3.9s
saveRDS(raw2_alk, file="../benchmarks/GCMS/GCMS_data2_alk_full.rds")
raw2_alk <- raw2_alk[as.numeric(rownames(raw2_alk))>6,,drop=F] # we leave out the overloaded hexane/C6 in the beginning
plotmat(raw2_alk, subsamp=10) # whole chromatogram from RT of 6 mins, with C7 to C22
samplefile2_alk_TIC = rowSums(raw2_alk) # total ion chromatogram
names(samplefile2_alk_TIC) = rownames(raw2_alk) # RTs
plot(x=as.numeric(names(samplefile2_alk_TIC)), y=samplefile2_alk_TIC, type="l", 
     xlab="Retention time (min)", ylab="Total ion intensity", main="C7 to C22 alkane ladder") 
GCMS_data2 <- raw2_sect
attr(GCMS_data2, "alkanes_TIC") <- samplefile2_alk_TIC # we add the alkane ladder TIC as an attribute to data(GCMS_data2)
save(GCMS_data2, file="./data/GCMS_data2.rda")
