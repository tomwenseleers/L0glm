# POSSIBLE DIMENSIONALITY REDUCTION TECHNIQUE FOR MULTICHANNEL CASE ####

library(L0glm)
dat=simulate_spike_train(n=10000, p=10000, k=100, family="poisson", multichannel=T)
Y=dat$y
weights=1/(Y+0.1) # approx Poisson 1/variance weights
Yw = Y*sqrt(weights) # to take into account approx 1/variance Poisson weights
L2n = sqrt(colSums(Yw^2)) # L2 norm normalize rows
Yw = sweep(Yw, 2, L2n, FUN="/")
Yw[is.na(Yw)] = 0
colnames(Yw) = NULL
plotmat(Yw, gamma=0.3)


## calculate matrix with weighted dot products (using 1/(Y+0.1) as approx 1/variance weights) ####
# system.time(weighted_dot_matrix <- weightedDotProduct(Y, 1/(Y + 0.1))) # 0.05s (this function is coded in Rcpp & uses multithreading via OpenMP)

## Cluster co-correlated outcome variables, here using model-based clustering with automatic detection of optimal nr of clusters (BEST) ####
library(mclust)
system.time(clustering_result <- Mclust(t(Yw[seq(1, nrow(Yw), by=30),]),  
                                        G=1:round(ncol(Yw)/2), verbose=F)) # 7s

# Optimal number of clusters based on BIC criterion
optimal_clusters <- clustering_result$G
optimal_clusters
# 9

# Clusters
clusters <- clustering_result$classification

# Sets of co-correlated outcome variables
lapply(1:optimal_clusters, function (group) which(clusters %in% group))


## Sum signal of these co-correlated outcome variables #
# note: the sum of Poisson distributed variables would still be Poisson

# Number of clusters
num_clusters <- max(clustering_result$classification)

# Initialize an aggregate matrix
aggregate_matrix <- matrix(0, nrow(Y), num_clusters)

# Sum columns by cluster
for (cluster in 1:num_clusters) {
  # Find columns that belong to the current cluster
  cluster_columns <- which(clustering_result$classification == cluster)
  # print(cluster_columns)
  
  # Sum these columns and store the result in the aggregate matrix
  aggregate_matrix[, cluster] <- rowSums(Y[, cluster_columns, drop=F])
}

aggregate_matrix = aggregate_matrix[,order(colSums(aggregate_matrix), decreasing=T)] # order by amount of signal

plotmat(aggregate_matrix, ylab="Mass (m/z) cluster")
plot(rowSums(aggregate_matrix)^0.5, type="l")
for (m in 1:ncol(aggregate_matrix)) {
  plot(aggregate_matrix[,m]^0.5, type="l")
}




# TEST WITH ACTUAL GC/MS SAMPLE ####

library(xcms)
readcdf = function(file, mz.ignore = NULL) {
  f=xcmsRaw(file,profstep=1,profmethod="intlin") # set profstep to 0 to read in the raw data
  mzrange=f@mzrange
  mzvals=mzrange[[1]]:mzrange[[2]]
  scantimes=f@scantime/60 # we return scantimes in mins
  M=t(f@env$profile) # data matrix
  M[M<=0] = 0
  colnames(M)=mzvals
  rownames(M)=scantimes
  if(!is.null(mz.ignore)){
    M = M[,-which(colnames(M) %in% mz.ignore)] # Remove these ion traces
  }
  M[is.na(M)]=0
  M[M<0]=0
  return(M)
}

# function to plot matrix as image, also taking sparse matrix as input
plotmat <- function(M,subsamp=10,gamma=0.12,xlab="Retention time (mins)",ylab="Mass (m/z)", ...) { 
  colpalette=colorRampPalette(c("black", "black","blue","green","yellow","orange", "red")) 
  M=M[seq(1,nrow(M),subsamp),]
  M=as.matrix(M) # convert sparse to dense matrix
  M[is.na(M)]=0
  M[M<0]=0
  if (length(rownames(M))!=0) { xvals = as.numeric(rownames(M)) } else { xvals = 1:nrow(M) }
  xrange = range(xvals)
  image(x=1:nrow(M),
        y=1:ncol(M),
        z=M^gamma, 
        col = colpalette(1000),
        useRaster=T,
        xlab=xlab,ylab=ylab,
        axes=FALSE, 
        ...)    # xaxt="n" # ylim=c(180,350) to zoom in on mass spec
  xlabls=as.character(round(xvals))
  xsel=((as.numeric(xlabls) %% 10)==0)&(c(TRUE,diff(as.numeric(xlabls))!=0))
  ylabls=as.character(round(as.numeric(colnames(M))))
  ysel=((as.numeric(ylabls) %% 100)==0)&(c(TRUE,diff(as.numeric(ylabls))!=0))
  box()
  axis(1, at = (1:nrow(M))[xsel], labels=xlabls[xsel], tick=TRUE, tcl=-0.4)
  axis(2, at = (1:ncol(M))[ysel], labels=ylabls[ysel], tick=TRUE, las=1, tcl=-0.4)
  # plot(x=xvals, y=rowSums(M)+1E6, type="l",
  #     xlab="Retention time (mins)", ylab="Relative intensity", xaxs="i", log="y")
}


samplefile = "C:/GCMS/lucanus_cervus_vliegend_hert/4_lucanus_cervus_vliegend_hert_headspace_11feb2021/F_thorax/20_F3_thorax.cdf"
raw = readcdf(file=samplefile)
plotmat(raw)

Y=raw
weights=1/(raw+0.1) # approx Poisson 1/variance weights
Yw = Y*sqrt(weights) # to take into account approx 1/variance Poisson weights
L2n = sqrt(colSums(Yw^2)) # L2 norm normalize rows
Yw = sweep(Yw, 2, L2n, FUN="/")
Yw[is.na(Yw)] = 0
plotmat(Yw, gamma=0.3)

library(L0glm)
## calculate matrix with weighted dot products (using 1/(Y+0.1) as approx 1/variance weights) ####
# system.time(weighted_dot_matrix <- weightedDotProduct(Y, 1/(Y + 0.1))) # 0.75s (this function is coded in Rcpp & uses multithreading via OpenMP)

## Cluster co-correlated outcome variables, here using model-based clustering with automatic detection of optimal nr of clusters (BEST) ####
library(mclust)
system.time(clustering_result <- Mclust(t(Yw[seq(1, nrow(Yw), by=30),]), 
                                        G=1:round(ncol(Y)/20), verbose=F)) # 224s

# Optimal number of clusters based on BIC criterion
optimal_clusters <- clustering_result$G
optimal_clusters
# 4

# Clusters
clusters <- clustering_result$classification

# Sets of co-correlated outcome variables
clustersets <- lapply(1:optimal_clusters, function (group) which(clusters %in% group))


## Sum signal of these co-correlated outcome variables #
# note: the sum of Poisson distributed variables would still be Poisson

# Number of clusters
num_clusters <- max(clustering_result$classification)

# Initialize an aggregate matrix
aggregate_matrix <- matrix(0, nrow(Y), num_clusters)

# Sum columns by cluster
for (cluster in 1:num_clusters) {
  # Find columns that belong to the current cluster
  cluster_columns <- which(clustering_result$classification == cluster)
  
  # Sum these columns and store the result in the aggregate matrix
  aggregate_matrix[, cluster] <- rowSums(Y[, cluster_columns])
}

ord = order(colSums(aggregate_matrix), decreasing=T) # order by amount of signal
clustersets = clustersets[ord]
aggregate_matrix = aggregate_matrix[,ord]

plotmat(aggregate_matrix, ylab="Mass (m/z) cluster")
plot(rowSums(aggregate_matrix)^0.5, type="l")
for (m in 1:ncol(aggregate_matrix)) {
plot(aggregate_matrix[,m]^0.5, type="l")
}


# most of the signal seems to get aggregated in the first cluster ####
# maybe this cluster could still benefit from being subclustered? ####
plot(aggregate_matrix[,1]^0.5, type="l")

raw2 = raw[,clustersets[[1]]]
Y2 = raw2
weights2 = 1/(Y2+0.1) # approx Poisson 1/variance weights

Yw2 = Y2*sqrt(weights2)
L2n2 = sqrt(colSums(Yw2^2)) # L2 norm normalize rows
Yw2 = sweep(Yw2, 2, L2n2, FUN="/")
Yw2[is.na(Yw2)] = 0
plotmat(Yw2, gamma=0.3)

library(L0glm)
## calculate matrix with weighted dot products (using 1/(Y+0.1) as approx 1/variance weights) ####
# system.time(weighted_dot_matrix <- weightedDotProduct(Y, 1/(Y + 0.1))) # 0.03s (this function is coded in Rcpp & uses multithreading via OpenMP)

## Cluster co-correlated outcome variables, here using model-based clustering with automatic detection of optimal nr of clusters (BEST) ####
library(mclust)
system.time(clustering_result2 <- Mclust(t(Yw2[seq(1, nrow(Yw2), by=30),]),  
                                         G=1:round(ncol(Yw2)/20), verbose=F)) # 0.8s

# Optimal number of clusters based on BIC criterion
optimal_clusters2 <- clustering_result2$G
optimal_clusters2
# 9

# Clusters
clusters2 <- clustering_result2$classification

# Sets of co-correlated outcome variables
clustersets2 <- lapply(1:optimal_clusters2, function (group) which(clusters2 %in% group))


## Sum signal of these co-correlated outcome variables #
# note: the sum of Poisson distributed variables would still be Poisson

# Number of clusters
num_clusters2 <- max(clustering_result2$classification)

# Initialize an aggregate matrix
aggregate_matrix2 <- matrix(0, nrow(Y2), num_clusters2)

# Sum columns by cluster
for (cluster in 1:num_clusters2) {
  # Find columns that belong to the current cluster
  cluster_columns <- which(clustering_result2$classification == cluster)
  # print(cluster_columns)
  
  # Sum these columns and store the result in the aggregate matrix
  aggregate_matrix2[, cluster] <- rowSums(Y[, cluster_columns, drop=F])
}

ord2 = order(colSums(aggregate_matrix2), decreasing=T)
clustersets2 = clustersets2[ord2]
aggregate_matrix2 = aggregate_matrix2[,ord2] # order by amount of signal
aggregate_matrix2b = cbind(aggregate_matrix[,-1], aggregate_matrix2)
ord2b = order(colSums(aggregate_matrix2b), decreasing=T)
aggregate_matrix2b = aggregate_matrix2b[,ord2b] # order by amount of signal

plotmat(aggregate_matrix2b, ylab="Mass (m/z) cluster")
plot(rowSums(aggregate_matrix2b), type="l", log="y")
for (m in 1:ncol(aggregate_matrix2b)) {
  plot(aggregate_matrix2b[,m]^0.5, type="l", log="y")
}
