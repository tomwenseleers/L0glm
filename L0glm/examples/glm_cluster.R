# Native likelihood clustering of response channels
set.seed(12)
Y <- matrix(rpois(120*12, rep(c(2, 12), each=120*6)), 120, 12)
cl <- glm_cluster(Y, k.max=6L, nthreads=2L)
cl$history
cl$clusterids
dim(cl$Y_agg)

# A fixed-variance Gaussian example uses exact Bregman likelihood merges.
Z <- cbind(matrix(rnorm(60*5, -2, .5),60),
           matrix(rnorm(60*5,  2, .5),60))
normal <- glm_cluster(Z, family=gaussian(), dispersion=.25, k.max=5L)
normal$clusterids

\donttest{
data(GCMS_data2)
Y <- as.matrix(GCMS_data2)
# Ordered scans allow coarsening for the clustering stage.
ward <- glm_cluster(Y, method="poissondist", linkage="ward",
                    size_fac=colSums(Y), rowbins=5L, k.max=32L)
bregman <- glm_cluster(Y, method="poissondist", linkage="bregman",
                       size_fac=colSums(Y), rowbins=5L, k.max=32L)
matplot(ward$Y_agg^.2,type="l",lty=1,xlab="Scan",ylab="Summed intensity ^ 0.2")
}
