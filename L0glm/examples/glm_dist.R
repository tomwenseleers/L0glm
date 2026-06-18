## example: calculating distance among different m/z mass traces in GC/MS datasets
data(GCMS_data2) # section of a GC/MS chromatogram of beer volatiles
Y <- GCMS_data2
dim(Y) # 1764 timepoints x 518 m/z mass traces
# we aggregate together every bin of 10 subsequent timepoints for extra performance:
Y2 <- t(glm_aggregate(t(Y), 
                    family = poisson(identity), 
                    clusterids = rep(1L:ceiling(nrow(Y)/10), each=10)[1:nrow(Y)])$Y_agg)
dim(Y2) # 177 time bins x 518 m/z mass channels
D1 <- glm_dist(Y = Y2, family  = poisson(identity), method="poissondist")
# identical to D1b <- PoiClaClu::PoissonDistance(t(Y2), type = "mle", transform = FALSE)$dd
# but faster
