\donttest{
data(GCMS_data2)                        # GC/MS beer-volatile window
Y <- GCMS_data2                         # 1764 scans x 518 m/z traces

##--- 1. Ward-Bregman clustering with Poisson size-factor distance ---
## We bin every 5 successive scans to speed up the distance calculation
out_pois <- glm_cluster(Y,
                        family   = poisson(identity),
                        method   = "poissondist",
                        rowbins  = 5)

dim(out_pois$Y_agg)                    # 1764 scans x 4 trace-clusters
matplot(out_pois$Y_agg ^ 0.2, type = "l", lty = 1,
        xlab = "Scan index",
        ylab = "Cluster intensity ^ 0.2",
        main = paste("518 m/z traces collapsed into 4 clusters",
                     "Witten Poisson distance, rowbins = 5",
                     sep = "\n"))

##--- 2. Same with symmetrised deviance (Jensen-Shannon) ---
out_js  <- glm_cluster(Y,
                       family   = poisson(identity),
                       method   = "deviance",
                       rowbins  = 5)

dim(out_js$Y_agg)                      # 1764 scans x 5 trace-clusters
matplot(out_js$Y_agg ^ 0.2, type = "l", lty = 1,
        xlab = "Scan index",
        ylab = "Cluster intensity ^ 0.2",
        main = paste("518 m/z traces collapsed into 5 clusters",
                     "Jensen-Shannon deviance, rowbins = 5",
                     sep = "\n"))
}
