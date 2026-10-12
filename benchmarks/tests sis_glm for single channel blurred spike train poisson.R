n <- p <- 10000 
npeaks <- 1000
sim.pois <- simul(method = 2, n = n, p = p, s = npeaks/p, sparse = T)
beta_true <- as.numeric(sim.pois$beta_true)
beta_true_set <- which(beta_true!=0)
X <- sim.pois$X
y <- sim.pois$y
weights <- 1/(y+0.1)
nonnegative <- TRUE
family <- poisson(identity)
nsis <- round(n / log(n))
system.time(SIS_glm <- sis_glm(X = X, y = y, weights = weights,
        family = family,
        offset = NULL,
        nsis = nsis, return_stats = FALSE)) # 0.01s
length(SIS_glm$indices) # 1086
length(intersect(SIS_glm$indices, which(beta_true!=0))) # 265

IC <- "gic"
lambda_preset = preset.lambda(IC = IC, # options aic, bic, gic, mebic, ebic, mbic
                              p = ncol(X), # note that we use original problem size
                              n = nrow(X))  
# iterating SIS & L0glm:
beta_start <- rep(0, ncol(X))
S <- NULL
timing_SISL0glm <- system.time(for (iter in 1:10) {
  sis <- sis_glm(X = X, y = y, family = family, 
                 # weights = weights, 
                 offset = X %*% beta_start,
                 nsis = n, # n-1, i.e. no prescreening works best here; round(n/log(n)), # min(n, round(n/(log(log(n))*log(p)))), # max support size in abess = default
                 return_stats = FALSE)
  S   <- union(S, sis$indices)
  print(paste("Iteration", iter))
  print(paste("Size of SIS set:", length(S)))
  for (j in 1:2) { # iterated twice to remove some false positives
    L0glmfit_updated <- L0glm.fit(X = X[,S,drop=F], 
                                  y = y, 
                                  family = family,
                                  weights = NULL,
                                  lambda = lambda_preset,  
                                  nonnegative = nonnegative,
                                  start = rep(1, length(S)), # if (iter==1) rep(1, length(S)) else beta_start[S], # better solution than NULL 
                                  warm.start = 1,
                                  solver = "eigen", # , # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                                  method = 1L, # Cholesky
                                  gamma = 2 # taking gamma slightly higher than 2 for p>>n case
    )  # 
    beta_start = rep(0, ncol(X))
    beta_start[S] = L0glmfit_updated$coefficients
    S <- which(beta_start!=0)
  }
  print(paste("Number of nonzero coefficients:", sum(beta_start!=0)))
  print(paste("Number of true positives:", length(intersect(S, beta_true_set))))
  
})[[3]] 
timing_SISL0glm
# 28s

system.time(L0glmfit_final <- L0glm.fit(X = X[,S,drop=F], 
                                        y = y, 
                                        family = family, 
                                        lambda = lambda_preset,  
                                        nonnegative = nonnegative,
                                        start = beta_start[S],  
                                        warm.start = 0.99,
                                        solver = "eigen", # , # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                                        method = 1L, # Cholesky
                                        gamma = 2
))  # 0.09s
beta_final = rep(0, ncol(X))
beta_final[S] = L0glmfit_final$coefficients
# beta_final = L0glmfit_final$coefficients
sum(beta_final!=0) # 725

beta_SISL0glm = as.numeric(beta_final) # beta_start

plot(x=beta_SISL0glm[beta_true!=0&beta_SISL0glm!=0], 
     y=beta_true[beta_true!=0&beta_SISL0glm!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main=paste0("ISIS-L0glm best subset fit (GIC tuned)")) # true positives
abline(a=0, b=1, col="lightgrey")
points(x=beta_SISL0glm[beta_true!=0&beta_SISL0glm==0], 
       y=beta_true[beta_true!=0&beta_SISL0glm==0], 
       pch=16, col="purple") # false negatives
points(x=beta_SISL0glm[beta_true==0&beta_SISL0glm!=0], 
       y=beta_true[beta_true==0&beta_SISL0glm!=0], 
       pch=16, col="red") # false positives
points(x=beta_SISL0glm[beta_true==0&beta_SISL0glm==0][[1]], 
       y=beta_true[beta_true==0&beta_SISL0glm==0][[1]], 
       pch=16, col="green") # true negatives
method="ISIS_L0glm"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

table(beta_SISL0glm!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))

cor(beta_SISL0glm, beta_true)     # 0.9679615
sum(abs(beta_SISL0glm-beta_true)) # 5.539483

fitted = X %*% beta_SISL0glm
cor(y_true, fitted)^2 # population R2=0.9480623

# check metrics & quality of ISIS-L0glm solution
method = "ISIS_L0glm"
system.time(metrics_ISISL0glm <- data.frame(method=method,
                                            timing_s=timing_SISL0glm,
                                            do.call(cbind, get.metrics(beta = beta_SISL0glm,  # some metrics
                                                                       X = X, 
                                                                       y = y,
                                                                       family = gaussian(identity), 
                                                                       beta_true = beta_true)))) 
metrics_ISISL0glm


