

# A) Bonferroni-FWER selector
sis_select_fwer <- function(Tj, alpha = 0.05) {
  p <- length(Tj)
  thr <- stats::qchisq(1 - alpha / p, df = 1)
  idx <- which(Tj >= thr)
  list(indices = idx[order(Tj[idx], decreasing = TRUE)],
       threshold = thr, rule = "FWER")
}

# B) BH-FDR selector (with optional cap)
sis_select_fdr <- function(Tj, q = 0.1, cap = Inf) {
  pvals <- stats::pchisq(Tj, df = 1, lower.tail = FALSE)
  pBH   <- stats::p.adjust(pvals, method = "BH")
  idx   <- which(pBH <= q)
  if (length(idx) == 0L) {
    # fallback: take the top one to keep ISIS alive
    o1 <- which.max(Tj)
    idx <- o1
  }
  # order by strength and cap
  o <- idx[order(Tj[idx], decreasing = TRUE)]
  if (is.finite(cap) && length(o) > cap) o <- o[seq_len(cap)]
  list(indices = o, threshold = max(Tj[o]), rule = "BH-FDR", q = q)
}

# C) Permutation/parametric-null max-threshold
sis_select_perm <- function(X, y, family, weights = NULL, offset = NULL,
                            B = 20, alpha = 0.05, seed = 1) {
  set.seed(seed)
  base <- sis_score_stats(X, y, family, weights, offset)
  Tj   <- base$T
  
  # Simulate null responses with the *null* GLM (parametric)
  n <- NROW(X)
  mu0 <- stats::glm.fit(x = matrix(1, n, 1), y = y,
                        weights = weights, offset = offset,
                        family = family, intercept = TRUE)$fitted.values
  
  maxTs <- numeric(B)
  for (b in seq_len(B)) {
    y0 <- switch(family$family,
                 "binomial" = stats::rbinom(n, size = 1, prob = mu0),
                 "poisson"  = stats::rpois(n, lambda = mu0),
                 "gaussian" = stats::rnorm(n, mean = mu0, sd = sqrt(stats::summary(glm(y ~ 1, family = gaussian()))$dispersion)),
                 { # generic fallback via family$rng if available
                   if (!is.null(family$rng)) family$rng(n, mu0) else stop("Add a null generator for this family.")
                 }
    )
    maxTs[b] <- max(sis_score_stats(X, y0, family, weights, offset)$T)
  }
  thr <- stats::quantile(maxTs, probs = 1 - alpha, names = FALSE)
  idx <- which(Tj >= thr)
  list(indices = idx[order(Tj[idx], decreasing = TRUE)],
       threshold = thr, rule = "perm-max", alpha = alpha, T = Tj)
}

# Knee (elbow) selector to decide on best nsis with cap/floor
sis_select_knee <- function(Tj, floor =  min(n, round(n/(log(log(n))*log(p)))), cap = n ) {
  o   <- order(Tj, decreasing = TRUE)
  srt <- Tj[o]
  # work on semilog-scale for stability
  z   <- log1p(srt)
  # discrete second derivative
  d2  <- c(NA, diff(diff(z)), NA)
  # plot(d2[1:10000], type="h")
  k   <- which.max(pmax(d2[!is.na(d2)][-c(1:floor)], 0))+floor
  k   <- max(floor, min(k, length(o)))
  if (is.finite(cap)) k <- min(k, cap)
  list(indices = o[seq_len(k)], threshold = srt[k], rule = "knee")
}



# test with standard normal dense design matrix with p=1E6, n=500, k=60
n <- 500; p <- 1e6; k <- 60
set.seed(1)
X <- matrix(rnorm(n * p), n, p)
beta_true <- numeric(p); beta_true[sample(p, k)] <- runif(k, -0.8, 0.8)
beta_true_set <- which(beta_true!=0)
y_true <- as.numeric(X %*% beta_true) 
y <- y_true + rnorm(n)
family <- gaussian(identity)

IC <- "gic"
lambda_preset = preset.lambda(IC = IC, # options aic, bic, gic, ebic, mbic
                              p = ncol(X), # note that we use original problem size
                              n = nrow(X))  
nonnegative <- FALSE

system.time(sis <- sis_glm(X, y, family, 
               weights = NULL, 
               offset = NULL, 
               nsis = n, # min(n, round(n/(log(log(n))*log(p)))), # max support size in abess = default
               return_stats = TRUE)) # 2s
Tj=sis$stats$stat
Tj_signed=Tj*sis$stats$sign
pvals <- stats::pchisq(Tj, df = 1, lower.tail = FALSE)
pBH   <- stats::p.adjust(pvals, method = "BH")
pBonf <- stats::p.adjust(pvals, method = "bonferroni")
sel_fwer = sis_select_fwer(Tj, alpha=0.1)
length(sel_fwer$indices) # 111522
sel_fdr = sis_select_fdr(Tj, q = 0.01, cap = Inf)
length(sel_fdr$indices) # 388824
sel_knee = sis_select_knee(Tj, cap = n, floor = min(n, round(n/(log(log(n))*log(p)))) )
length(sel_knee$indices) # 35 -> works best??
length(intersect(sel_knee$indices, beta_true_set)) # TPs: 9
pBHcutoff <- 0.01
sel_knee_withFDR <- sel_knee$indices[pBH[sel_knee$indices]<pBHcutoff]
length(sel_knee_withFDR) # 35 -> does not help but also does not hurt
length(intersect(sel_knee_withFDR, beta_true_set)) # TPs: 9
pBonfcutoff <- 0.1
sel_knee_withBonf <- sel_knee$indices[pBonf[sel_knee$indices]<pBonfcutoff]
length(sel_knee_withBonf) # 35 -> does not help but also does not hurt
length(intersect(sel_knee_withBonf, beta_true_set)) # TPs: 9

# automatic way to decide on nsis?
floor = min(n, round(n/(log(log(n))*log(p)))) 
floor # 20
ceiling = n
ceiling # 500
best = which.min(abs(length(sel_knee$indices)-c(floor, ceiling)))
# SIS = (best==1) # apply SIS or not? -> here TRUE
# SIS
nsis = c(floor, ceiling)[best]
nsis # 20

# iterating SIS & L0glm:
beta_start <- rep(0, ncol(X))
S <- Snew <- Sknee <- NULL
# S <- which(beta_start!=0)
# priormean <- rep(0, p)
disp <- NULL
timing_SISL0glm <- system.time(for (iter in 1:1) {
  sis <- sis_glm(X, y, family, 
                 weights = NULL, 
                 offset = X %*% beta_start,
                 nsis = n, # min(n, round(n/(log(log(n))*log(p)))), # max support size in abess = default
                 return_stats = TRUE)
  Tj=sis$stats$stat
  Tj_signed=Tj*sis$stats$sign
  pvals <- stats::pchisq(Tj, df = 1, lower.tail = FALSE)
  pBH   <- stats::p.adjust(pvals, method = "BH")
  pBonf <- stats::p.adjust(pvals, method = "bonferroni")
  
  sel_knee = sis_select_knee(Tj, 
                             cap = n, 
                             floor = min(n, round(n/(log(log(n))*log(p)))) )
  Sknee <- sel_knee$indices
  print(paste("Iteration", iter))
  print(paste("Size of knee selected SIS set:", length(Sknee)))

  # pBonfcutoff <- 0.1
  # Snew <- Sknee[pBonf[Sknee]<pBonfcutoff]
  pBHcutoff <- 0.01
  Snew <- Sknee[pBH[Sknee]<pBHcutoff]
  
  print(paste("Size of knee selected SIS set with pBH<pBHcutoff:", length(Sknee[pBH[Sknee]<0.001])))
  print(paste("Size of knee selected SIS set with pBonf<pBonfcutoff:", length(Sknee[pBonf[Sknee]<0.001])))
  
  # Snew <- sis$indices
  # TODO fix Warning message:
  # In Matrix::Matrix(offset, n, 1, sparse = TRUE) :
  #  'nrow', 'ncol', 'byrow' disregarded for [mM]atrix 'data'
  S   <- union(S, Snew)
  print(paste("Size of new SIS set:", length(Snew)))
  print(paste("Number of true positives in new SIS set:", length(intersect(Snew, beta_true_set))))
  print(paste("Total size of current SIS set:", length(S)))
  print(paste("Number of true positives in current total SIS set:", length(intersect(S, beta_true_set))))
#  for (j in 1:1) { # iterated twice to remove some false positives
    L0glmfit_updated <- L0glm.fit(X = X[,S,drop=F], 
                                  y = y, 
                                  family = family, 
                                  lambda = lambda_preset,  
                                  nonnegative = nonnegative,
                                  start = rep(1, length(S)), # if (iter==1) rep(1, length(S)) else beta_start[S], # better solution than NULL 
                                  # prior.mean = priormean[S],
                                  # dispersion = disp,
                                  solver = "eigen", # , # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                                  method = 1L, # Cholesky
                                  thresh = 1E-5,
                                  gamma = 2.1 # taking gamma slightly higher than 2 for p>>n case
    )  # 
    sum(L0glmfit_updated$coefficients!=0) # 8
    beta_start = rep(0, ncol(X))
    beta_start[S] = L0glmfit_updated$coefficients
    S <- which(beta_start!=0)
    # largestcoefset <- which(abs(beta_start)>median(abs(beta_start[S]))) 
    # priormean[S] <- beta_start[largestcoefset]
    # disp <- L0glmfit_updated$dispersion
#  }
  print(paste("Number of nonzero coefficients:", sum(beta_start!=0)))
  print(paste("Number of true positives:", length(intersect(S, beta_true_set))))
  
})[[3]] 
timing_SISL0glm
# 33s - we should get 31 TPs
# gic & hbic return 31 TPs & 1 FP, ebic returns 29 TPs & 1 FP, while mbic returns only 16 TPs
# TODO check if it helps to set non-top n coefs to zero in L0glm.fit.internal

# system.time(L0glmfit_updated <- L0glm.fit(X = X[,S,drop=F], 
#                                         y = y, 
#                                         family = family, 
#                                         lambda = lambda_preset,  
#                                         nonnegative = nonnegative,
#                                         start = beta_start[S],  
#                                         # prior.mean = beta_start, # [S],
#                                         solver = "eigen", # , # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
#                                         method = 1L, # Cholesky
#                                         thresh = 1E-5,
#                                         gamma = 2.1
# ))  # 0.09s

beta_final = beta_start # rep(0, ncol(X))
# beta_final[S] = L0glmfit_updated$coefficients

# beta_final[S] = L0glmfit_updated$coefficients[L0glmfit_updated$coefficients!=0]

sum(beta_final!=0) # 32

beta_SISL0glm = as.numeric(beta_final) # beta_start
table(beta_SISL0glm!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))

fitted = X %*% beta_final
summary(lm(y_true ~ 0 + fitted))$r.squared # population R2=0.7547788
plot(x = fitted, y = y_true, pch = 16, col = "steelblue")




# test with sparse Toeplitz design
n <- p <- 10000 
npeaks <- 1000
sim.pois <- simul(method = 2, n = n, p = p, s = npeaks/p, sparse = T)
beta_true <- as.numeric(sim.pois$beta_true)
beta_true_set <- which(beta_true!=0)
X <- sim.pois$X
y <- sim.pois$y
family <- poisson(identity)
nonnegative <- TRUE

system.time(sis <- sis_glm(X, y, family, 
                           weights = NULL, 
                           offset = NULL, 
                           nsis = n, # min(n, round(n/(log(log(n))*log(p)))), # max support size in abess = default
                           return_stats = TRUE)) # 0.02s
Tj=sis$stats$stat
Tj_signed=Tj*sis$stats$sign
pvals <- stats::pchisq(Tj, df = 1, lower.tail = FALSE)
pBH   <- stats::p.adjust(pvals, method = "BH")
pBonf <- stats::p.adjust(pvals, method = "bonferroni")
sel_fwer = sis_select_fwer(Tj, alpha=0.1)
length(sel_fwer$indices) # 9276
sel_fdr = sis_select_fdr(Tj, q = 0.01, cap = Inf)
length(sel_fdr$indices) # 9336
sel_knee = sis_select_knee(Tj, cap = n, floor = min(n, round(n/(log(log(n))*log(p)))) )
length(sel_knee$indices) # 9364 -> works best??
length(intersect(sel_knee$indices, beta_true_set)) # TPs: 1000
pBHcutoff <- 0.01
sel_knee_withFDR <- sel_knee$indices[pBH[sel_knee$indices]<pBHcutoff]
length(sel_knee_withFDR) # 9336 -> does not help but also does not hurt
length(intersect(sel_knee_withFDR, beta_true_set)) # TPs: 1000
pBonfcutoff <- 0.1
sel_knee_withBonf <- sel_knee$indices[pBonf[sel_knee$indices]<pBonfcutoff]
length(sel_knee_withBonf) # 9276 -> does not help but also does not hurt
length(intersect(sel_knee_withBonf, beta_true_set)) # TPs: 1000

# automatic way to decide whether or not to apply (I)SIS?
floor = min(n, round(n/(log(log(n))*log(p)))) 
floor # 489
ceiling = n
ceiling # 10000
best = which.min(abs(length(sel_knee$indices)-c(floor, ceiling)))
SIS = (best==1) # apply SIS or not? -> here FALSE
# ALT if length(sel_knee$indices) is within 10% of ceiling use nsis=ceiling
nsis = c(floor, ceiling)[best]
nsis # 10000

# iterating SIS & L0glm:
beta_start <- rep(0, ncol(X))
S <- Snew <- Sknee <- NULL
# S <- which(beta_start!=0)
# priormean <- rep(0, p)
disp <- NULL
timing_SISL0glm <- system.time(for (iter in 1:10) {
  system.time(sis <- sis_glm(X, y, family, 
                 weights = NULL, 
                 offset = X %*% beta_start,
                 nsis = min(n, round(n/(log(log(n))*log(p)))), # max support size in abess = default
                 return_stats = TRUE))
  Tj=sis$stats$stat
  if (nonnegative) { Tj[sis$stats$sign==-1] = 0 }
  Tj_signed=Tj*sis$stats$sign
  pvals <- stats::pchisq(Tj, df = 1, lower.tail = FALSE)
  pBH   <- stats::p.adjust(pvals, method = "BH")
  pBonf <- stats::p.adjust(pvals, method = "bonferroni")
  
  sel_knee = sis_select_knee(Tj, 
                             cap = n, 
                             floor = min(n, round(n/(log(log(n))*log(p)))) )
  Sknee <- sel_knee$indices
  print(paste("Iteration", iter))
  print(paste("Size of knee selected SIS set:", length(Sknee)))
  
  # pBonfcutoff <- 0.1
  # Snew <- Sknee[pBonf[Sknee]<pBonfcutoff]
  pBHcutoff <- 0.01
  Snew <- Sknee[pBH[Sknee]<pBHcutoff]
  
  print(paste("Size of knee selected SIS set with pBH<pBHcutoff:", length(Sknee[pBH[Sknee]<0.001])))
  print(paste("Size of knee selected SIS set with pBonf<pBonfcutoff:", length(Sknee[pBonf[Sknee]<0.001])))
  
  # Snew <- sis$indices
  # TODO fix Warning message:
  # In Matrix::Matrix(offset, n, 1, sparse = TRUE) :
  #  'nrow', 'ncol', 'byrow' disregarded for [mM]atrix 'data'
  S   <- union(S, Snew)
  # S <- beta_true_set
  print(paste("Size of new SIS set:", length(Snew)))
  print(paste("Number of true positives in new SIS set:", length(intersect(Snew, beta_true_set))))
  print(paste("Total size of current SIS set:", length(S)))
  print(paste("Number of true positives in current total SIS set:", length(intersect(S, beta_true_set))))
  # for (j in 1:1) { # iterated twice to remove some false positives
    L0glmfit_updated <- L0glm.fit(X = X[,S,drop=F], 
                                  y = y, 
                                  family = family, 
                                  lambda = IC, # lambda_preset,  
                                  nonnegative = nonnegative,
                                  start = NULL, # rep(1, length(S)), #  sqrt(Tj_signed[S]), # rep(1, length(S)), #  # # if (iter==1) rep(1, length(S)) else beta_start[S], # better solution than NULL 
                                  maxit = 100,
                                  # prior.mean = priormean[S],
                                  # dispersion = disp,
                                  solver = "eigen", # , # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                                  method = 1L, # Cholesky
                                  thresh = 1E-5,
                                  gamma = 2 # .1 # taking gamma slightly higher than 2 for p>>n case
    )  # 
    sum(L0glmfit_updated$coefficients!=0) # 8
    beta_start = rep(0, ncol(X))
    beta_start[S] = L0glmfit_updated$coefficients
    S <- which(beta_start!=0)
    # largestcoefset <- which(abs(beta_start)>median(abs(beta_start[S]))) 
    # priormean[S] <- beta_start[largestcoefset]
    # disp <- L0glmfit_updated$dispersion
#  }
  print(paste("Number of nonzero coefficients:", sum(beta_start!=0)))
  print(paste("Number of true positives:", length(intersect(S, beta_true_set))))
  
})[[3]] 
timing_SISL0glm
# 6.3s
# gic & hbic return 31 TPs & 1 FP, ebic returns 29 TPs & 1 FP, while mbic returns only 16 TPs
# TODO check if it helps to set non-top n coefs to zero in L0glm.fit.internal

# system.time(L0glmfit_updated <- L0glm.fit(X = X[,S,drop=F], 
#                                         y = y, 
#                                         family = family, 
#                                         lambda = lambda_preset,  
#                                         nonnegative = nonnegative,
#                                         start = beta_start[S],  
#                                         # prior.mean = beta_start, # [S],
#                                         solver = "eigen", # , # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
#                                         method = 1L, # Cholesky
#                                         thresh = 1E-5,
#                                         gamma = 2.1
# ))  # 0.09s

# beta_final = rep(0, ncol(X))
# beta_final[S] = L0glmfit_updated$coefficients

# beta_final[S] = L0glmfit_updated$coefficients[L0glmfit_updated$coefficients!=0]

beta_final = beta_start
sum(beta_final!=0) # 1005

beta_SISL0glm = as.numeric(beta_final) # beta_start
table(beta_SISL0glm!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))


L0glm.chol.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                             family = poisson(identity), 
                             lambda = IC, # or tune_zoomed$lambda_best
                             nonnegative = TRUE,
                             relaxed = FALSE,
                             thresh = 1E-5,
                             # verbose = TRUE,
                             solver = "eigen",
                             method = 1L)
table(as.numeric(L0glm.chol.pois$coefficients)!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))


rbind(data.frame(method="L0glm_ISIS", do.call(cbind, get.metrics(beta = Matrix(beta_final, sparse=T),  # some metrics
                                                                     X = X, y = y, 
                                                                     family = poisson(identity),
                                                                     beta_true = sim.pois$beta_true) )),
      data.frame(method="L0glm_cholesky", do.call(cbind, get.metrics(beta = L0glm.chol.pois$coefficients,  # some metrics
                                                               X = X, y = y, 
                                                               family = poisson(identity),
                                                               beta_true = sim.pois$beta_true) )))
#           method     p     n m deviance deviance_true    loglik    MD_coefs
# 1     L0glm_ISIS 10000 10000 1 9120.703      2748.035 -40921.83 -0.02052188
# 2 L0glm_cholesky 10000 10000 1 9348.514      2987.828 -41035.74 -0.02862166
#   MAD_coefs sMAPE_coefs bias_coefs  R2_coefs MD_selcoefs MAD_selcoefs
# 1  74.26700  0.07082828  0.9952626 0.9192988    254.3827     484.3882
# 2  80.26889  0.07160071  0.9928768 0.9093261    262.8916     529.9963
#   sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1      0.3863510     0.9952626   0.9465471              1000
# 2      0.3992165     0.9928768   0.9399447              1000
#   support_size_estimated estimated_over_true_support_size  TP  FP   TN  FN
# 1                   1005                            1.005 840 165 8835 160
# 2                   1012                            1.012 844 168 8832 156
# precision recall specificity accuracy  F1_score      gic     ebic      bic
# 1 0.8358209  0.840   0.9816667   0.9675 0.8379052 102395.9 100356.4 91100.05
# 2 0.8339921  0.844   0.9813333   0.9676 0.8389662 102766.8 100713.2 91392.34
# aic     mbic
# 1 83853.66 106826.4
# 2 84095.47 107228.2

coefs=beta_final #  # as.numeric(L0glm.chol.pois$coefficients) # 
plot(x=coefs[beta_true!=0&coefs!=0], 
     y=beta_true[beta_true!=0&coefs!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main=paste0("ISIS-L0glm best subset fit (GIC tuned)")) # true positives
abline(a=0, b=1, col="lightgrey")
points(x=coefs[beta_true!=0&coefs==0], 
       y=beta_true[beta_true!=0&coefs==0], 
       pch=16, col="purple") # false negatives
points(x=coefs[beta_true==0&coefs!=0], 
       y=beta_true[beta_true==0&coefs!=0], 
       pch=16, col="red") # false positives
points(x=coefs[beta_true==0&coefs==0][[1]], 
       y=beta_true[beta_true==0&coefs==0][[1]], 
       pch=16, col="green") # true negatives
method="ISIS_L0glm"



